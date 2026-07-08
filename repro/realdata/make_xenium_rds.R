#!/usr/bin/env Rscript
# Reshape 10x Xenium FFPE Human Breast Cancer Rep1 into the run_singlecell.R contract:
#   list(counts = cells x genes integer, coords = cells x 2 [x,y], disease = 0/1 per cell)
# Prereq: bash repro/realdata/fetch_xenium_breast.sh   (downloads matrix + cells.csv + analysis)
#
# Env overrides: XENIUM_DIR (default /tmp/xenium_breast), OUT_RDS, LABEL_CSV (authoritative labels),
#   HVG (keep only the top-K highly-variable genes; 0 = keep all). fit_pgcavi is O(N*P^2)/iter, so
#   for the 313-gene Xenium panel set e.g. HVG=50 (as run_her2st.R does) to keep the fit fast.
suppressMessages(library(Matrix))
DEST <- Sys.getenv("XENIUM_DIR", "/tmp/xenium_breast")
OUT  <- Sys.getenv("OUT_RDS", file.path(DEST, "xenium_breast_rep1.rds"))
mdir <- file.path(DEST, "cell_feature_matrix")
ccsv <- file.path(DEST, "Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.csv.gz")
stopifnot(dir.exists(mdir), file.exists(ccsv))

## --- counts: 10x MTX (genes x cells) -> keep Gene Expression rows only ---
mat  <- readMM(gzfile(file.path(mdir, "matrix.mtx.gz")))
feat <- read.delim(gzfile(file.path(mdir, "features.tsv.gz")), header = FALSE)  # V1 id, V2 name, V3 type
bc   <- read.delim(gzfile(file.path(mdir, "barcodes.tsv.gz")), header = FALSE)$V1
rownames(mat) <- make.unique(feat$V2); colnames(mat) <- as.character(bc)
mat <- mat[feat$V3 == "Gene Expression", , drop = FALSE]           # drop control/blank probes
cat(sprintf("counts: %d genes x %d cells (Gene Expression only)\n", nrow(mat), ncol(mat)))

## --- coords: per-cell centroids, aligned to matrix cell order ---
cells <- read.csv(gzfile(ccsv))
rownames(cells) <- as.character(cells$cell_id)
cells <- cells[colnames(mat), ]
stopifnot(!any(is.na(cells$x_centroid)))
coords <- as.matrix(cells[, c("x_centroid", "y_centroid")])

## --- disease label ---
# Authoritative route: supervised cell-type annotation (Janesick 2023: DCIS/invasive tumor -> 1),
# supplied as a 2-column CSV (cell_id, label) via LABEL_CSV. Otherwise a self-contained
# tumor-epithelial MARKER-SCORE proxy so the pipeline runs out-of-the-box.
LABEL_CSV <- Sys.getenv("LABEL_CSV", "")
if (nzchar(LABEL_CSV) && file.exists(LABEL_CSV)) {
  lab <- read.csv(LABEL_CSV); rownames(lab) <- as.character(lab[[1]])
  ct  <- as.character(lab[colnames(mat), 2])                        # cell-type per cell (aligned)
  # tumor cell types -> disease=1 (Janesick supervised: DCIS #1/#2, Invasive Tumor, Prolif Inv Tumor)
  disease <- as.integer(grepl("DCIS|Invasive|Tumou?r|Carcinoma", ct, ignore.case = TRUE) |
                          ct %in% c("1", "tumor", "TRUE"))
  disease[is.na(disease)] <- 0L
  cat(sprintf("disease: supervised annotation %s -> %.0f%% disease=1\n", LABEL_CSV, 100*mean(disease)))
} else {
  markers <- intersect(c("ERBB2","EPCAM","KRT8","KRT18","FASN","ELF3","CDH1","TACSTD2"), rownames(mat))
  libsize <- pmax(Matrix::colSums(mat), 1)
  cpm     <- mat %*% Matrix::Diagonal(x = 1/libsize) * 100           # per-cell normalized
  score   <- Matrix::colSums(cpm[markers, , drop = FALSE])
  disease <- as.integer(score > median(score[score > 0]))            # top-half epithelial -> 1
  cat(sprintf("disease: MARKER-SCORE proxy on {%s} -> %.0f%% disease=1\n",
              paste(markers, collapse = ","), 100*mean(disease)))
  cat("  NOTE: heuristic proxy for scale/pipeline demo. For the headline biological result,\n")
  cat("        pass LABEL_CSV with the Janesick supervised tumor annotation.\n")
}

## --- optional HVG reduction (after disease, so tumor markers were available above) ---
# fit_pgcavi is O(N*P^2)/iter; the 313-gene panel is slow at N=1.7e5. Keep top-K HVGs.
HVG <- as.integer(Sys.getenv("HVG", "0"))
if (!is.na(HVG) && HVG > 0 && HVG < nrow(mat)) {
  lib <- pmax(Matrix::colSums(mat), 1)
  lc  <- log1p(mat %*% Matrix::Diagonal(x = 1e4 / lib))     # log-CPM, stays sparse
  gv  <- Matrix::rowMeans(lc^2) - Matrix::rowMeans(lc)^2     # per-gene variance
  keep <- order(gv, decreasing = TRUE)[1:HVG]
  mat <- mat[keep, , drop = FALSE]
  cat(sprintf("HVG: kept top %d/%d genes by log-CPM variance\n", HVG, length(gv)))
}

counts <- t(as.matrix(mat))                                          # cells x genes integer
storage.mode(counts) <- "integer"
saveRDS(list(counts = counts, coords = coords, disease = disease), OUT)
cat(sprintf("wrote %s  (%d cells x %d genes, %.0f%% disease)\n",
            OUT, nrow(counts), ncol(counts), 100*mean(disease)))
cat(sprintf("Next: DATA_RDS=%s Rscript repro/realdata/run_singlecell.R\n", OUT))
