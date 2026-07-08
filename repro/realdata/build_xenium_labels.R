#!/usr/bin/env Rscript
# Build an INDEPENDENT (cluster-based) tumor label for Xenium Breast Rep1 -> a LABEL_CSV that
# make_xenium_rds.R consumes. Cell identity comes from 10x's UNSUPERVISED graph clustering over the
# FULL 313-gene transcriptome; each cluster is assigned to a lineage (epithelial / immune / stromal)
# by canonical marker signatures, and epithelial-dominant clusters -> disease=1 (tumor region proxy).
# Far less circular than a per-cell marker threshold: the label is a whole-transcriptome cluster
# identity, annotated with a few markers, not a threshold on the covariates themselves.
#
# Prereq: bash repro/realdata/fetch_xenium_breast.sh   (needs cell_feature_matrix/ + analysis/)
suppressMessages(library(Matrix))
DEST <- Sys.getenv("XENIUM_DIR", "/tmp/xenium_breast")
OUT  <- Sys.getenv("LABEL_OUT", file.path(DEST, "xenium_breast_rep1_labels.csv"))
clf  <- file.path(DEST, "analysis/clustering/gene_expression_graphclust/clusters.csv")
mdir <- file.path(DEST, "cell_feature_matrix")
stopifnot(file.exists(clf), dir.exists(mdir))

## counts (genes x cells), Gene Expression only
mat  <- readMM(gzfile(file.path(mdir, "matrix.mtx.gz")))
feat <- read.delim(gzfile(file.path(mdir, "features.tsv.gz")), header = FALSE)
bc   <- read.delim(gzfile(file.path(mdir, "barcodes.tsv.gz")), header = FALSE)$V1
rownames(mat) <- make.unique(feat$V2); colnames(mat) <- as.character(bc)
mat  <- mat[feat$V3 == "Gene Expression", , drop = FALSE]
lc   <- log1p(mat %*% Matrix::Diagonal(x = 1e4 / pmax(Matrix::colSums(mat), 1)))   # log-CPM

## graph clusters, aligned to matrix cell order
cls <- read.csv(clf); names(cls) <- c("cell_id", "cluster")
cls <- cls[match(colnames(mat), as.character(cls$cell_id)), ]
cl  <- factor(cls$cluster)

## canonical lineage signatures (restricted to genes present in the panel)
sig <- lapply(list(
  epithelial = c("EPCAM","KRT8","KRT7","ELF3","CDH1","FASN","TACSTD2","ERBB2","FOXA1","GATA3","ANKRD30A","CEACAM6","ESR1"),
  immune     = c("PTPRC","CD3D","CD3E","CD8A","CD68","CD14","LYZ","MS4A1","NKG7","ITGAX"),
  stromal    = c("PECAM1","VWF","PDGFRB","PDGFRA","LUM","ACTA2","MYH11")
), intersect, y = rownames(lc))

## per-cell lineage score -> per-cluster mean -> z across clusters -> dominant lineage
cell_score <- sapply(sig, function(g) Matrix::colMeans(lc[g, , drop = FALSE]))
clmean <- apply(cell_score, 2, function(s) tapply(s, cl, mean))
clz    <- scale(clmean)
lineage <- colnames(clz)[max.col(clz)]; names(lineage) <- rownames(clz)
# tumor = epithelial-dominant AND above-average epithelial (z>0), so low-signal/ambiguous
# clusters whose "epithelial" is merely the least-negative lineage are NOT called tumor.
tumor_clusters <- names(lineage)[lineage == "epithelial" & clz[, "epithelial"] > 0]

disease <- as.integer(as.character(cl) %in% tumor_clusters)
write.csv(data.frame(cell_id = colnames(mat), label = disease), OUT, row.names = FALSE)

cat(sprintf("graph clusters: %d | epithelial(tumor): {%s}\n", nlevels(cl), paste(tumor_clusters, collapse=",")))
cat(sprintf("cells disease=1: %d/%d (%.0f%%)\n", sum(disease), length(disease), 100*mean(disease)))
cat("wrote", OUT, "\n\n")
print(data.frame(cluster = rownames(clz), lineage = lineage,
                 epi = round(clz[,"epithelial"],2), imm = round(clz[,"immune"],2), str = round(clz[,"stromal"],2),
                 n = as.integer(table(cl)[rownames(clz)])), row.names = FALSE)
