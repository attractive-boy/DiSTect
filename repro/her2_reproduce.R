# Reproduce Fig 2 pipeline of DiSTect (btaf530) on REAL HER2+ data (her2st, Andersson 2021)
# Scaled reproduction: 4 patients (A-D), one pathologist-labeled section each,
# top-HVG subset, joint multi-slice model via ADVI.
# (Paper's full run = 8 patients x multiple slices x all genes = ~25 h; here minutes.)
suppressMessages({ library(rstan); library(ggplot2); library(dplyr); library(GGally); library(network) })
options(mc.cores = parallel::detectCores()); set.seed(128); seed <- 128

root  <- "/Users/licheng/Documents/DiSTect"
her2  <- "/tmp/her2st/data"
outdir<- file.path(root, "repro"); dir.create(outdir, showWarnings = FALSE)
for (f in c("dsgd.R","prediction.R","coef_plot.R","network_plot.R")) source(file.path(root,"R",f))

SECTIONS <- c(A="A1", B="B1", C="C1", D="D1")   # one labeled section per patient
N_HVG    <- 40
sep <- function(t) cat("\n\n========== ", t, " ==========\n")

## ---- loader: counts + pathologist labels -> per-spot (counts, coords, disease) ----
load_section <- function(sec) {
  cnt <- read.delim(gzfile(file.path(her2,"ST-cnts", paste0(sec,".tsv.gz"))),
                     row.names = 1, check.names = FALSE)          # spots x genes
  lab <- read.delim(file.path(her2,"ST-pat","lbl",
                              paste0(sec,"_labeled_coordinates.tsv")), check.names = FALSE)
  lab$key <- paste0(round(lab$x), "x", round(lab$y))               # array-coord key
  lab <- lab[!duplicated(lab$key), ]
  rownames(lab) <- lab$key
  common <- intersect(rownames(cnt), lab$key)                     # spots with a label
  cnt <- cnt[common, , drop = FALSE]
  lab <- lab[common, , drop = FALSE]
  keep <- !grepl("undetermined", lab$label, ignore.case = TRUE)   # drop undetermined
  cnt <- cnt[keep, ]; lab <- lab[keep, ]
  disease <- as.integer(grepl("cancer", lab$label, ignore.case = TRUE))  # cancer vs non-cancer
  coords  <- do.call(rbind, strsplit(rownames(cnt), "x"))
  list(counts  = as.matrix(cnt),
       coords  = data.frame(x = as.numeric(coords[,1]), y = as.numeric(coords[,2])),
       disease = disease)
}

sep("LOAD 4 patient sections + pathologist annotations")
secs <- lapply(SECTIONS, load_section)
for (i in seq_along(secs))
  cat(sprintf("  %s (%s): %d spots, %d cancer (%.0f%%)\n",
              names(SECTIONS)[i], SECTIONS[i], nrow(secs[[i]]$counts),
              sum(secs[[i]]$disease), 100*mean(secs[[i]]$disease)))

## ---- gene filtering (paper): drop total count <300 per section, common subset ----
gene_ok <- lapply(secs, function(s) colnames(s$counts)[colSums(s$counts) >= 300])
common_genes <- Reduce(intersect, gene_ok)
cat(sprintf("\nGenes passing total-count>=300 in ALL sections (common subset): %d\n",
            length(common_genes)))

## ---- combine, log(X+1), pick top-HVG, standardize ----
Xlist <- lapply(secs, function(s) log1p(s$counts[, common_genes, drop = FALSE]))
Xall  <- do.call(rbind, Xlist)
y     <- unlist(lapply(secs, `[[`, "disease"))
label <- rep(seq_along(secs), sapply(secs, function(s) nrow(s$counts)))
coord <- do.call(rbind, lapply(secs, `[[`, "coords"))
hvg   <- names(sort(apply(Xall, 2, var), decreasing = TRUE))[1:N_HVG]
Xh    <- scale(Xall[, hvg])                                    # z-score genes
matrix_x <- cbind(as.data.frame(Xh), x = coord$x, y = coord$y)
cat(sprintf("Design: %d spots x %d HVGs (+2 coords), %d slices/patients. Cancer rate %.0f%%\n",
            nrow(matrix_x), N_HVG, length(unique(label)), 100*mean(y)))

## ---- FIT joint multi-slice model (Fig 2a) ----
sep("FIT joint multi-slice DiSTect (ADVI)")
t_fit <- system.time(fit <- dsgd(list_y = y, matrix_x = as.matrix(matrix_x), label_list = label))
cat(sprintf("Fit time: %.1f s\n", t_fit["elapsed"]))
s <- rstan::summary(fit)$summary
b <- s[grep("^beta\\[", rownames(s)), ]
rank <- data.frame(gene = hvg, coef = b[,"mean"], sd = b[,"sd"],
                   std_effect = abs(b[,"mean"]/b[,"sd"])) |> arrange(desc(std_effect))
cat("\nTop 15 disease-associated genes (by |mean/sd|):\n")
print(head(transform(rank, coef=round(coef,3), sd=round(sd,3), std_effect=round(std_effect,2)), 15),
      row.names = FALSE)
cat(sprintf("\nSpatial autocorrelation eta = %.3f  (95%% CrI [%.3f, %.3f])\n",
            s["eta","mean"], s["eta","2.5%"], s["eta","97.5%"]))
try(ggsave(file.path(outdir,"fig_her2_coef.png"),
           plot_coef(fit, as.matrix(matrix_x), n = 30), width = 9, height = 4.5, dpi = 120))

## ---- Gene-gene interaction network among top genes (Fig 2c) ----
sep("Interaction network among top-8 genes")
top8 <- head(rank$gene, 8)
fit_int <- dsgd(list_y = y, matrix_x = as.matrix(matrix_x[, c(top8,"x","y")]),
                label_list = label, interaction = top8)
try(ggsave(file.path(outdir,"fig_her2_network.png"),
           plot_network(fit_int, top8), width = 6, height = 6, dpi = 120))
si <- rstan::summary(fit_int)$summary
bi <- si[grep("^beta\\[", rownames(si)), ]
np <- length(top8); npair <- choose(np,2)
pair_lbl <- combn(top8, 2, FUN = function(p) paste(p, collapse="*"))
intdf <- data.frame(pair = pair_lbl,
                    mean = bi[(np+1):(np+npair),"mean"],
                    std_effect = abs(bi[(np+1):(np+npair),"mean"]/bi[(np+1):(np+npair),"sd"]))
cat("Interactions retained (|std effect|>1.96):\n")
print(transform(intdf[intdf$std_effect>1.96,], mean=round(mean,3), std_effect=round(std_effect,2)),
      row.names = FALSE)

## ---- Leave-one-patient-out prediction (Fig 2f) ----
sep("Leave-one-patient-out prediction: train A,B,C -> test D")
tr <- which(label != 4); te <- which(label == 4)
fit_tr <- dsgd(list_y = y[tr], matrix_x = as.matrix(matrix_x[tr,]), label_list = label[tr])
pred <- predict(fit_tr, matrix_x[te, ], sweep = 60)
acc  <- mean(pred == y[te]); base <- max(mean(y[te]), 1-mean(y[te]))
cat(sprintf("Test patient D: %d spots | DiSTect acc = %.3f | majority-class baseline = %.3f\n",
            length(te), acc, base))
print(table(truth = y[te], predicted = pred))

sep("HER2 REPRODUCTION DONE")
cat("Figures: fig_her2_coef.png, fig_her2_network.png in", outdir, "\n")
