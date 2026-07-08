# M0 equivalence check: sparse O(N) fit vs original dense O(N^2) fit on real her2st data.
suppressMessages({ library(rstan) })
options(mc.cores = parallel::detectCores()); set.seed(128)
root <- "/Users/licheng/Documents/DiSTect"; her2 <- "/tmp/her2st/data"
source(file.path(root, "repro/m0/neighbors.R"))
source(file.path(root, "repro/m0/dsgd_sparse.R"))
source(file.path(root, "R/dsgd.R"))            # original dense dsgd_single

## --- minimal her2 loader (same logic as repro/her2_reproduce.R::load_section) ---
load_section <- function(sec) {
  cnt <- read.delim(gzfile(file.path(her2,"ST-cnts",paste0(sec,".tsv.gz"))),
                    row.names = 1, check.names = FALSE)
  lab <- read.delim(file.path(her2,"ST-pat","lbl",paste0(sec,"_labeled_coordinates.tsv")),
                    check.names = FALSE)
  lab$key <- paste0(round(lab$x),"x",round(lab$y)); lab <- lab[!duplicated(lab$key),]
  rownames(lab) <- lab$key
  common <- intersect(rownames(cnt), lab$key); cnt <- cnt[common,]; lab <- lab[common,]
  keep <- !grepl("undetermined", lab$label, ignore.case = TRUE); cnt <- cnt[keep,]; lab <- lab[keep,]
  co <- do.call(rbind, strsplit(rownames(cnt),"x"))
  list(counts = as.matrix(cnt),
       coords = data.frame(x = as.numeric(co[,1]), y = as.numeric(co[,2])),
       disease = as.integer(grepl("cancer", lab$label, ignore.case = TRUE)))
}

cat("== Load her2st section B1 ==\n")
s <- load_section("B1")
X <- log1p(s$counts[, colSums(s$counts) >= 300, drop = FALSE])
hvg <- names(sort(apply(X, 2, var), decreasing = TRUE))[1:15]
matrix_x <- as.matrix(cbind(scale(X[, hvg]), x = s$coords$x, y = s$coords$y))
y <- s$disease
cat(sprintf("  %d spots, %d HVGs, cancer rate %.0f%%\n", nrow(matrix_x), length(hvg), 100*mean(y)))

## --- (1) neighbor engine correctness: grid-hash vs brute-force dist<=1 (dense rule) ---
cat("\n== (1) neighbor_sum vs brute-force dense rule ==\n")
brute <- sapply(seq_along(y), function(i) {
  d <- sqrt((s$coords$x - s$coords$x[i])^2 + (s$coords$y - s$coords$y[i])^2)
  sum(y[d <= 1 & seq_along(y) != i])
})
c_hash <- neighbor_sum(s$coords, y)
cat(sprintf("  max |c_hash - c_brute| = %g   (must be 0)\n", max(abs(c_hash - brute))))
stopifnot(all(c_hash == brute))
cat("  PASS: grid-hash reproduces the dense dist<=1 neighbor sum exactly.\n")

## --- (2) VI: speed + a note on dense premature convergence ---
cat("\n== (2) VI fit: dense vs sparse (speed) ==\n")
t_dense  <- system.time(fit_d <- dsgd_single(y, matrix_x))["elapsed"]
t_sparse <- system.time(fit_s <- dsgd_sparse_single(y, matrix_x))["elapsed"]
cat(sprintf("  wall-clock (compile+fit): dense = %.1fs  sparse = %.1fs  (%.2fx)\n",
            t_dense, t_sparse, t_dense/t_sparse))
cat("  NOTE: dense's N^2-repeated uniform priors inflate the ELBO magnitude (~1.5e5),\n")
cat("        so VI's *relative* tol triggers PREMATURE convergence; sparse ELBO (~1.2e2)\n")
cat("        optimizes honestly. Hence VI point estimates need not match -> use NUTS below.\n")

## --- (3) EQUIVALENCE via NUTS (exact posterior, free of ELBO-scale artifact) ---
cat("\n== (3) NUTS equivalence: dense vs sparse target the SAME posterior ==\n")
seed <- 128                                   # dsgd_single's NUTS branch reads global `seed`
fit_dN <- dsgd_single(y, matrix_x, method = "NUTS", niter = 600, nwarmup = 300, nchain = 2)
fit_sN <- dsgd_sparse_single(y, matrix_x, method = "NUTS", iter = 600, nwarmup = 300, nchain = 2)
sdN <- rstan::summary(fit_dN)$summary; ssN <- rstan::summary(fit_sN)$summary
bdN <- sdN[grep("^beta\\[", rownames(sdN)), "mean"]
bsN <- ssN[grep("^beta\\[", rownames(ssN)), "mean"]
cmp <- data.frame(gene = hvg, dense_NUTS = round(bdN,3), sparse_NUTS = round(bsN,3))
print(cmp, row.names = FALSE)
cat(sprintf("\n  eta:  dense = %.3f   sparse = %.3f\n", sdN["eta","mean"], ssN["eta","mean"]))
cat(sprintf("  max|d beta| = %.4f | cor(beta) = %.4f | top10 rank agreement = %d/10\n",
            max(abs(bdN-bsN)), cor(bdN,bsN),
            length(intersect(hvg[order(-abs(bdN))][1:10], hvg[order(-abs(bsN))][1:10]))))
ok <- max(abs(bdN-bsN)) < 0.15 && cor(bdN,bsN) > 0.98
cat(sprintf("\n== EQUIVALENCE %s ==\n",
            ifelse(ok, "PASS: sparse == dense within MCMC error", "CHECK: gap exceeds MCMC noise")))
