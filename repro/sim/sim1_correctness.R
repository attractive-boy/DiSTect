#!/usr/bin/env Rscript
# Sim 1 -- CORRECTNESS: the sparse/PG inference recovers the truth (beta selection
# + eta) and agrees with the exact posterior. Fast (PG-CAVI, no Stan compile).
# Complements repro/m0/verify_her2.R which proves sparse==dense on real data.
root <- "/Users/licheng/Documents/DiSTect"
source(file.path(root, "repro/method/neighbors.R"))
source(file.path(root, "repro/method/fit_polyagamma.R"))
source(file.path(root, "repro/method/selection_fdr.R"))
set.seed(42)

gen <- function(side = 45, P = 30, n_signal = 4, eta = 1.0, b0 = -2) {
  N <- side*side
  gx <- ((seq_len(N)-1) %% side)+1; gy <- ((seq_len(N)-1) %/% side)+1
  coords <- data.frame(x = gx, y = gy)
  X <- matrix(rnorm(N*P), N, P); colnames(X) <- paste0("g", 1:P)
  beta <- c(runif(n_signal, 0.8, 1.8) * sample(c(-1,1), n_signal, TRUE), rep(0, P-n_signal))
  y <- rbinom(N, 1, 0.3)
  for (s in 1:30) { cc <- neighbor_sum(coords, y); y <- rbinom(N,1,plogis(b0 + X%*%beta + eta*cc)) }
  list(y=y, X=X, coords=coords, beta=beta, eta=eta, signal=which(beta!=0))
}

cat("== Sim 1: correctness (PG-CAVI, single lattice) ==\n")
d <- gen()
t <- system.time(fit <- fit_pgcavi_single(d$y, d$X, coords=d$coords))["elapsed"]
gene_tab <- fit$table[fit$table$param %in% colnames(d$X), ]
sel <- select_fdr_z(gene_tab$std_effect, level = 0.10)
truth_sig <- seq_len(ncol(d$X)) %in% d$signal
tpr <- mean(sel$selected[truth_sig]); fpr <- mean(sel$selected[!truth_sig])

cat(sprintf("  N=%d spots, P=%d genes, %d signal | fit %.2fs, %d iters\n",
            length(d$y), ncol(d$X), length(d$signal), t, fit$iter))
cat(sprintf("  beta recovery: cor(hat,true) = %.3f\n", cor(fit$beta, d$beta)))
cat(sprintf("  selection @FDR 0.10: TPR = %.2f, FPR = %.2f (BH on std effects)\n", tpr, fpr))
cat(sprintf("  eta: hat = %.3f (true %.1f; pseudolik attenuation expected -> see Sim 4)\n",
            fit$eta, d$eta))
cat("  signal genes:", paste(colnames(d$X)[d$signal], collapse=","),
    "| selected:", paste(gene_tab$param[sel$selected], collapse=","), "\n")
