#!/usr/bin/env Rscript
# Sim 1b -- SPIKE-AND-SLAB selection. PG-CAVI with an NMIG spike-and-slab prior on the gene
# coefficients yields per-gene inclusion probabilities gamma_j and Bayesian-FDR selection,
# tightening the liberal FPR of the mean-field z-BH route (STORY §7: FPR ~ 0.19). Same data as Sim 1.
root <- "/Users/licheng/Documents/DiSTect"
source(file.path(root, "repro/method/neighbors.R"))
source(file.path(root, "repro/method/fit_polyagamma.R"))
source(file.path(root, "repro/method/selection_fdr.R"))
set.seed(42)

gen <- function(side = 45, P = 30, n_signal = 4, eta = 1.0, b0 = -2) {   # identical to sim1
  N <- side*side
  gx <- ((seq_len(N)-1) %% side)+1; gy <- ((seq_len(N)-1) %/% side)+1
  coords <- data.frame(x = gx, y = gy)
  X <- matrix(rnorm(N*P), N, P); colnames(X) <- paste0("g", 1:P)
  beta <- c(runif(n_signal, 0.8, 1.8) * sample(c(-1,1), n_signal, TRUE), rep(0, P-n_signal))
  y <- rbinom(N, 1, 0.3)
  for (s in 1:30) { cc <- neighbor_sum(coords, y); y <- rbinom(N,1,plogis(b0 + X%*%beta + eta*cc)) }
  list(y = y, X = X, coords = coords, beta = beta, signal = which(beta != 0))
}

d <- gen(); truth <- seq_len(ncol(d$X)) %in% d$signal
rate <- function(sel) c(TPR = mean(sel[truth]), FPR = mean(sel[!truth]))

## (A) ridge + z-BH  (current default)
fr <- fit_pgcavi_single(d$y, d$X, coords = d$coords)
gr <- fr$table[fr$table$param %in% colnames(d$X), ]
sr <- select_fdr_z(gr$std_effect, level = 0.10)

## (B) spike-and-slab + Bayesian FDR
fs <- fit_pgcavi_single(d$y, d$X, coords = d$coords, spike_slab = TRUE)
ss <- select_fdr_bayes(fs$gamma, level = 0.10)

cat("== Sim 1b: ridge/z-BH  vs  spike-and-slab/Bayesian-FDR  (target FDR 0.10) ==\n")
cat(sprintf("  signal genes: %s\n", paste(colnames(d$X)[d$signal], collapse = ",")))
cat(sprintf("  ridge  z-BH   : TPR=%.2f  FPR=%.2f  | selected: %s\n",
            rate(sr$selected)["TPR"], rate(sr$selected)["FPR"],
            paste(gr$param[sr$selected], collapse = ",")))
cat(sprintf("  s&s  Bayes-FDR: TPR=%.2f  FPR=%.2f  | selected: %s\n",
            rate(ss$selected)["TPR"], rate(ss$selected)["FPR"],
            paste(colnames(d$X)[ss$selected], collapse = ",")))
cat(sprintf("  inclusion prob gamma: signal mean=%.3f | null mean=%.3f (max null=%.3f)\n",
            mean(fs$gamma[truth]), mean(fs$gamma[!truth]), max(fs$gamma[!truth])))
cat(sprintf("  slab weight w=%.3f (true signal fraction=%.3f) | s&s iters=%d\n",
            fs$raw$w, mean(truth), fs$iter))
cat("Expect: s&s inclusion probs separate signal (~1) from null (~0), and Bayesian-FDR\n")
cat("        cuts the mean-field z-BH FPR while keeping TPR high.\n")
