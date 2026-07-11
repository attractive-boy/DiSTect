#!/usr/bin/env Rscript

# Falsify-myself experiment.
#
# Claim under test (MY claim): the paper's eta CV ~0.05 is below the Cramer-Rao
# lower bound for a 30x30 single-slice autologistic with the stated design, so no
# estimator can achieve it. If instead some REASONABLE data-generation setting
# pushes the CRLB down to ~0.05, then I am wrong and the paper is fine.
#
# eta's Fisher information = Var of the sufficient statistic T(y) = concordant
# edges, under the model. Strong covariates pin the field and starve eta of
# information. So the way to rescue the paper is a WEAKER effective covariate
# signal. We sweep the covariate SD (equivalently signal-beta scale) and, at each
# setting, compute the exact-as-MC CRLB for eta via the sufficient-stat
# covariance (same machinery as MCMLE Fisher info). We also report eta's implied
# best-possible CV.
#
# Decisive question: is there any plausible xsd where CRLB CV(eta) drops from
# ~0.7 (my current design) toward the paper's ~0.05, i.e. Fisher info rises ~200x?

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
source(file.path(root, "repro/paper_reproduce/paper_models.R"))

outdir <- file.path(root, "repro/paper_reproduce/output")

eta_sum   <- as.numeric(Sys.getenv("ETA_SUM", "0.11"))  # sum-scale eta (0.11 ~ avg 0.4)
side      <- as.integer(Sys.getenv("SIDE", "30"))
p         <- as.integer(Sys.getenv("P", "20"))
m_sim     <- as.integer(Sys.getenv("M_SIM", "3000"))
burn      <- as.integer(Sys.getenv("BURN", "1000"))
thin      <- as.integer(Sys.getenv("THIN", "3"))
xsd_grid  <- as.numeric(strsplit(Sys.getenv("XSD", "1.0,0.5,0.25,0.1,0.0"), ",")[[1]])
beta_scale<- as.numeric(Sys.getenv("BETA_SCALE", "1.0"))

coords <- grid_coords(side); n <- nrow(coords)
A <- adjacency_matrix(coords)
parity <- (as.integer(round(coords[, 1])) + as.integer(round(coords[, 2]))) %% 2L
blocks <- list(which(parity == 0L), which(parity == 1L))
beta_base <- c(1, 2, 3, -4, -5, rep(0, p - 5L)) * beta_scale

sim_suffstats <- function(X, beta, eta, m, burn, thin) {
  xb <- as.numeric(X %*% beta)
  y <- rbinom(n, 1L, 0.5)
  step <- function(y) { for (idx in blocks) { csum <- as.numeric(A %*% y); y[idx] <- rbinom(length(idx), 1L, plogis(xb[idx] + eta * csum[idx])) }; y }
  for (s in seq_len(burn)) y <- step(y)
  S <- matrix(0, m, p + 1L); rate <- numeric(m)
  for (k in seq_len(m)) {
    for (s in seq_len(thin)) y <- step(y)
    S[k, 1:p] <- as.numeric(crossprod(X, y)); S[k, p + 1L] <- 0.5 * sum(y * (A %*% y)); rate[k] <- mean(y)
  }
  list(S = S, rate = mean(rate))
}

cat(sprintf("== CRLB(eta) vs covariate scale | eta_sum=%.2f beta_scale=%.2f m=%d ==\n", eta_sum, beta_scale, m_sim))
cat("Paper eta CV target ~0.05. My current-design CRLB CV ~0.7.\n\n")
rows <- list()
for (xsd in xsd_grid) {
  set.seed(12345L)
  X <- matrix(rnorm(n * p, sd = xsd), n, p)
  ss <- sim_suffstats(X, beta_base, eta_sum, m_sim, burn, thin)
  I <- cov(ss$S)                       # Fisher information (P+1 x P+1)
  Vinv <- tryCatch(solve(I), error = function(e) MASS::ginv(I))
  eta_var_crlb <- diag(Vinv)[p + 1L]   # best-possible Var(eta_hat)
  eta_se_crlb <- sqrt(max(eta_var_crlb, 0))
  fisher_eta <- I[p + 1L, p + 1L]      # marginal Var(T(y)) = raw info scale
  lp_sd <- xsd * sqrt(sum((beta_base[1:5])^2))
  rows[[length(rows)+1L]] <- data.frame(
    xsd = xsd, lp_sd = round(lp_sd, 2), disease_rate = round(ss$rate, 3),
    VarT = round(fisher_eta, 1),
    eta_se_CRLB = round(eta_se_crlb, 4),
    eta_CV_CRLB = round(eta_se_crlb / abs(eta_sum), 3),
    stringsAsFactors = FALSE)
  cat(sprintf("  xsd=%.2f (LPsd %.2f) rate=%.2f | Var(T)=%.0f | CRLB SE(eta)=%.4f | CRLB CV=%.3f\n",
              xsd, lp_sd, ss$rate, fisher_eta, eta_se_crlb, eta_se_crlb/abs(eta_sum)))
}
raw <- do.call(rbind, rows)
raw_path <- file.path(outdir, sprintf("crlb_vs_covscale_eta%s.csv", gsub("\\.","",as.character(eta_sum))))
write.csv(raw, raw_path, row.names = FALSE)
cat("\nWrote", raw_path, "\n")
cat("\nInterpretation: if CRLB CV never approaches 0.05 for any plausible xsd, the\n")
cat("paper's eta precision is below the information bound for this design => not an\n")
cat("estimator choice. If some xsd hits ~0.05, I was wrong about the design.\n")
