#!/usr/bin/env Rscript

# Decisive test of the pseudolikelihood-vs-full-likelihood efficiency hypothesis.
#
# Diagnosis so far: the local plug-in conditional (pseudo)likelihood -- the model
# in the paper text and the released BayModDSGD package, fit by GLM/ADVI/NUTS --
# gives an eta estimate with coefficient-of-variation ~0.73, ~15x looser than the
# paper's ~0.05. This gap is method-independent (GLM n200, ADVI n200, div=0 NUTS
# all agree) and is NOT closed by sampler tuning, priors, covariate scale, or the
# sum-vs-avg convention (which only rescales eta).
#
# Hypothesis: the paper's table-generating "proposed" fit is a FULL-likelihood
# autologistic inference, which is far more efficient for the dependence
# parameter eta than pseudolikelihood -- but which Stan CANNOT do (intractable
# normalizing constant). This script fits the proper sum-convention autologistic
# joint by Geyer-Thompson Monte Carlo MLE and compares SE(eta) to the GLM
# pseudolikelihood SE(eta) on the SAME data.
#
#   Joint:  P(y | beta, eta) proportional to exp( beta . (X^T y) + eta * T(y) )
#   T(y) = sum over edges (i~j) of y_i y_j   (concordant-edge count)
#   Conditional: P(y_i=1 | rest) = logistic( x_i . beta + eta * sum_{j in N(i)} y_j )
#
# If SE(eta)_MCMLE << SE(eta)_PL (CV approaching ~0.05), the paper's tight eta is
# explained by full-likelihood efficiency, confirming the tables cannot come from
# the Stan pseudolikelihood the paper describes.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
source(file.path(root, "repro/paper_reproduce/paper_models.R"))  # grid_coords, adjacency_matrix

outdir <- file.path(root, "repro/paper_reproduce/output")

n_data   <- as.integer(Sys.getenv("N_DATA", "1"))     # number of independent datasets
eta_true <- as.numeric(Sys.getenv("ETA_TRUE", "0.4")) # eta on the SUM scale
side     <- as.integer(Sys.getenv("SIDE", "30"))
p        <- as.integer(Sys.getenv("P", "20"))
m_sim    <- as.integer(Sys.getenv("M_SIM", "1500"))   # MC fields per MCMLE iter
mcmle_it <- as.integer(Sys.getenv("MCMLE_IT", "5"))
burn     <- as.integer(Sys.getenv("BURN", "1000"))
thin     <- as.integer(Sys.getenv("THIN", "3"))

beta_true <- c(1, 2, 3, -4, -5, rep(0, p - 5L))
coords <- grid_coords(side)
n <- nrow(coords)
A <- adjacency_matrix(coords)                 # symmetric 0/1 adjacency (each edge twice)

# Gibbs sampler for the SUM-convention joint at (beta, eta): returns final field.
gibbs_sum <- function(X, beta, eta, sweeps, y0 = NULL) {
  xb <- as.numeric(X %*% beta)
  y <- if (is.null(y0)) rbinom(n, 1L, 0.5) else y0
  parity <- (as.integer(round(coords[, 1])) + as.integer(round(coords[, 2]))) %% 2L
  blocks <- list(which(parity == 0L), which(parity == 1L))
  for (s in seq_len(sweeps)) {
    for (idx in blocks) {
      csum <- as.numeric(A %*% y)
      y[idx] <- rbinom(length(idx), 1L, plogis(xb[idx] + eta * csum[idx]))
    }
  }
  y
}

# Simulate m fields at theta=(beta,eta); return matrix of sufficient stats (m x (P+1)).
sim_suffstats <- function(X, beta, eta, m, burn, thin) {
  xb <- as.numeric(X %*% beta)
  parity <- (as.integer(round(coords[, 1])) + as.integer(round(coords[, 2]))) %% 2L
  blocks <- list(which(parity == 0L), which(parity == 1L))
  y <- rbinom(n, 1L, 0.5)
  step <- function(y) { for (idx in blocks) { csum <- as.numeric(A %*% y); y[idx] <- rbinom(length(idx), 1L, plogis(xb[idx] + eta * csum[idx])) }; y }
  for (s in seq_len(burn)) y <- step(y)
  S <- matrix(0, m, p + 1L)
  for (k in seq_len(m)) {
    for (s in seq_len(thin)) y <- step(y)
    S[k, 1:p] <- as.numeric(crossprod(X, y))
    S[k, p + 1L] <- 0.5 * sum(y * (A %*% y))    # T(y) = concordant edges
  }
  S
}

suff_obs <- function(X, y) c(as.numeric(crossprod(X, y)), 0.5 * sum(y * (A %*% y)))

# Geyer-Thompson MCMLE starting from theta0 (pseudolikelihood estimate).
mcmle_fit <- function(X, y_obs, theta0) {
  Sobs <- suff_obs(X, y_obs)
  theta_cur <- theta0
  for (it in seq_len(mcmle_it)) {
    S <- sim_suffstats(X, theta_cur[1:p], theta_cur[p + 1L], m_sim, burn, thin)
    Sc <- colMeans(S)
    # negative log MC-likelihood-ratio r(theta) and its gradient
    negr <- function(theta) {
      d <- theta - theta_cur
      lw <- as.numeric(S %*% d)                       # (theta-theta_cur).S_k
      lse <- max(lw) + log(mean(exp(lw - max(lw))))   # log mean exp
      -(sum(d * Sobs) - lse)
    }
    negg <- function(theta) {
      d <- theta - theta_cur
      lw <- as.numeric(S %*% d); w <- exp(lw - max(lw)); w <- w / sum(w)
      -(Sobs - as.numeric(crossprod(S, w)))
    }
    opt <- optim(theta_cur, negr, negg, method = "BFGS",
                 control = list(maxit = 200, reltol = 1e-10))
    theta_cur <- opt$par
  }
  # Fisher information at theta_hat = Cov of sufficient stats (re-simulate clean)
  Sfin <- sim_suffstats(X, theta_cur[1:p], theta_cur[p + 1L], max(m_sim, 3000L), burn, thin)
  I <- cov(Sfin)
  Vinv <- tryCatch(solve(I), error = function(e) MASS::ginv(I))
  se <- sqrt(pmax(diag(Vinv), 0))
  list(theta = theta_cur, se = se, eta = theta_cur[p + 1L], eta_se = se[p + 1L])
}

rows <- list()
cat(sprintf("== Full-likelihood MCMLE vs pseudolikelihood | eta_true(sum)=%.2f n_data=%d m=%d it=%d ==\n\n",
            eta_true, n_data, m_sim, mcmle_it))

for (dd in seq_len(n_data)) {
  set.seed(900000L + dd)
  X <- matrix(rnorm(n * p), n, p)
  y <- gibbs_sum(X, beta_true, eta_true, sweeps = 2000L)
  drate <- mean(y)
  # pseudolikelihood (GLM, sum convention)
  csum <- as.numeric(A %*% y)
  dat <- data.frame(y = y, X, spatial = csum)
  gl <- suppressWarnings(glm(y ~ . - 1, data = dat, family = binomial()))
  co <- coef(summary(gl))
  eta_pl <- unname(coef(gl)["spatial"]); eta_pl_se <- unname(co["spatial", "Std. Error"])
  theta0 <- c(unname(coef(gl)[paste0("X", 1:p)]), eta_pl)
  t0 <- Sys.time()
  mm <- mcmle_fit(X, y, theta0)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rows[[dd]] <- data.frame(
    dataset = dd, disease_rate = round(drate, 3),
    eta_true = eta_true,
    eta_PL = round(eta_pl, 4), eta_PL_se = round(eta_pl_se, 4), CV_PL = round(eta_pl_se / abs(eta_pl), 3),
    eta_MCMLE = round(mm$eta, 4), eta_MCMLE_se = round(mm$eta_se, 4), CV_MCMLE = round(mm$eta_se / abs(mm$eta), 3),
    SE_ratio_PL_over_MCMLE = round(eta_pl_se / mm$eta_se, 2), secs = round(secs, 1),
    stringsAsFactors = FALSE
  )
  cat(sprintf("  d%d rate=%.2f | PL eta=%.3f se=%.3f (CV %.2f) | MCMLE eta=%.3f se=%.3f (CV %.2f) | SEratio=%.1fx | %.0fs\n",
              dd, drate, eta_pl, eta_pl_se, eta_pl_se/abs(eta_pl),
              mm$eta, mm$eta_se, mm$eta_se/abs(mm$eta), eta_pl_se/mm$eta_se, secs))
}

raw <- do.call(rbind, rows)
raw_path <- file.path(outdir, sprintf("mcmle_vs_pl_eta%s_n%s.csv", gsub("\\.", "", as.character(eta_true)), n_data))
write.csv(raw, raw_path, row.names = FALSE)
cat("\nWrote", raw_path, "\n\n")
cat("Paper eta CV reference (Table C1): ~0.05 (SD 0.021 / eta 0.42). Local PL CV ~0.73.\n")
print(raw, row.names = FALSE)
