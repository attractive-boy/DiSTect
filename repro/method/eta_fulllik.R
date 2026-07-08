# Full-likelihood de-biasing of the spatial autocorrelation eta (Direction B).
#
# DiSTect fits by conditional pseudolikelihood (neighbor labels = observed y),
# which is fast but ATTENUATES eta (we observe ~0.4 recovered vs 1.0 true). This
# module corrects eta by Monte-Carlo moment matching on the spatial sufficient
# statistic S = sum_i y_i * c_i under the FULL autologistic model, holding beta fixed.
#
# Requires neighbor_sum() and build_adjacency() from method/neighbors.R.

# one Gibbs sweep of the autologistic given (Xbeta, eta) using sparse adjacency A
.gibbs_sweep <- function(y, Xbeta, eta, A) {
  # random-scan update; vectorized approximate sweep via current neighbor sums
  cc <- as.numeric(A %*% y)
  p  <- plogis(Xbeta + eta * cc)
  rbinom(length(y), 1, p)
}

# E[S(eta)] by simulating the field to (approx) equilibrium, averaged over draws
.expected_S <- function(Xbeta, eta, A, coords, label, burn = 30, draws = 20) {
  n <- length(Xbeta); y <- rbinom(n, 1, plogis(Xbeta))
  for (s in 1:burn) y <- .gibbs_sweep(y, Xbeta, eta, A)
  S <- numeric(draws)
  for (d in 1:draws) {
    y <- .gibbs_sweep(y, Xbeta, eta, A)
    S[d] <- sum(y * neighbor_sum(coords, y, label = label))
  }
  mean(S)
}

# Debias: solve E[S(eta)] = S_obs by bisection on eta in [0, eta_max].
# `offset` = fitted intercept (baseline log-odds); pass it so the simulated field
# has the correct disease rate.
debias_eta <- function(y, X, coords, beta, offset = 0, label = NULL,
                       eta_max = 8, tol = 0.02, maxit = 25, seed = 1) {
  set.seed(seed)
  A     <- build_adjacency(coords, label = label)
  Xbeta <- as.numeric(as.matrix(X) %*% beta) + offset
  S_obs <- sum(y * neighbor_sum(coords, y, label = label))
  lo <- 0; hi <- eta_max
  f  <- function(e) .expected_S(Xbeta, e, A, coords, label) - S_obs
  flo <- f(lo); fhi <- f(hi)
  if (flo > 0) return(list(eta = 0, note = "S_obs below null; eta~0"))
  if (fhi < 0) return(list(eta = eta_max, note = "clipped at eta_max"))
  for (it in 1:maxit) {
    mid <- (lo + hi) / 2; fm <- f(mid)
    if (abs(fm) / max(S_obs, 1) < tol) break
    if (fm < 0) lo <- mid else hi <- mid
  }
  list(eta = (lo + hi) / 2, S_obs = S_obs, iters = it)
}
