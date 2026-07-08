# Calibrated spatial prediction (Direction B).
#
# R/prediction.R plugs in posterior-mean (beta_hat, eta_hat) and ignores posterior
# uncertainty and the random effect U -> point predictions only, and we observed
# modest, possibly mis-calibrated accuracy. This module produces POSTERIOR-PREDICTIVE
# probabilities by propagating posterior draws through the sparse Gibbs field, and
# reports calibration (Brier, ECE) alongside accuracy.
#
# Requires build_adjacency() from method/neighbors.R.

# draw S posterior samples of theta=(beta,eta) from a Gaussian q (PG-CAVI) or a
# stanfit; run the sparse field for each; average predictive probability per spot.
predict_posterior <- function(theta_mean, theta_cov, X, A, S = 100, sweep = 40) {
  d  <- length(theta_mean); P <- ncol(X)
  L  <- chol(theta_cov + diag(1e-8, d))
  n  <- nrow(X)
  prob <- numeric(n)
  for (s in 1:S) {
    th   <- theta_mean + as.numeric(crossprod(L, rnorm(d)))
    beta <- th[1:P]; eta <- th[P + 1]
    Xb   <- as.numeric(as.matrix(X) %*% beta)
    y    <- rbinom(n, 1, plogis(Xb))
    for (k in 1:sweep) y <- rbinom(n, 1, plogis(Xb + eta * as.numeric(A %*% y)))
    prob <- prob + plogis(Xb + eta * as.numeric(A %*% y))
  }
  prob / S
}

# calibration + accuracy metrics
calibration <- function(prob, truth, bins = 10) {
  brier <- mean((prob - truth)^2)
  acc   <- mean((prob > 0.5) == (truth == 1))
  # expected calibration error
  br <- cut(prob, breaks = seq(0, 1, length.out = bins + 1), include.lowest = TRUE)
  ece <- 0
  for (b in levels(br)) {
    idx <- br == b; if (!any(idx)) next
    ece <- ece + (sum(idx)/length(prob)) * abs(mean(prob[idx]) - mean(truth[idx]))
  }
  data.frame(accuracy = acc, brier = brier, ece = ece)
}
