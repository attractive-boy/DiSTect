# Polya-Gamma CAVI for the sparse autologistic disease model.
#
# M0 showed the fit reduces to a logistic regression with design z_i = [x_i, c_i]
# (c_i = precomputed neighbor sum). Polya-Gamma augmentation makes that logistic
# model conditionally Gaussian, so we get CLOSED-FORM coordinate-ascent VI
# (Durante & Rigon 2019) -- no HMC/ADVI gradient ascent, O(N*P^2) per iteration,
# and numerically stable (monotone ELBO). This is the "optimization route" that
# also sidesteps the ADVI instabilities (Pareto-k, premature convergence) seen in M0.
#
# Prior here is a ridge (Normal) prior -- the demonstrator for speed/stability.
# Spike-and-slab selection layers on by replacing the fixed prior precision `p0`
# with per-coordinate precisions updated from inclusion probabilities (see
# `pgcavi_spikeslab_hook` note at the bottom); left as the documented extension.

# Core CAVI: y in {0,1}, Z is N x d design, Normal(0, prior_sd^2) prior on all coefs
# (prior_sd may be a length-d vector). Returns posterior mean/sd + ELBO trace.
fit_pgcavi <- function(y, Z, prior_sd = 10, max_iter = 500, tol = 1e-6) {
  N <- nrow(Z); d <- ncol(Z)
  if (length(prior_sd) == 1) prior_sd <- rep(prior_sd, d)
  P0    <- diag(1 / prior_sd^2, d)          # prior precision
  kappa <- y - 0.5
  Ztk   <- crossprod(Z, kappa)              # Z' (y - 1/2)
  mu    <- rep(0, d); xi <- rep(1, N)
  elbo_prev <- -Inf; trace <- numeric(0)
  lambda <- function(x) {                   # tanh(x/2)/(2x), stable at x->0
    out <- tanh(x / 2) / (2 * x); out[!is.finite(out)] <- 0.25; out
  }
  for (it in 1:max_iter) {
    a     <- lambda(xi)                      # E[omega_i]
    Prec  <- crossprod(Z, Z * a) + P0        # Z' diag(a) Z + P0
    Sigma <- chol2inv(chol(Prec))
    mu    <- as.numeric(Sigma %*% Ztk)
    # update variational PG params: xi_i^2 = z_i (Sigma + mu mu') z_i'
    ZS    <- Z %*% Sigma
    quad  <- rowSums(ZS * Z) + as.numeric(Z %*% mu)^2
    xi    <- sqrt(pmax(quad, 1e-12))
    # ELBO (up to constants) for convergence monitoring
    elbo  <- sum(y * (Z %*% mu)) - 0.5 * as.numeric(t(mu) %*% P0 %*% mu) +
             0.5 * determinant(Sigma, logarithm = TRUE)$modulus
    trace <- c(trace, elbo)
    if (abs(elbo - elbo_prev) < tol * abs(elbo_prev)) break
    elbo_prev <- elbo
  }
  list(mean = mu, sd = sqrt(diag(Sigma)), Sigma = Sigma,
       elbo = trace, iter = length(trace))
}

# Convenience wrapper matching the sparse model: covariates X (N x P) + neighbor
# sum c -> design [X, c]; last coefficient is eta (spatial autocorrelation).
# `neighbor_sum()` from method/neighbors.R must be sourced by the caller when
# passing coords; or pass c_vec directly.
fit_pgcavi_single <- function(y, X, c_vec, coords = NULL, label = NULL,
                              intercept = TRUE,
                              prior_sd_beta = 5, prior_sd_eta = 10,
                              prior_sd_int = 10, ...) {
  if (missing(c_vec) || is.null(c_vec)) {
    stopifnot(!is.null(coords))
    c_vec <- neighbor_sum(coords, y, label = label)
  }
  P  <- ncol(X)
  gene_names <- colnames(X) %||% paste0("beta", 1:P)
  Z  <- cbind(as.matrix(X), c = c_vec)
  psd <- c(rep(prior_sd_beta, P), prior_sd_eta)
  if (intercept) { Z <- cbind(intercept = 1, Z); psd <- c(prior_sd_int, psd) }
  fit <- fit_pgcavi(y, Z, prior_sd = psd, ...)
  nm  <- if (intercept) c("intercept", gene_names, "eta") else c(gene_names, "eta")
  tab <- data.frame(param = nm, mean = fit$mean, sd = fit$sd,
                    std_effect = abs(fit$mean / fit$sd))
  bidx <- if (intercept) 2:(P + 1) else 1:P
  list(table = tab, beta = fit$mean[bidx], eta = fit$mean[length(fit$mean)],
       intercept = if (intercept) fit$mean[1] else 0,
       iter = fit$iter, elbo = fit$elbo, raw = fit)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- Spike-and-slab extension (documented hook, not yet enabled) -------------
# To turn the ridge prior into NMIG spike-and-slab selection: after the Sigma/mu
# update, for each gene coordinate j compute inclusion prob
#   gamma_j = sigmoid( logit(w) + 0.5*log(v_spike/v_slab)
#                      + 0.5*(mu_j^2 + Sigma_jj)*(1/v_spike - 1/v_slab) )
# then set prior precision p0_j = gamma_j/v_slab + (1-gamma_j)/v_spike and update
# w ~ Beta from sum(gamma). This keeps every update closed-form (still O(N*P^2)).
