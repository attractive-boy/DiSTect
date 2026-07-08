# Count-likelihood layer (Direction 2): make the covariate representation valid at
# single-cell / low-depth resolution, where log(X+1) + Gaussian breaks down.
#
# Provides a fast, per-gene Negative-Binomial normalization (sctransform-lite):
# size-factor offset + analytic Pearson residuals, which stabilize the mean-variance
# relationship and de-sparsify low counts. Drop-in replacement for log1p() as the
# covariate builder feeding the sparse fitters / PG-CAVI.
#
# `counts`: spots (rows) x genes (cols), raw integer UMI counts.

# median-of-ratios-ish size factors (robust to a few high-count genes)
size_factors <- function(counts) {
  tot <- rowSums(counts)
  sf  <- tot / exp(mean(log(pmax(tot, 1))))   # geometric-mean normalized depth
  pmax(sf, 1e-8)
}

# Per-gene NB Pearson residuals with a shared-ish dispersion (method-of-moments).
# mu_ig = s_i * lambda_g ; Var = mu + mu^2/theta_g ; resid = (x-mu)/sqrt(Var).
nb_pearson <- function(counts, clip = sqrt(nrow(counts))) {
  counts <- as.matrix(counts)
  s   <- size_factors(counts)
  lam <- colSums(counts) / sum(s)                    # gene rate
  Mu  <- outer(s, lam)                               # N x G expected
  # method-of-moments dispersion per gene: Var_g = mean(mu) + mean(mu)^2/theta
  m   <- colMeans(Mu)
  v   <- apply(counts, 2, var)
  theta <- ifelse(v > m, m^2 / pmax(v - m, 1e-8), 1e6)  # large theta -> Poisson
  Var <- sweep(Mu + sweep(Mu^2, 2, theta, "/"), 2, 1, "*")
  R   <- (counts - Mu) / sqrt(pmax(Var, 1e-8))
  R[R >  clip] <-  clip
  R[R < -clip] <- -clip
  R
}

# Convenience: build a design matrix from raw counts using either transform,
# so pipelines can switch method="nb" vs method="log" for the count-regime sim.
build_covariates <- function(counts, method = c("nb", "log"), scale = TRUE) {
  method <- match.arg(method)
  M <- if (method == "nb") nb_pearson(counts) else log1p(as.matrix(counts))
  if (scale) M <- scale(M)
  M[is.na(M)] <- 0
  M
}

# --- Full latent-NB measurement model (M2 headline, documented) --------------
# The residual route above is the fast/robust version. The premium version treats
# expression as latent: x_ig ~ NB(mu = s_i * exp(z_ig), theta_g), z_ig the latent
# log-expression entering the disease logistic. With PG augmentation for the NB
# (Polson et al. 2013) the z-update stays conditionally Gaussian, so it composes
# with fit_polyagamma.R at O(N*G). Implemented as an extension in M2.
