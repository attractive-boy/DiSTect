# M0 sparse O(N) fitters -- drop-in equivalents of dsgd_single / dsgd_multiple.
# Identical priors and vb() call as R/dsgd.R; only the mean is reformulated:
#   dense:  mu[i] = x[i,1:P]*beta + eta * sum_{j in N(i)} y[j]   (rebuilt via NxN loop)
#   sparse: mu    = X*beta + eta * c   with c precomputed once (constant vector)
# The uniform priors on w/eta that the original adds inside its loops are constant
# in log-density, so adding them once here leaves the posterior identical.
suppressMessages(library(rstan))
# NOTE: driver scripts must source("repro/m0/neighbors.R") before this file
# (provides neighbor_sum() and build_adjacency()).

.SPARSE_SINGLE <- "
data {
  int<lower=0> N; int<lower=0> P;
  matrix[N,P] x;
  int y[N];
  vector[N] c;                      // precomputed neighbor sum
}
parameters {
  vector[P] beta;
  real<lower=0,upper=8> eta;
  vector<lower=0>[P] beta_gamma;
  real<lower=0,upper=1> w;
}
model {
  vector[N] mu;
  w ~ uniform(0,1);
  eta ~ uniform(0,8);
  for (i in 1:P) {
    beta_gamma[i] ~ inv_gamma(5,50);
    target += log_sum_exp(log(1-w)+normal_lpdf(beta[i]|0,0.000001*beta_gamma[i]),
                          log(w)  +normal_lpdf(beta[i]|0,beta_gamma[i]));
  }
  mu = x*beta + eta*c;
  y ~ bernoulli_logit(mu);
}"

.SPARSE_MULTI <- "
data {
  int<lower=0> N; int<lower=0> P; int<lower=0> S;
  matrix[N,P] x;
  int y[N];
  vector[N] c;
  int label[N];
  vector[S] u_mu;
}
parameters {
  vector[P] beta;
  real<lower=0,upper=8> eta;
  vector<lower=0>[P] beta_gamma;
  real<lower=0,upper=1> w;
  vector<lower=-1,upper=1>[S] U;
  real<lower=0,upper=1> sigma_square;
  real<lower=0,upper=1> rho;
}
model {
  vector[N] mu;
  matrix[S,S] covar_matrix;
  rho ~ uniform(0,1);
  sigma_square ~ inv_gamma(5,50);
  for (i in 1:S) for (j in 1:S)
    covar_matrix[i,j] = (i==j) ? sigma_square : sigma_square*rho;
  U ~ multi_normal(u_mu, covar_matrix);
  w ~ uniform(0,1);
  eta ~ uniform(0,8);
  for (i in 1:P) {
    beta_gamma[i] ~ inv_gamma(5,50);
    target += log_sum_exp(log(1-w)+normal_lpdf(beta[i]|0,0.000001*beta_gamma[i]),
                          log(w)  +normal_lpdf(beta[i]|0,beta_gamma[i]));
  }
  for (i in 1:N) mu[i] = dot_product(x[i], beta) + eta*c[i] + U[label[i]];
  y ~ bernoulli_logit(mu);
}"

# cache compiled models across calls in the benchmark
.sm_cache <- new.env()
.get_model <- function(code, name) {
  if (is.null(.sm_cache[[name]])) .sm_cache[[name]] <- stan_model(model_code = code)
  .sm_cache[[name]]
}

# matrix_x: covariates in cols 1:P, coordinates in the last two columns (as in dsgd()).
dsgd_sparse_single <- function(list_y, matrix_x, label = NULL, method = "VI",
                               iter = 3000, nwarmup = 400, nchain = 2) {
  P      <- ncol(matrix_x) - 2
  coords <- matrix_x[, (P+1):(P+2), drop = FALSE]
  X      <- matrix_x[, 1:P, drop = FALSE]
  c_vec  <- neighbor_sum(coords, list_y, label = label)
  data   <- list(N = nrow(X), P = P, x = X, y = list_y, c = c_vec)
  m      <- .get_model(.SPARSE_SINGLE, "single")
  if (method == "NUTS")
    sampling(m, data = data, iter = iter, warmup = nwarmup, chains = nchain, seed = 128, refresh = 0)
  else
    vb(m, data = data, algorithm = "fullrank", seed = 128, iter = iter, tol_rel_obj = 0.00001)
}

dsgd_sparse_multiple <- function(list_y, matrix_x, label_list, method = "VI",
                                 iter = 3000, nwarmup = 400, nchain = 2) {
  P      <- ncol(matrix_x) - 2
  coords <- matrix_x[, (P+1):(P+2), drop = FALSE]
  X      <- matrix_x[, 1:P, drop = FALSE]
  c_vec  <- neighbor_sum(coords, list_y, label = label_list)
  S      <- length(unique(label_list))
  lab    <- as.integer(factor(label_list))
  data   <- list(N = nrow(X), P = P, S = S, x = X, y = list_y, c = c_vec,
                 label = lab, u_mu = rep(0, S))
  m      <- .get_model(.SPARSE_MULTI, "multi")
  if (method == "NUTS")
    sampling(m, data = data, iter = iter, warmup = nwarmup, chains = nchain, seed = 128, refresh = 0)
  else
    vb(m, data = data, algorithm = "fullrank", seed = 128, iter = iter, tol_rel_obj = 0.00001)
}

# O(N*k*sweep) Gibbs predictor: y_new is updated each sweep, so we use the sparse
# adjacency A (contrast with R/prediction.R which rebuilds an NxN matrix per sweep).
predict_sparse <- function(fit, data, A, sweep = 100) {
  n     <- nrow(data)
  P     <- ncol(data) - 2
  cf    <- rstan::summary(fit)$summary[, "mean"]
  beta  <- cf[1:P]; eta <- cf[P + 1]
  Xbeta <- as.numeric(as.matrix(data[, 1:P]) %*% beta)
  y_new <- rbinom(n, 1, 0.5)
  for (s in 1:sweep) {
    cvec <- as.numeric(A %*% y_new)          # neighbor sum for all spots at once
    p    <- plogis(Xbeta + eta * cvec)
    y_new <- rbinom(n, 1, p)
  }
  y_new
}
