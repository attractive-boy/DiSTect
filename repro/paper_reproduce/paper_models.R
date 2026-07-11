#!/usr/bin/env Rscript

suppressMessages({
  library(Matrix)
  library(rstan)
})

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rook_offsets <- list(c(1L, 0L), c(-1L, 0L), c(0L, 1L), c(0L, -1L))

grid_coords <- function(side = 30) {
  n <- side * side
  data.frame(
    x = ((seq_len(n) - 1L) %% side) + 1L,
    y = ((seq_len(n) - 1L) %/% side) + 1L
  )
}

neighbor_sum_count <- function(coords, y = NULL, label = NULL, radius = 1L) {
  n <- nrow(coords)
  ix <- as.integer(round(coords[, 1]))
  iy <- as.integer(round(coords[, 2]))
  if (is.null(label)) label <- rep(1L, n)
  key <- paste(ix, iy, label, sep = ",")
  lookup <- setNames(seq_len(n), key)
  c_vec <- numeric(n)
  n_vec <- integer(n)
  for (off in rook_offsets) {
    nkey <- paste(ix + off[1], iy + off[2], label, sep = ",")
    j <- lookup[nkey]
    hit <- !is.na(j)
    n_vec[hit] <- n_vec[hit] + 1L
    if (!is.null(y)) c_vec[hit] <- c_vec[hit] + y[j[hit]]
  }
  list(sum = c_vec, count = n_vec, avg = ifelse(n_vec > 0L, c_vec / n_vec, 0))
}

neighbor_index_list <- function(coords, label = NULL) {
  n <- nrow(coords)
  ix <- as.integer(round(coords[, 1]))
  iy <- as.integer(round(coords[, 2]))
  if (is.null(label)) label <- rep(1L, n)
  key <- paste(ix, iy, label, sep = ",")
  lookup <- setNames(seq_len(n), key)
  out <- vector("list", n)
  for (i in seq_len(n)) out[[i]] <- integer(0)
  for (off in rook_offsets) {
    nkey <- paste(ix + off[1], iy + off[2], label, sep = ",")
    j <- lookup[nkey]
    hit <- which(!is.na(j))
    for (i in hit) out[[i]] <- c(out[[i]], as.integer(j[i]))
  }
  out
}

adjacency_matrix <- function(coords, label = NULL) {
  n <- nrow(coords)
  ix <- as.integer(round(coords[, 1]))
  iy <- as.integer(round(coords[, 2]))
  if (is.null(label)) label <- rep(1L, n)
  key <- paste(ix, iy, label, sep = ",")
  lookup <- setNames(seq_len(n), key)
  from <- integer(0)
  to <- integer(0)
  for (off in rook_offsets) {
    nkey <- paste(ix + off[1], iy + off[2], label, sep = ",")
    j <- lookup[nkey]
    hit <- which(!is.na(j))
    from <- c(from, hit)
    to <- c(to, as.integer(j[hit]))
  }
  sparseMatrix(i = from, j = to, x = 1, dims = c(n, n))
}

gibbs_autologistic <- function(X, coords, beta, eta, sweeps = 2000L, label = NULL,
                               offset = 0, seed = NULL, update = "checkerboard") {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  y <- rbinom(n, 1L, 0.5)
  xb <- as.numeric(X %*% beta) + offset
  if (update == "single_site") {
    neigh <- neighbor_index_list(coords, label = label)
    for (s in seq_len(sweeps)) {
      for (i in sample.int(n)) {
        ni <- neigh[[i]]
        c_avg <- if (length(ni) == 0L) 0 else mean(y[ni])
        y[i] <- rbinom(1L, 1L, plogis(xb[i] + eta * c_avg))
      }
    }
  } else if (update == "checkerboard") {
    parity <- (as.integer(round(coords[, 1])) + as.integer(round(coords[, 2]))) %% 2L
    blocks <- list(which(parity == 0L), which(parity == 1L))
    A <- adjacency_matrix(coords, label = label)
    deg <- pmax(Matrix::rowSums(A), 1)
    for (s in seq_len(sweeps)) {
      for (idx in blocks) {
        c_avg <- as.numeric(A %*% y) / deg
        y[idx] <- rbinom(length(idx), 1L, plogis(xb[idx] + eta * c_avg[idx]))
      }
    }
  } else if (update == "parallel") {
    A <- adjacency_matrix(coords, label = label)
    deg <- pmax(Matrix::rowSums(A), 1)
    for (s in seq_len(sweeps)) {
      c_avg <- as.numeric(A %*% y) / deg
      y <- rbinom(n, 1L, plogis(xb + eta * c_avg))
    }
  } else {
    stop("unknown Gibbs update mode")
  }
  y
}

sim1_dataset <- function(side = 30L, p = 20L, eta = 0.4, seed = 1L) {
  set.seed(seed)
  coords <- grid_coords(side)
  n <- nrow(coords)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("beta", seq_len(p))
  beta <- c(1, 2, 3, -4, -5, rep(0, p - 5L))
  y <- gibbs_autologistic(X, coords, beta, eta, sweeps = 2000L)
  list(X = X, coords = coords, y = y, beta = beta, eta = eta)
}

sim2_dataset <- function(side = 30L, p = 20L, g = 6L, eta = 1.6,
                         sigma2 = 0.1, rho = 0.1, seed = 1L) {
  set.seed(seed)
  beta <- c(1, 2, 3, -4, -5, rep(0, p - 5L))
  Sigma <- matrix(sigma2 * rho, g, g)
  diag(Sigma) <- sigma2
  U <- as.numeric(t(chol(Sigma)) %*% rnorm(g))
  X_all <- NULL
  coord_all <- NULL
  y_all <- integer(0)
  label <- integer(0)
  for (slice in seq_len(g)) {
    coords <- grid_coords(side)
    X <- matrix(rnorm(nrow(coords) * p), nrow(coords), p)
    colnames(X) <- paste0("beta", seq_len(p))
    y <- gibbs_autologistic(X, coords, beta, eta, sweeps = 2000L, offset = U[slice])
    X_all <- rbind(X_all, X)
    coord_all <- rbind(coord_all, coords)
    y_all <- c(y_all, y)
    label <- c(label, rep(slice, nrow(coords)))
  }
  colnames(X_all) <- paste0("beta", seq_len(p))
  list(X = X_all, coords = coord_all, y = y_all, label = label, beta = beta,
       eta = eta, sigma2 = sigma2, rho = rho, U = U)
}

stan_cache <- new.env(parent = emptyenv())

get_stan_model <- function(name, code) {
  if (is.null(stan_cache[[name]])) {
    stan_cache[[name]] <- stan_model(model_code = code)
  }
  stan_cache[[name]]
}

paper_single_stan <- "
data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N,P] x;
  int<lower=0,upper=1> y[N];
  vector[N] c_avg;
}
parameters {
  vector[P] beta;
  real<lower=0,upper=8> eta;
  vector<lower=0>[P] beta_gamma;
  real<lower=0,upper=1> w;
}
model {
  vector[N] mu;
  w ~ uniform(0, 1);
  eta ~ uniform(0, 8);
  for (i in 1:P) {
    beta_gamma[i] ~ inv_gamma(5, 50);
    target += log_sum_exp(log(1-w) + normal_lpdf(beta[i] | 0, 0.000001 * beta_gamma[i]),
                          log(w) + normal_lpdf(beta[i] | 0, beta_gamma[i]));
  }
  mu = x * beta + eta * c_avg;
  y ~ bernoulli_logit(mu);
}
"

paper_naive_stan <- "
data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N,P] x;
  int<lower=0,upper=1> y[N];
}
parameters {
  vector[P] beta;
  vector<lower=0>[P] beta_gamma;
  real<lower=0,upper=1> w;
}
model {
  vector[N] mu;
  w ~ uniform(0, 1);
  for (i in 1:P) {
    beta_gamma[i] ~ inv_gamma(5, 50);
    target += log_sum_exp(log(1-w) + normal_lpdf(beta[i] | 0, 0.000001 * beta_gamma[i]),
                          log(w) + normal_lpdf(beta[i] | 0, beta_gamma[i]));
  }
  mu = x * beta;
  y ~ bernoulli_logit(mu);
}
"

paper_multiple_stan <- "
data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> S;
  matrix[N,P] x;
  int<lower=0,upper=1> y[N];
  vector[N] c_avg;
  int<lower=1,upper=S> label[N];
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
  rho ~ uniform(0, 1);
  sigma_square ~ inv_gamma(5, 50);
  for (i in 1:S) {
    for (j in 1:S) {
      covar_matrix[i,j] = (i == j) ? sigma_square : sigma_square * rho;
    }
  }
  U ~ multi_normal(u_mu, covar_matrix);
  w ~ uniform(0, 1);
  eta ~ uniform(0, 8);
  for (i in 1:P) {
    beta_gamma[i] ~ inv_gamma(5, 50);
    target += log_sum_exp(log(1-w) + normal_lpdf(beta[i] | 0, 0.000001 * beta_gamma[i]),
                          log(w) + normal_lpdf(beta[i] | 0, beta_gamma[i]));
  }
  for (i in 1:N) {
    mu[i] = dot_product(x[i], beta) + eta * c_avg[i] + U[label[i]];
  }
  y ~ bernoulli_logit(mu);
}
"

fit_paper_single <- function(y, X, coords, method = "VI", iter = 3000L,
                             warmup = 1000L, chains = 2L, seed = 128L) {
  c_avg <- neighbor_sum_count(coords, y)$avg
  data <- list(N = nrow(X), P = ncol(X), x = X, y = as.integer(y), c_avg = c_avg)
  model <- get_stan_model("paper_single", paper_single_stan)
  if (method == "NUTS") {
    sampling(model, data = data, iter = iter, warmup = warmup, chains = chains,
             seed = seed, refresh = 0)
  } else {
    vb(model, data = data, algorithm = "fullrank", seed = seed, iter = iter,
       tol_rel_obj = 0.00001, refresh = 0)
  }
}

fit_paper_single_cavg <- function(y, X, c_avg, method = "VI", iter = 3000L,
                                  warmup = 1000L, chains = 2L, seed = 128L) {
  data <- list(N = nrow(X), P = ncol(X), x = X, y = as.integer(y), c_avg = c_avg)
  model <- get_stan_model("paper_single", paper_single_stan)
  if (method == "NUTS") {
    sampling(model, data = data, iter = iter, warmup = warmup, chains = chains,
             seed = seed, refresh = 0)
  } else {
    vb(model, data = data, algorithm = "fullrank", seed = seed, iter = iter,
       tol_rel_obj = 0.00001, refresh = 0)
  }
}

fit_paper_naive <- function(y, X, method = "VI", iter = 3000L,
                            warmup = 1000L, chains = 2L, seed = 128L) {
  data <- list(N = nrow(X), P = ncol(X), x = X, y = as.integer(y))
  model <- get_stan_model("paper_naive", paper_naive_stan)
  if (method == "NUTS") {
    sampling(model, data = data, iter = iter, warmup = warmup, chains = chains,
             seed = seed, refresh = 0)
  } else {
    vb(model, data = data, algorithm = "fullrank", seed = seed, iter = iter,
       tol_rel_obj = 0.00001, refresh = 0)
  }
}

interior_indices <- function(coords) {
  ix <- coords[, 1]
  iy <- coords[, 2]
  which(ix > min(ix) & ix < max(ix) & iy > min(iy) & iy < max(iy))
}

mask_ignorable <- function(coords, n_missing, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  miss <- sample(interior_indices(coords), n_missing)
  observed <- rep(TRUE, nrow(coords))
  observed[miss] <- FALSE
  observed
}

gibbs_missing_indicator <- function(y, coords, gamma, sweeps = 2000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(y)
  m <- rbinom(n, 1L, 0.02) # 1 = missing
  A <- adjacency_matrix(coords)
  deg <- pmax(Matrix::rowSums(A), 1)
  parity <- (as.integer(round(coords[, 1])) + as.integer(round(coords[, 2]))) %% 2L
  blocks <- list(which(parity == 0L), which(parity == 1L))
  for (s in seq_len(sweeps)) {
    for (idx in blocks) {
      m_avg <- as.numeric(A %*% m) / deg
      p <- plogis(gamma[1] + gamma[2] * y[idx] + gamma[3] * m_avg[idx])
      m[idx] <- rbinom(length(idx), 1L, p)
    }
  }
  observed <- m == 0L
  observed
}

fit_missing_firstpass <- function(y_full, X, coords, observed, impute_iter = 1L,
                                  method = "VI", seed = 128L) {
  y_fill <- as.numeric(y_full)
  y_fill[!observed] <- mean(y_full[observed])
  fit <- NULL
  for (it in seq_len(max(1L, impute_iter))) {
    c_all <- neighbor_sum_count(coords, y_fill)$avg
    fit <- fit_paper_single_cavg(y_full[observed], X[observed, , drop = FALSE],
                                 c_all[observed], method = method, seed = seed)
    s <- rstan::summary(fit)$summary
    beta_rows <- grep("^beta\\[", rownames(s), value = TRUE)
    beta_hat <- s[beta_rows, "mean"]
    eta_hat <- s["eta", "mean"]
    xb <- as.numeric(X %*% beta_hat)
    c_all <- neighbor_sum_count(coords, y_fill)$avg
    y_fill[!observed] <- plogis(xb[!observed] + eta_hat * c_all[!observed])
  }
  fit
}

fit_paper_multiple <- function(y, X, coords, label, method = "VI", iter = 3000L,
                               warmup = 1000L, chains = 2L, seed = 128L) {
  lab <- as.integer(factor(label))
  c_avg <- neighbor_sum_count(coords, y, label = lab)$avg
  s <- length(unique(lab))
  data <- list(N = nrow(X), P = ncol(X), S = s, x = X, y = as.integer(y),
               c_avg = c_avg, label = lab, u_mu = rep(0, s))
  model <- get_stan_model("paper_multiple", paper_multiple_stan)
  if (method == "NUTS") {
    sampling(model, data = data, iter = iter, warmup = warmup, chains = chains,
             seed = seed, refresh = 0)
  } else {
    vb(model, data = data, algorithm = "fullrank", seed = seed, iter = iter,
       tol_rel_obj = 0.00001, refresh = 0)
  }
}

fit_summary_table <- function(fit, truth_beta, truth_eta = NULL) {
  s <- rstan::summary(fit)$summary
  beta_rows <- grep("^beta\\[", rownames(s), value = TRUE)
  out <- data.frame(
    param = paste0("beta", seq_along(beta_rows)),
    truth = truth_beta,
    mean = s[beta_rows, "mean"],
    sd = s[beta_rows, "sd"],
    q025 = s[beta_rows, "2.5%"],
    q975 = s[beta_rows, "97.5%"],
    stringsAsFactors = FALSE
  )
  if (!is.null(truth_eta) && "eta" %in% rownames(s)) {
    out <- rbind(out, data.frame(
      param = "eta", truth = truth_eta, mean = s["eta", "mean"], sd = s["eta", "sd"],
      q025 = s["eta", "2.5%"], q975 = s["eta", "97.5%"], stringsAsFactors = FALSE
    ))
  }
  out$bias <- out$mean - out$truth
  out$abs_bias <- abs(out$bias)
  out$covered <- out$q025 <= out$truth & out$q975 >= out$truth
  out
}

aggregate_sim_table <- function(rows) {
  all_rows <- do.call(rbind, rows)
  groups <- split(all_rows, list(all_rows$method, all_rows$eta_setting, all_rows$param),
                  drop = TRUE)
  out <- lapply(groups, function(d) {
    data.frame(
      method = d$method[1],
      eta_setting = d$eta_setting[1],
      param = d$param[1],
      avgBias = mean(d$abs_bias),
      avgSEE = sd(d$mean),
      avgSEM = mean(d$sd),
      avgCI = mean(d$covered),
      n = nrow(d),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}
