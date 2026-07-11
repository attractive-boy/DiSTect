#!/usr/bin/env Rscript

# Small targeted ADVI grid. This separates single-slice algorithm/scale/X tests
# from multiple-slice covariance tests and intentionally avoids a full Cartesian
# product.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "repro/paper_reproduce/run_targeted_advi_grid.R"
root <- normalizePath(file.path(dirname(script_file), "../.."))
source(file.path(root, "repro/paper_reproduce/paper_models.R"))

outdir <- file.path(root, "repro/paper_reproduce/output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_rep <- as.integer(Sys.getenv("N_REP", "10"))
iter <- as.integer(Sys.getenv("ADVI_ITER", "3000"))
run_single <- tolower(Sys.getenv("RUN_SINGLE", "true")) %in% c("1", "true", "yes")
run_multi <- tolower(Sys.getenv("RUN_MULTI", "true")) %in% c("1", "true", "yes")
algorithms <- strsplit(Sys.getenv("ALGORITHMS", "meanfield,fullrank"), ",")[[1]]

targeted_gibbs <- function(X, coords, beta, eta, scale = "average", sweeps = 2000L,
                           offset = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  y <- rbinom(n, 1L, 0.5)
  xb <- as.numeric(X %*% beta) + offset
  A <- adjacency_matrix(coords)
  deg <- pmax(Matrix::rowSums(A), 1)
  parity <- (as.integer(coords[, 1]) + as.integer(coords[, 2])) %% 2L
  blocks <- list(which(parity == 0L), which(parity == 1L))
  for (s in seq_len(sweeps)) {
    for (idx in blocks) {
      neigh <- as.numeric(A %*% y)
      if (scale == "average") neigh <- neigh / deg
      y[idx] <- rbinom(length(idx), 1L, plogis(xb[idx] + eta * neigh[idx]))
    }
  }
  y
}

fit_single_grid <- function(y, X, coords, algorithm, scale, seed) {
  neigh <- neighbor_sum_count(coords, y)
  c_term <- if (scale == "average") neigh$avg else neigh$sum
  dat <- list(N = nrow(X), P = ncol(X), x = X, y = as.integer(y), c_avg = c_term)
  model <- get_stan_model("paper_single", paper_single_stan)
  vb(model, data = dat, algorithm = algorithm, seed = seed, iter = iter,
     tol_rel_obj = 1e-5, refresh = 0)
}

fit_multi_grid <- function(d, algorithm, seed) {
  lab <- as.integer(factor(d$label))
  c_term <- neighbor_sum_count(d$coords, d$y, label = lab)$avg
  dat <- list(N = nrow(d$X), P = ncol(d$X), S = length(unique(lab)), x = d$X,
              y = as.integer(d$y), c_avg = c_term, label = lab,
              u_mu = rep(0, length(unique(lab))))
  model <- get_stan_model("paper_multiple", paper_multiple_stan)
  vb(model, data = dat, algorithm = algorithm, seed = seed, iter = iter,
     tol_rel_obj = 1e-5, refresh = 0)
}

eta_summary <- function(fit, truth) {
  s <- rstan::summary(fit)$summary
  eta_draws <- as.numeric(rstan::extract(fit, pars = "eta", permuted = TRUE)$eta)
  data.frame(
    eta_mean = s["eta", "mean"],
    posterior_sd = s["eta", "sd"],
    stan_mcse = if ("se_mean" %in% colnames(s)) s["eta", "se_mean"] else NA_real_,
    draw_mcse = sd(eta_draws) / sqrt(length(eta_draws)),
    n_draws = length(eta_draws),
    q025 = s["eta", "2.5%"],
    q975 = s["eta", "97.5%"],
    abs_bias = abs(s["eta", "mean"] - truth),
    covered = s["eta", "2.5%"] <= truth && s["eta", "97.5%"] >= truth
  )
}

single_rows <- list()
if (run_single) {
  side <- 30L
  p <- 20L
  eta <- 0.4
  beta <- c(1, 2, 3, -4, -5, rep(0, p - 5L))
  coords <- grid_coords(side)
  set.seed(910001L)
  fixed_X <- matrix(rnorm(nrow(coords) * p), nrow(coords), p)
  colnames(fixed_X) <- paste0("beta", seq_len(p))
  configs <- expand.grid(algorithm = algorithms,
                         neighbor_scale = c("average", "sum"),
                         x_mode = c("fixed", "regenerated"),
                         stringsAsFactors = FALSE)
  cat(sprintf("Single-slice grid: %d configs x %d reps\n", nrow(configs), n_rep))
  for (k in seq_len(nrow(configs))) {
    cfg <- configs[k, ]
    for (rep_id in seq_len(n_rep)) {
      seed <- 911000L + k * 100L + rep_id
      if (cfg$x_mode == "fixed") {
        X <- fixed_X
      } else {
        set.seed(seed)
        X <- matrix(rnorm(nrow(coords) * p), nrow(coords), p)
        colnames(X) <- paste0("beta", seq_len(p))
      }
      y <- targeted_gibbs(X, coords, beta, eta, scale = cfg$neighbor_scale,
                          seed = seed + 50000L)
      elapsed <- system.time(
        fit <- fit_single_grid(y, X, coords, cfg$algorithm, cfg$neighbor_scale,
                               seed = seed + 90000L)
      )["elapsed"]
      row <- cbind(cfg, rep = rep_id, eta_summary(fit, eta), seconds = as.numeric(elapsed))
      single_rows[[length(single_rows) + 1L]] <- row
      cat(sprintf("  single %s/%s/%s rep %d: %.1fs\n", cfg$algorithm,
                  cfg$neighbor_scale, cfg$x_mode, rep_id, elapsed))
    }
  }
}

multi_rows <- list()
if (run_multi) {
  settings <- data.frame(
    sigma2 = c(0.1, 0.1, 0.4, 0.4),
    rho = c(0.1, 0.4, 0.1, 0.4)
  )
  configs <- merge(settings, data.frame(algorithm = algorithms), all = TRUE)
  cat(sprintf("Multiple-slice grid: %d configs x %d reps\n", nrow(configs), n_rep))
  for (k in seq_len(nrow(configs))) {
    cfg <- configs[k, ]
    for (rep_id in seq_len(n_rep)) {
      seed <- 921000L + k * 100L + rep_id
      d <- sim2_dataset(sigma2 = cfg$sigma2, rho = cfg$rho, seed = seed)
      elapsed <- system.time(
        fit <- fit_multi_grid(d, cfg$algorithm, seed = seed + 90000L)
      )["elapsed"]
      row <- cbind(cfg, rep = rep_id, eta_summary(fit, d$eta), seconds = as.numeric(elapsed))
      multi_rows[[length(multi_rows) + 1L]] <- row
      cat(sprintf("  multi %s sigma2=%.1f rho=%.1f rep %d: %.1fs\n",
                  cfg$algorithm, cfg$sigma2, cfg$rho, rep_id, elapsed))
    }
  }
}

aggregate_rows <- function(raw, keys) {
  groups <- split(raw, interaction(raw[keys], drop = TRUE, lex.order = TRUE))
  out <- do.call(rbind, lapply(groups, function(d) {
    id <- d[1, keys, drop = FALSE]
    cbind(id,
      n = nrow(d),
      eta_bias = mean(d$abs_bias),
      empirical_sd = sd(d$eta_mean),
      mean_posterior_sd = mean(d$posterior_sd),
      mean_stan_mcse = mean(d$stan_mcse, na.rm = TRUE),
      mean_draw_mcse = mean(d$draw_mcse),
      coverage = mean(d$covered),
      mean_seconds = mean(d$seconds))
  }))
  rownames(out) <- NULL
  out
}

if (length(single_rows)) {
  single_raw <- do.call(rbind, single_rows)
  single_agg <- aggregate_rows(single_raw, c("algorithm", "neighbor_scale", "x_mode"))
  write.csv(single_raw, file.path(outdir, sprintf("targeted_single_advi_raw_n%s.csv", n_rep)), row.names = FALSE)
  write.csv(single_agg, file.path(outdir, sprintf("targeted_single_advi_aggregate_n%s.csv", n_rep)), row.names = FALSE)
  cat("\nSingle-slice aggregate:\n")
  print(single_agg, row.names = FALSE)
}

if (length(multi_rows)) {
  multi_raw <- do.call(rbind, multi_rows)
  multi_agg <- aggregate_rows(multi_raw, c("algorithm", "sigma2", "rho"))
  write.csv(multi_raw, file.path(outdir, sprintf("targeted_multi_advi_raw_n%s.csv", n_rep)), row.names = FALSE)
  write.csv(multi_agg, file.path(outdir, sprintf("targeted_multi_advi_aggregate_n%s.csv", n_rep)), row.names = FALSE)
  cat("\nMultiple-slice aggregate:\n")
  print(multi_agg, row.names = FALSE)
}
