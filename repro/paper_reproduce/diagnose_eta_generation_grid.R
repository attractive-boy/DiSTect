#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
source(file.path(root, "repro/paper_reproduce/paper_models.R"))

outdir <- file.path(root, "repro/paper_reproduce/output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_rep <- as.integer(Sys.getenv("N_REP", "50"))
eta_values <- as.numeric(strsplit(Sys.getenv("ETA_VALUES", "0.4,1.6,2.8"), ",")[[1]])
sweeps_grid <- as.integer(strsplit(Sys.getenv("SWEEPS", "20,100,500,2000"), ",")[[1]])
updates <- strsplit(Sys.getenv("UPDATES", "parallel,checkerboard"), ",")[[1]]

gen_dataset_custom <- function(eta, seed, sweeps, update, scale = c("avg", "sum")) {
  scale <- match.arg(scale)
  set.seed(seed)
  coords <- grid_coords(30L)
  n <- nrow(coords)
  p <- 20L
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("beta", seq_len(p))
  beta <- c(1, 2, 3, -4, -5, rep(0, p - 5L))
  y <- rbinom(n, 1L, 0.5)
  xb <- as.numeric(X %*% beta)
  A <- adjacency_matrix(coords)
  deg <- pmax(Matrix::rowSums(A), 1)
  if (update == "parallel") {
    for (s in seq_len(sweeps)) {
      cvec <- as.numeric(A %*% y)
      spatial <- if (scale == "avg") cvec / deg else cvec
      y <- rbinom(n, 1L, plogis(xb + eta * spatial))
    }
  } else if (update == "checkerboard") {
    parity <- (as.integer(round(coords[, 1])) + as.integer(round(coords[, 2]))) %% 2L
    blocks <- list(which(parity == 0L), which(parity == 1L))
    for (s in seq_len(sweeps)) {
      for (idx in blocks) {
        cvec <- as.numeric(A %*% y)
        spatial <- if (scale == "avg") cvec / deg else cvec
        y[idx] <- rbinom(length(idx), 1L, plogis(xb[idx] + eta * spatial[idx]))
      }
    }
  } else {
    stop("unsupported update")
  }
  list(X = X, coords = coords, y = y, beta = beta, eta = eta)
}

fit_eta_glm <- function(d, scale = c("avg", "sum")) {
  scale <- match.arg(scale)
  ns <- neighbor_sum_count(d$coords, d$y)
  spatial <- if (scale == "avg") ns$avg else ns$sum
  dat <- data.frame(y = d$y, d$X, spatial = spatial)
  fit <- suppressWarnings(glm(y ~ . - 1, data = dat, family = binomial()))
  unname(coef(fit)["spatial"])
}

rows <- list()
cat("== Eta generation grid diagnostic ==\n")
cat(sprintf("replicates=%d | eta=%s | sweeps=%s | updates=%s\n",
            n_rep, paste(eta_values, collapse = ","), paste(sweeps_grid, collapse = ","),
            paste(updates, collapse = ",")))

for (eta in eta_values) {
  for (sweeps in sweeps_grid) {
    for (update in updates) {
      for (gen_scale in c("avg", "sum")) {
        fit_scale <- gen_scale
        for (rep_id in seq_len(n_rep)) {
          seed <- 400000L + as.integer(round(eta * 1000)) + sweeps + rep_id
          d <- gen_dataset_custom(eta, seed, sweeps, update, scale = gen_scale)
          eh <- fit_eta_glm(d, scale = fit_scale)
          rows[[length(rows) + 1L]] <- data.frame(
            eta_true = eta,
            sweeps = sweeps,
            update = update,
            gen_scale = gen_scale,
            fit_scale = fit_scale,
            rep = rep_id,
            eta_hat = eh,
            eta_abs_bias = abs(eh - eta),
            disease_rate = mean(d$y),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

raw <- do.call(rbind, rows)
agg <- aggregate(cbind(eta_hat, eta_abs_bias, disease_rate) ~ eta_true + sweeps + update + gen_scale + fit_scale,
                 raw, function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))

raw_path <- file.path(outdir, sprintf("eta_generation_grid_raw_n%s.csv", n_rep))
agg_path <- file.path(outdir, sprintf("eta_generation_grid_aggregate_n%s.csv", n_rep))
write.csv(raw, raw_path, row.names = FALSE)
write.csv(agg, agg_path, row.names = FALSE)

cat("Wrote", raw_path, "\n")
cat("Wrote", agg_path, "\n\n")
print(agg, row.names = FALSE)
