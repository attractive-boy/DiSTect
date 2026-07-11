#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
source(file.path(root, "repro/paper_reproduce/paper_models.R"))

outdir <- file.path(root, "repro/paper_reproduce/output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_rep <- as.integer(Sys.getenv("N_REP", "200"))
eta_values <- as.numeric(strsplit(Sys.getenv("ETA_VALUES", "0.4,1.6,2.8"), ",")[[1]])

fit_glm_eta <- function(d, covariate = c("avg", "sum")) {
  covariate <- match.arg(covariate)
  ns <- neighbor_sum_count(d$coords, d$y)
  cvec <- if (covariate == "avg") ns$avg else ns$sum
  dat <- data.frame(y = d$y, d$X, spatial = cvec)
  fit <- suppressWarnings(glm(y ~ . - 1, data = dat, family = binomial()))
  co <- coef(summary(fit))
  est <- coef(fit)
  list(
    eta = unname(est["spatial"]),
    eta_se = if ("spatial" %in% rownames(co)) unname(co["spatial", "Std. Error"]) else NA_real_,
    beta = unname(est[paste0("beta", seq_len(ncol(d$X)))])
  )
}

rows <- list()
cat("== Eta bias diagnostic: GLM conditional pseudolikelihood ==\n")
cat(sprintf("replicates=%d | eta=%s\n", n_rep, paste(eta_values, collapse = ",")))

for (eta in eta_values) {
  for (rep_id in seq_len(n_rep)) {
    seed <- 100000L + as.integer(round(eta * 1000)) + rep_id
    d <- sim1_dataset(eta = eta, seed = seed)
    for (covariate in c("avg", "sum")) {
      f <- fit_glm_eta(d, covariate = covariate)
      rows[[length(rows) + 1L]] <- data.frame(
        eta_true = eta,
        rep = rep_id,
        covariate = covariate,
        eta_hat = f$eta,
        eta_abs_bias = abs(f$eta - eta),
        beta1_hat = f$beta[1],
        beta2_hat = f$beta[2],
        beta3_hat = f$beta[3],
        beta4_hat = f$beta[4],
        beta5_hat = f$beta[5],
        stringsAsFactors = FALSE
      )
    }
  }
}

raw <- do.call(rbind, rows)
agg <- aggregate(cbind(eta_hat, eta_abs_bias, beta1_hat, beta2_hat, beta3_hat,
                       beta4_hat, beta5_hat) ~ eta_true + covariate, raw,
                 function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))

raw_path <- file.path(outdir, sprintf("eta_bias_glm_raw_n%s.csv", n_rep))
agg_path <- file.path(outdir, sprintf("eta_bias_glm_aggregate_n%s.csv", n_rep))
write.csv(raw, raw_path, row.names = FALSE)
write.csv(agg, agg_path, row.names = FALSE)

cat("Wrote", raw_path, "\n")
cat("Wrote", agg_path, "\n\n")
print(agg, row.names = FALSE)
