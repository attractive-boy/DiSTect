#!/usr/bin/env Rscript

# Decisive check: does the LOCAL Simulation-1 model, fit with NUTS (the paper's
# gold-standard sampler), reproduce the paper's Proposed (NUTS) column?
#
# Paper Table C1 (eta = 0.4), Proposed (NUTS): eta avgBias 0.019, SEE 0.021;
#   signal-beta avgBias ~ 0.088-0.11 (flat).
# Paper Table C1 (eta = 0.4), Proposed (ADVI): eta avgBias 0.087, SEE 0.044.
# Local fullrank-ADVI reproduction: eta avgBias ~0.231, SEE ~0.245 (5x too wide).
#
# If local NUTS lands near the paper NUTS column, the local model is correct and
# the entire ADVI-column mismatch is an ADVI convergence/implementation artifact
# (fullrank vs meanfield, tolerance, spike-funnel instability) rather than a
# data/model/prior mismatch. This would flip the "cannot reproduce" conclusion.
#
# Runtime: paper reports ~23h per 200-rep NUTS study => ~7 min/fit. The stiff
# spike-slab (v0 = 1e-6) may push this higher and produce divergences. This
# script records wall-time and divergence counts per fit.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
source(file.path(root, "repro/paper_reproduce/paper_models.R"))

outdir <- file.path(root, "repro/paper_reproduce/output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_rep     <- as.integer(Sys.getenv("N_REP", "3"))
eta_values<- as.numeric(strsplit(Sys.getenv("ETA_VALUES", "0.4"), ",")[[1]])
iter      <- as.integer(Sys.getenv("ITER", "1500"))
warmup    <- as.integer(Sys.getenv("WARMUP", "750"))
chains    <- as.integer(Sys.getenv("CHAINS", "2"))

paper_nuts <- list("0.4" = list(eta = 0.019, beta = c(0.089, 0.092, 0.095, 0.11, 0.088)),
                   "1.6" = list(eta = NA, beta = rep(NA, 5)),
                   "2.8" = list(eta = NA, beta = rep(NA, 5)))

rows <- list()
cat(sprintf("== Local NUTS check on Sim1 | reps=%d eta=%s iter=%d warmup=%d chains=%d ==\n\n",
            n_rep, paste(eta_values, collapse = ","), iter, warmup, chains))

for (eta in eta_values) {
  for (rep_id in seq_len(n_rep)) {
    seed <- 700000L + as.integer(round(eta * 1000)) + rep_id
    d <- sim1_dataset(eta = eta, seed = seed)
    t0 <- Sys.time()
    fit <- fit_paper_single(d$y, d$X, d$coords, method = "NUTS",
                            iter = iter, warmup = warmup, chains = chains, seed = 128L)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    tab <- fit_summary_table(fit, d$beta, truth_eta = eta)
    # divergences
    sp <- tryCatch(rstan::get_sampler_params(fit, inc_warmup = FALSE), error = function(e) NULL)
    ndiv <- if (!is.null(sp)) sum(sapply(sp, function(x) sum(x[, "divergent__"]))) else NA
    eta_row <- tab[tab$param == "eta", ]
    b <- tab[tab$param %in% paste0("beta", 1:5), ]
    rows[[length(rows) + 1L]] <- data.frame(
      eta_true = eta, rep = rep_id, secs = round(secs, 1), divergences = ndiv,
      eta_hat = round(eta_row$mean, 3), eta_sd = round(eta_row$sd, 3),
      eta_absbias = round(abs(eta_row$mean - eta), 3),
      b1 = round(b$mean[1], 2), b2 = round(b$mean[2], 2), b3 = round(b$mean[3], 2),
      b4 = round(b$mean[4], 2), b5 = round(b$mean[5], 2),
      beta_absbias_mean = round(mean(abs(b$mean - b$truth)), 3),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  eta=%.1f rep=%d | %.0fs | div=%s | eta_hat=%.3f (bias %.3f) | betaBias=%.3f\n",
                eta, rep_id, secs, as.character(ndiv), eta_row$mean,
                abs(eta_row$mean - eta), mean(abs(b$mean - b$truth))))
  }
}

raw <- do.call(rbind, rows)
raw_path <- file.path(outdir, sprintf("nuts_sim1_check_n%s.csv", n_rep))
write.csv(raw, raw_path, row.names = FALSE)

cat("\nWrote", raw_path, "\n\n")
cat("Reference (paper Table C1, eta=0.4):\n")
cat("  Proposed NUTS : eta bias 0.019, SEE 0.021 ; signal-beta bias ~0.09 (flat)\n")
cat("  Proposed ADVI : eta bias 0.087, SEE 0.044\n")
cat("  Local fullrank ADVI (n200): eta bias ~0.231, SEE ~0.245\n\n")
agg <- aggregate(cbind(eta_hat, eta_absbias, eta_sd, beta_absbias_mean, secs, divergences) ~ eta_true,
                 raw, function(x) round(mean(x), 3))
print(agg, row.names = FALSE)
