#!/usr/bin/env Rscript

# Follow-up to nuts_sim1_check: the naive NUTS fit produced 67 divergences and
# eta bias 0.417 -- the spike-slab funnel (v0 = 1e-6) is not being sampled
# correctly. Test whether standard funnel remedies (higher adapt_delta, longer
# warmup, deeper trees) eliminate divergences and move eta toward the paper's
# Proposed (NUTS) target of ~0.019. If so, the paper's 23h/200-rep NUTS cost and
# clean numbers are both explained by expensive well-adapted sampling, and the
# local model is confirmed correct.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
source(file.path(root, "repro/paper_reproduce/paper_models.R"))

outdir <- file.path(root, "repro/paper_reproduce/output")

n_rep      <- as.integer(Sys.getenv("N_REP", "2"))
eta        <- as.numeric(Sys.getenv("ETA", "0.4"))
iter       <- as.integer(Sys.getenv("ITER", "2000"))
warmup     <- as.integer(Sys.getenv("WARMUP", "1000"))
chains     <- as.integer(Sys.getenv("CHAINS", "2"))
ad_grid    <- as.numeric(strsplit(Sys.getenv("ADAPT_DELTA", "0.99,0.999"), ",")[[1]])
max_td     <- as.integer(Sys.getenv("MAX_TREEDEPTH", "13"))

model <- get_stan_model("paper_single", paper_single_stan)

rows <- list()
cat(sprintf("== NUTS tuning on Sim1 | eta=%.1f reps=%d iter=%d warmup=%d chains=%d td=%d ==\n",
            eta, n_rep, iter, warmup, chains, max_td))
cat("Paper Proposed(NUTS) eta bias target: 0.019 (SEE 0.021); betaBias ~0.09\n\n")

for (ad in ad_grid) {
  for (rep_id in seq_len(n_rep)) {
    seed <- 700000L + as.integer(round(eta * 1000)) + rep_id
    d <- sim1_dataset(eta = eta, seed = seed)
    c_avg <- neighbor_sum_count(d$coords, d$y)$avg
    data <- list(N = nrow(d$X), P = ncol(d$X), x = d$X, y = as.integer(d$y), c_avg = c_avg)
    t0 <- Sys.time()
    fit <- sampling(model, data = data, iter = iter, warmup = warmup, chains = chains,
                    seed = 128L, refresh = 0,
                    control = list(adapt_delta = ad, max_treedepth = max_td))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    tab <- fit_summary_table(fit, d$beta, truth_eta = eta)
    sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    ndiv <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
    er <- tab[tab$param == "eta", ]
    b  <- tab[tab$param %in% paste0("beta", 1:5), ]
    rows[[length(rows) + 1L]] <- data.frame(
      adapt_delta = ad, rep = rep_id, secs = round(secs, 1), divergences = ndiv,
      eta_hat = round(er$mean, 3), eta_sd = round(er$sd, 3),
      eta_absbias = round(abs(er$mean - eta), 3),
      beta_absbias_mean = round(mean(abs(b$mean - b$truth)), 3),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  ad=%.3f rep=%d | %.0fs | div=%d | eta_hat=%.3f (bias %.3f) | betaBias=%.3f\n",
                ad, rep_id, secs, ndiv, er$mean, abs(er$mean - eta), mean(abs(b$mean - b$truth))))
  }
}

raw <- do.call(rbind, rows)
raw_path <- file.path(outdir, "nuts_sim1_tuned.csv")
write.csv(raw, raw_path, row.names = FALSE)
cat("\nWrote", raw_path, "\n\n")
agg <- aggregate(cbind(divergences, eta_hat, eta_absbias, eta_sd, beta_absbias_mean, secs) ~ adapt_delta,
                 raw, function(x) round(mean(x), 3))
print(agg, row.names = FALSE)
