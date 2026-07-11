#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
source(file.path(root, "repro/paper_reproduce/paper_models.R"))

outdir <- file.path(root, "repro/paper_reproduce/output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_rep <- as.integer(Sys.getenv("N_REP", "1"))
impute_iter <- as.integer(Sys.getenv("IMPUTE_ITER", "1"))
eta <- as.numeric(Sys.getenv("ETA", "1.6"))

settings <- list(
  list(kind = "ignorable", setting = "10 Missing Spots", n_missing = 10L, gamma = NA),
  list(kind = "ignorable", setting = "30 Missing Spots", n_missing = 30L, gamma = NA),
  list(kind = "nonignorable", setting = "gamma=(-6,1,4)", n_missing = NA, gamma = c(-6, 1, 4)),
  list(kind = "nonignorable", setting = "gamma=(-5,1,1.6)", n_missing = NA, gamma = c(-5, 1, 1.6))
)

cat("== btaf530 Supplement A2.3 / Simulation 3 ==\n")
cat(sprintf("replicates=%d | eta=%.1f | impute_iter=%d | ADVI missing first-pass\n",
            n_rep, eta, impute_iter))

rows <- list()
timings <- data.frame()
miss_counts <- data.frame()

for (cfg in settings) {
  for (rep_id in seq_len(n_rep)) {
    seed <- 300000L + match(cfg$setting, vapply(settings, `[[`, "", "setting")) * 1000L + rep_id
    cat(sprintf("\n-- %s replicate %d/%d (seed=%d) --\n", cfg$setting, rep_id, n_rep, seed))
    d <- sim1_dataset(eta = eta, seed = seed)
    if (cfg$kind == "ignorable") {
      observed <- mask_ignorable(d$coords, cfg$n_missing, seed = seed + 17L)
    } else {
      observed <- gibbs_missing_indicator(d$y, d$coords, cfg$gamma, seed = seed + 17L)
    }
    n_miss <- sum(!observed)
    miss_counts <- rbind(miss_counts, data.frame(setting = cfg$setting, kind = cfg$kind,
                                                 rep = rep_id, missing = n_miss))
    cat(sprintf("   missing spots: %d\n", n_miss))

    t <- system.time(
      fit <- fit_missing_firstpass(d$y, d$X, d$coords, observed,
                                   impute_iter = impute_iter, method = "VI")
    )["elapsed"]
    tab <- fit_summary_table(fit, d$beta, d$eta)
    tab$method <- "Missing first-pass (ADVI)"
    tab$setting <- cfg$setting
    tab$kind <- cfg$kind
    tab$rep <- rep_id
    tab$missing <- n_miss
    rows[[length(rows) + 1L]] <- tab
    timings <- rbind(timings, data.frame(setting = cfg$setting, kind = cfg$kind,
                                         rep = rep_id, method = "Missing first-pass (ADVI)",
                                         seconds = as.numeric(t)))
    cat(sprintf("   missing ADVI %.1fs\n", t))
  }
}

raw <- do.call(rbind, rows)
groups <- split(raw, list(raw$setting, raw$param), drop = TRUE)
agg <- do.call(rbind, lapply(groups, function(d) {
  data.frame(
    method = d$method[1],
    setting = d$setting[1],
    kind = d$kind[1],
    param = d$param[1],
    avgBias = mean(d$abs_bias),
    avgSEE = sd(d$mean),
    avgSEM = mean(d$sd),
    avgCI = mean(d$covered),
    meanMissing = mean(d$missing),
    n = nrow(d),
    stringsAsFactors = FALSE
  )
}))
rownames(agg) <- NULL
setting_order <- vapply(settings, `[[`, "", "setting")
agg <- agg[order(match(agg$setting, setting_order),
                 match(agg$param, c(paste0("beta", 1:20), "eta"))), ]
raw <- raw[, c("method", "setting", "kind", "rep", "missing", "param", "truth", "mean",
               "sd", "q025", "q975", "bias", "abs_bias", "covered")]

raw_path <- file.path(outdir, sprintf("sim3_raw_n%s.csv", n_rep))
agg_path <- file.path(outdir, sprintf("sim3_aggregate_n%s.csv", n_rep))
time_path <- file.path(outdir, sprintf("sim3_timings_n%s.csv", n_rep))
miss_path <- file.path(outdir, sprintf("sim3_missing_counts_n%s.csv", n_rep))
write.csv(raw, raw_path, row.names = FALSE)
write.csv(agg, agg_path, row.names = FALSE)
write.csv(timings, time_path, row.names = FALSE)
write.csv(miss_counts, miss_path, row.names = FALSE)

cat("\n== Simulation 3 outputs ==\n")
cat(raw_path, "\n")
cat(agg_path, "\n")
cat(time_path, "\n")
cat(miss_path, "\n")
cat("\nMissing counts:\n")
print(aggregate(missing ~ setting, miss_counts, function(x) c(mean = mean(x), min = min(x), max = max(x))))
cat("\nAggregate preview:\n")
print(head(agg, 40), row.names = FALSE)
