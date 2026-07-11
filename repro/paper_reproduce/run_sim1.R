#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
source(file.path(root, "repro/paper_reproduce/paper_models.R"))

outdir <- file.path(root, "repro/paper_reproduce/output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_rep <- as.integer(Sys.getenv("N_REP", "1"))
eta_values <- as.numeric(strsplit(Sys.getenv("ETA_VALUES", "0.4,1.6,2.8"), ",")[[1]])
run_nuts <- tolower(Sys.getenv("RUN_NUTS", "false")) %in% c("1", "true", "yes")
nuts_iter <- as.integer(Sys.getenv("NUTS_ITER", "800"))
nuts_warmup <- as.integer(Sys.getenv("NUTS_WARMUP", "400"))

cat("== btaf530 Supplement A2.1 / Simulation 1 ==\n")
cat(sprintf("replicates=%d | eta=%s | ADVI=yes | NUTS=%s\n",
            n_rep, paste(eta_values, collapse = ","), run_nuts))
cat("Formula: logit(mu_i) = beta'X_i + eta * mean_neighbor_y_i\n")

rows <- list()
timings <- data.frame()

for (eta in eta_values) {
  for (rep_id in seq_len(n_rep)) {
    seed <- 100000L + as.integer(round(eta * 1000)) + rep_id
    cat(sprintf("\n-- eta %.1f replicate %d/%d (seed=%d) --\n", eta, rep_id, n_rep, seed))
    d <- sim1_dataset(eta = eta, seed = seed)

    t <- system.time(fit_naive_vi <- fit_paper_naive(d$y, d$X, method = "VI"))["elapsed"]
    tab <- fit_summary_table(fit_naive_vi, d$beta)
    tab$method <- "Naive Method (ADVI)"
    tab$eta_setting <- eta
    tab$rep <- rep_id
    rows[[length(rows) + 1L]] <- tab
    timings <- rbind(timings, data.frame(eta_setting = eta, rep = rep_id,
                                         method = "Naive Method (ADVI)",
                                         seconds = as.numeric(t)))
    cat(sprintf("   naive ADVI %.1fs\n", t))

    t <- system.time(fit_prop_vi <- fit_paper_single(d$y, d$X, d$coords, method = "VI"))["elapsed"]
    tab <- fit_summary_table(fit_prop_vi, d$beta, d$eta)
    tab$method <- "Proposed Method (ADVI)"
    tab$eta_setting <- eta
    tab$rep <- rep_id
    rows[[length(rows) + 1L]] <- tab
    timings <- rbind(timings, data.frame(eta_setting = eta, rep = rep_id,
                                         method = "Proposed Method (ADVI)",
                                         seconds = as.numeric(t)))
    cat(sprintf("   proposed ADVI %.1fs\n", t))

    if (run_nuts) {
      t <- system.time(fit_naive_nuts <- fit_paper_naive(
        d$y, d$X, method = "NUTS", iter = nuts_iter, warmup = nuts_warmup
      ))["elapsed"]
      tab <- fit_summary_table(fit_naive_nuts, d$beta)
      tab$method <- "Naive Method (NUTS)"
      tab$eta_setting <- eta
      tab$rep <- rep_id
      rows[[length(rows) + 1L]] <- tab
      timings <- rbind(timings, data.frame(eta_setting = eta, rep = rep_id,
                                           method = "Naive Method (NUTS)",
                                           seconds = as.numeric(t)))
      cat(sprintf("   naive NUTS %.1fs\n", t))

      t <- system.time(fit_prop_nuts <- fit_paper_single(
        d$y, d$X, d$coords, method = "NUTS", iter = nuts_iter, warmup = nuts_warmup
      ))["elapsed"]
      tab <- fit_summary_table(fit_prop_nuts, d$beta, d$eta)
      tab$method <- "Proposed Method (NUTS)"
      tab$eta_setting <- eta
      tab$rep <- rep_id
      rows[[length(rows) + 1L]] <- tab
      timings <- rbind(timings, data.frame(eta_setting = eta, rep = rep_id,
                                           method = "Proposed Method (NUTS)",
                                           seconds = as.numeric(t)))
      cat(sprintf("   proposed NUTS %.1fs\n", t))
    }
  }
}

raw <- do.call(rbind, rows)
agg <- aggregate_sim_table(rows)
raw <- raw[, c("method", "eta_setting", "rep", "param", "truth", "mean", "sd",
               "q025", "q975", "bias", "abs_bias", "covered")]
agg <- agg[order(agg$eta_setting, agg$method, match(agg$param, c(paste0("beta", 1:20), "eta"))), ]

raw_path <- file.path(outdir, sprintf("sim1_raw_n%s.csv", n_rep))
agg_path <- file.path(outdir, sprintf("sim1_aggregate_n%s.csv", n_rep))
time_path <- file.path(outdir, sprintf("sim1_timings_n%s.csv", n_rep))
write.csv(raw, raw_path, row.names = FALSE)
write.csv(agg, agg_path, row.names = FALSE)
write.csv(timings, time_path, row.names = FALSE)

cat("\n== Simulation 1 outputs ==\n")
cat(raw_path, "\n")
cat(agg_path, "\n")
cat(time_path, "\n")
cat("\nAggregate preview:\n")
print(head(agg, 30), row.names = FALSE)
