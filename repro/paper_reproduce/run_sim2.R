#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
source(file.path(root, "repro/paper_reproduce/paper_models.R"))

outdir <- file.path(root, "repro/paper_reproduce/output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_rep <- as.integer(Sys.getenv("N_REP", "1"))
settings <- data.frame(
  sigma2 = c(0.1, 0.1, 0.4, 0.4),
  rho = c(0.1, 0.4, 0.1, 0.4)
)

cat("== btaf530 Supplement A2.2 / Simulation 2 ==\n")
cat(sprintf("replicates=%d | six slices | ADVI proposed multiple-slice model\n", n_rep))

rows <- list()
timings <- data.frame()

for (k in seq_len(nrow(settings))) {
  sigma2 <- settings$sigma2[k]
  rho <- settings$rho[k]
  setting <- sprintf("sigma2=%.1f,rho=%.1f", sigma2, rho)
  for (rep_id in seq_len(n_rep)) {
    seed <- 200000L + k * 1000L + rep_id
    cat(sprintf("\n-- %s replicate %d/%d (seed=%d) --\n", setting, rep_id, n_rep, seed))
    d <- sim2_dataset(sigma2 = sigma2, rho = rho, seed = seed)
    t <- system.time(
      fit <- fit_paper_multiple(d$y, d$X, d$coords, d$label, method = "VI")
    )["elapsed"]
    tab <- fit_summary_table(fit, d$beta, d$eta)
    tab$method <- "Proposed Method (ADVI)"
    tab$setting <- setting
    tab$sigma2 <- sigma2
    tab$rho <- rho
    tab$rep <- rep_id
    rows[[length(rows) + 1L]] <- tab
    timings <- rbind(timings, data.frame(setting = setting, sigma2 = sigma2, rho = rho,
                                         rep = rep_id, method = "Proposed Method (ADVI)",
                                         seconds = as.numeric(t)))
    cat(sprintf("   proposed multiple-slice ADVI %.1fs\n", t))
  }
}

raw <- do.call(rbind, rows)
groups <- split(raw, list(raw$setting, raw$param), drop = TRUE)
agg <- do.call(rbind, lapply(groups, function(d) {
  data.frame(
    method = d$method[1],
    setting = d$setting[1],
    sigma2 = d$sigma2[1],
    rho = d$rho[1],
    param = d$param[1],
    avgBias = mean(d$abs_bias),
    avgSEE = sd(d$mean),
    avgSEM = mean(d$sd),
    avgCI = mean(d$covered),
    n = nrow(d),
    stringsAsFactors = FALSE
  )
}))
rownames(agg) <- NULL
agg <- agg[order(agg$sigma2, agg$rho, match(agg$param, c(paste0("beta", 1:20), "eta"))), ]
raw <- raw[, c("method", "setting", "sigma2", "rho", "rep", "param", "truth", "mean",
               "sd", "q025", "q975", "bias", "abs_bias", "covered")]

raw_path <- file.path(outdir, sprintf("sim2_raw_n%s.csv", n_rep))
agg_path <- file.path(outdir, sprintf("sim2_aggregate_n%s.csv", n_rep))
time_path <- file.path(outdir, sprintf("sim2_timings_n%s.csv", n_rep))
write.csv(raw, raw_path, row.names = FALSE)
write.csv(agg, agg_path, row.names = FALSE)
write.csv(timings, time_path, row.names = FALSE)

cat("\n== Simulation 2 outputs ==\n")
cat(raw_path, "\n")
cat(agg_path, "\n")
cat(time_path, "\n")
cat("\nAggregate preview:\n")
print(head(agg, 30), row.names = FALSE)
