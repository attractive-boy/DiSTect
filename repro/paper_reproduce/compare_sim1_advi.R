#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
outdir <- file.path(root, "repro/paper_reproduce/output")
agg <- read.csv(file.path(outdir, "sim1_aggregate_n200.csv"))

target <- data.frame(
  eta_setting = c(
    rep(0.4, 12), rep(1.6, 12), rep(2.8, 12)
  ),
  method = rep(c(rep("Naive Method (ADVI)", 5), rep("Proposed Method (ADVI)", 6), "Proposed Method (ADVI)"), 3),
  param = rep(c("beta1", "beta2", "beta3", "beta4", "beta5",
                "beta1", "beta2", "beta3", "beta4", "beta5", "eta", "eta"), 3),
  metric_group = rep(c(rep("naive_beta", 5), rep("proposed_beta", 5), "proposed_eta", "proposed_eta_duplicate"), 3),
  target_avgBias = c(
    0.401, 0.648, 1.830, 2.274, 3.296, 0.332, 0.353, 0.451, 0.445, 0.458, 0.087, 0.087,
    0.389, 0.689, 1.540, 2.386, 3.982, 0.353, 0.443, 0.355, 0.452, 0.447, 0.083, 0.083,
    0.423, 0.582, 1.673, 2.833, 3.786, 0.334, 0.467, 0.344, 0.356, 0.556, 0.095, 0.095
  ),
  target_avgSEE = c(
    0.019, 0.025, 0.038, 0.059, 0.064, 0.055, 0.047, 0.065, 0.047, 0.056, 0.044, 0.044,
    0.016, 0.025, 0.034, 0.067, 0.069, 0.043, 0.032, 0.033, 0.046, 0.038, 0.069, 0.069,
    0.017, 0.027, 0.038, 0.063, 0.072, 0.055, 0.047, 0.065, 0.047, 0.056, 0.057, 0.057
  ),
  target_avgSEM = c(
    0.018, 0.025, 0.038, 0.058, 0.066, 0.056, 0.047, 0.066, 0.047, 0.057, 0.044, 0.044,
    0.016, 0.027, 0.034, 0.067, 0.069, 0.043, 0.033, 0.034, 0.047, 0.039, 0.070, 0.070,
    0.017, 0.027, 0.037, 0.063, 0.073, 0.056, 0.047, 0.066, 0.047, 0.057, 0.056, 0.056
  ),
  target_avgCI = c(
    0.950, 0.955, 0.945, 0.960, 0.950, 0.955, 0.955, 0.950, 0.950, 0.945, 0.950, 0.950,
    0.955, 0.950, 0.940, 0.950, 0.955, 0.935, 0.955, 0.945, 0.960, 0.965, 0.950, 0.950,
    0.940, 0.945, 0.955, 0.970, 0.955, 0.945, 0.965, 0.955, 0.940, 0.940, 0.955, 0.955
  )
)

target <- subset(target, metric_group != "proposed_eta_duplicate")
cmp <- merge(target, agg, by = c("eta_setting", "method", "param"), all.x = TRUE)
cmp$delta_avgBias <- cmp$avgBias - cmp$target_avgBias
cmp$delta_avgSEE <- cmp$avgSEE - cmp$target_avgSEE
cmp$delta_avgSEM <- cmp$avgSEM - cmp$target_avgSEM
cmp$delta_avgCI <- cmp$avgCI - cmp$target_avgCI
cmp <- cmp[order(cmp$eta_setting, cmp$method, match(cmp$param, c("beta1","beta2","beta3","beta4","beta5","eta"))), ]

out <- file.path(outdir, "sim1_advi_key_comparison_n200.csv")
write.csv(cmp, out, row.names = FALSE)

cat("Wrote", out, "\n\n")
print(cmp[, c("eta_setting", "method", "param", "target_avgBias", "avgBias",
              "delta_avgBias", "target_avgCI", "avgCI", "delta_avgCI")], row.names = FALSE)
