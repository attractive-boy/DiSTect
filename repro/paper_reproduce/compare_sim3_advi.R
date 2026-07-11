#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
outdir <- file.path(root, "repro/paper_reproduce/output")
agg <- read.csv(file.path(outdir, "sim3_aggregate_n200.csv"))

params <- c("beta1", "beta2", "beta3", "beta4", "beta5", "eta")

target_values <- list(
  "10 Missing Spots" = list(
    avgBias = c(0.335, 0.449, 0.465, 0.536, 0.552, 0.089),
    avgSEE = c(0.043, 0.056, 0.048, 0.039, 0.055, 0.011),
    avgSEM = c(0.043, 0.057, 0.048, 0.039, 0.056, 0.012),
    avgCI = c(0.950, 0.950, 0.955, 0.965, 0.935, 0.955)
  ),
  "30 Missing Spots" = list(
    avgBias = c(0.535, 0.796, 0.664, 0.816, 0.787, 0.156),
    avgSEE = c(0.028, 0.028, 0.038, 0.048, 0.047, 0.029),
    avgSEM = c(0.028, 0.029, 0.038, 0.047, 0.047, 0.029),
    avgCI = c(0.970, 0.955, 0.950, 0.965, 0.955, 0.950)
  ),
  "gamma=(-6,1,4)" = list(
    avgBias = c(0.398, 0.349, 0.489, 0.501, 0.471, 0.062),
    avgSEE = c(0.041, 0.037, 0.048, 0.050, 0.059, 0.052),
    avgSEM = c(0.041, 0.037, 0.048, 0.051, 0.056, 0.052),
    avgCI = c(0.945, 0.940, 0.940, 0.955, 0.955, 0.955)
  ),
  "gamma=(-5,1,1.6)" = list(
    avgBias = c(0.438, 0.451, 0.511, 0.519, 0.507, 0.066),
    avgSEE = c(0.046, 0.048, 0.053, 0.042, 0.055, 0.055),
    avgSEM = c(0.045, 0.047, 0.053, 0.045, 0.055, 0.055),
    avgCI = c(0.950, 0.955, 0.960, 0.955, 0.925, 0.950)
  )
)

target <- do.call(rbind, lapply(names(target_values), function(setting) {
  vals <- target_values[[setting]]
  data.frame(
    setting = setting,
    method = "Missing first-pass (ADVI)",
    param = params,
    target_avgBias = vals$avgBias,
    target_avgSEE = vals$avgSEE,
    target_avgSEM = vals$avgSEM,
    target_avgCI = vals$avgCI,
    stringsAsFactors = FALSE
  )
}))

cmp <- merge(target, agg, by = c("setting", "method", "param"), all.x = TRUE)
cmp$delta_avgBias <- cmp$avgBias - cmp$target_avgBias
cmp$delta_avgSEE <- cmp$avgSEE - cmp$target_avgSEE
cmp$delta_avgSEM <- cmp$avgSEM - cmp$target_avgSEM
cmp$delta_avgCI <- cmp$avgCI - cmp$target_avgCI
cmp <- cmp[order(match(cmp$setting, names(target_values)), match(cmp$param, params)), ]

out <- file.path(outdir, "sim3_advi_key_comparison_n200.csv")
write.csv(cmp, out, row.names = FALSE)

cat("Wrote", out, "\n\n")
print(cmp[, c("setting", "param", "target_avgBias", "avgBias", "delta_avgBias",
              "target_avgCI", "avgCI", "delta_avgCI", "meanMissing")], row.names = FALSE)
