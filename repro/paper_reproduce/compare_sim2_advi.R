#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
outdir <- file.path(root, "repro/paper_reproduce/output")
agg <- read.csv(file.path(outdir, "sim2_aggregate_n200.csv"))

settings <- data.frame(
  sigma2 = c(0.1, 0.1, 0.4, 0.4),
  rho = c(0.1, 0.4, 0.1, 0.4)
)
params <- c("beta1", "beta2", "beta3", "beta4", "beta5", "eta")

target_values <- list(
  "0.1_0.1" = list(
    avgBias = c(0.303, 0.331, 0.398, 0.494, 0.469, 0.076),
    avgSEE = c(0.046, 0.032, 0.048, 0.057, 0.043, 0.025),
    avgSEM = c(0.047, 0.032, 0.048, 0.059, 0.043, 0.025),
    avgCI = c(0.955, 0.940, 0.955, 0.965, 0.950, 0.945)
  ),
  "0.1_0.4" = list(
    avgBias = c(0.362, 0.481, 0.452, 0.568, 0.524, 0.089),
    avgSEE = c(0.045, 0.054, 0.042, 0.076, 0.054, 0.031),
    avgSEM = c(0.043, 0.054, 0.043, 0.074, 0.054, 0.031),
    avgCI = c(0.940, 0.950, 0.945, 0.955, 0.935, 0.940)
  ),
  "0.4_0.1" = list(
    avgBias = c(0.288, 0.351, 0.424, 0.412, 0.459, 0.083),
    avgSEE = c(0.034, 0.045, 0.046, 0.057, 0.063, 0.022),
    avgSEM = c(0.034, 0.045, 0.045, 0.054, 0.063, 0.022),
    avgCI = c(0.950, 0.940, 0.940, 0.950, 0.955, 0.935)
  ),
  "0.4_0.4" = list(
    avgBias = c(0.371, 0.323, 0.511, 0.493, 0.551, 0.081),
    avgSEE = c(0.032, 0.059, 0.048, 0.049, 0.063, 0.029),
    avgSEM = c(0.033, 0.059, 0.047, 0.045, 0.059, 0.026),
    avgCI = c(0.945, 0.955, 0.930, 0.950, 0.945, 0.950)
  )
)

target <- do.call(rbind, lapply(seq_len(nrow(settings)), function(i) {
  key <- sprintf("%.1f_%.1f", settings$sigma2[i], settings$rho[i])
  vals <- target_values[[key]]
  data.frame(
    sigma2 = settings$sigma2[i],
    rho = settings$rho[i],
    method = "Proposed Method (ADVI)",
    param = params,
    target_avgBias = vals$avgBias,
    target_avgSEE = vals$avgSEE,
    target_avgSEM = vals$avgSEM,
    target_avgCI = vals$avgCI,
    stringsAsFactors = FALSE
  )
}))

cmp <- merge(target, agg, by = c("sigma2", "rho", "method", "param"), all.x = TRUE)
cmp$delta_avgBias <- cmp$avgBias - cmp$target_avgBias
cmp$delta_avgSEE <- cmp$avgSEE - cmp$target_avgSEE
cmp$delta_avgSEM <- cmp$avgSEM - cmp$target_avgSEM
cmp$delta_avgCI <- cmp$avgCI - cmp$target_avgCI
cmp <- cmp[order(cmp$sigma2, cmp$rho, match(cmp$param, params)), ]

out <- file.path(outdir, "sim2_advi_key_comparison_n200.csv")
write.csv(cmp, out, row.names = FALSE)

cat("Wrote", out, "\n\n")
print(cmp[, c("sigma2", "rho", "param", "target_avgBias", "avgBias",
              "delta_avgBias", "target_avgCI", "avgCI", "delta_avgCI")], row.names = FALSE)
