#!/usr/bin/env Rscript

# Test the hypothesis that Tables C1-C3 report Monte Carlo standard errors of
# the 200-replicate averages (SE / sqrt(200)) instead of estimator-level SEs.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
outdir <- file.path(root, "repro/paper_reproduce/output")

n_rep <- 200L
paper <- data.frame(
  eta = rep(c(0.4, 1.6, 2.8), each = 2L),
  method = rep(c("NUTS", "ADVI"), 3L),
  published_SEE = c(0.021, 0.044, 0.028, 0.069, 0.047, 0.057),
  published_SEM = c(0.022, 0.044, 0.030, 0.070, 0.045, 0.056),
  stringsAsFactors = FALSE
)

paper$implied_estimator_SEE <- paper$published_SEE * sqrt(n_rep)
paper$implied_estimator_SEM <- paper$published_SEM * sqrt(n_rep)

# Independent local references. GLM is included for all eta settings because it
# bypasses Stan/ADVI. NUTS references are available for eta=0.4 only.
glm <- read.csv(file.path(outdir, "eta_bias_glm_aggregate_n200.csv"))
glm_ref <- setNames(glm$eta_hat.sd, glm$eta_true)

advi <- read.csv(file.path(outdir, "sim1_aggregate_n200.csv"))
advi <- subset(advi, method == "Proposed Method (ADVI)" & param == "eta")
advi_see <- setNames(advi$avgSEE, advi$eta_setting)
advi_sem <- setNames(advi$avgSEM, advi$eta_setting)

nuts <- read.csv(file.path(outdir, "nuts_sim1_tuned.csv"))
nuts <- subset(nuts, adapt_delta == 0.99)

paper$local_GLM_empirical_SD <- unname(glm_ref[as.character(paper$eta)])
paper$local_method_SEE <- NA_real_
paper$local_method_SEM <- NA_real_
for (i in seq_len(nrow(paper))) {
  key <- as.character(paper$eta[i])
  if (paper$method[i] == "ADVI") {
    paper$local_method_SEE[i] <- advi_see[key]
    paper$local_method_SEM[i] <- advi_sem[key]
  } else if (paper$eta[i] == 0.4) {
    paper$local_method_SEE[i] <- sd(nuts$eta_hat)
    paper$local_method_SEM[i] <- mean(nuts$eta_sd)
  }
}

paper$implied_to_local_SEE <- paper$implied_estimator_SEE / paper$local_method_SEE
paper$implied_to_local_SEM <- paper$implied_estimator_SEM / paper$local_method_SEM

num <- vapply(paper, is.numeric, logical(1))
paper[num] <- lapply(paper[num], function(x) round(x, 4))

out <- file.path(outdir, "uncertainty_scaling_audit.csv")
write.csv(paper, out, row.names = FALSE)
cat("Wrote", out, "\n\n")
print(paper, row.names = FALSE)

cat("\nKey check:\n")
cat("Table C1 NUTS eta SEE: 0.021 * sqrt(200) =",
    round(0.021 * sqrt(200), 3), "\n")
cat("Local NUTS eta empirical SD =", round(sd(nuts$eta_hat), 3), "\n")
cat("Table C1 NUTS eta SEM: 0.022 * sqrt(200) =",
    round(0.022 * sqrt(200), 3), "\n")
cat("Local NUTS mean posterior SD =", round(mean(nuts$eta_sd), 3), "\n")
