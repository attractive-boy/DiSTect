#!/usr/bin/env Rscript

# Figures for ETA_UNCERTAINTY_FINDING.md
#   fig_eta_se_floor.png : eta SE vs covariate scale -> floor ~0.11 >> paper
#   fig_mcmle_vs_pl.png  : full-likelihood MCMLE SE == pseudolikelihood SE (ratio ~1)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(dirname(script_file))
od <- file.path(root, "output")

## ---- Figure 1: eta SE floor vs covariate scale ----
a <- read.csv(file.path(od, "covariate_scale_aggregate_n50.csv"))
a <- a[order(a$xsd), ]

png(file.path(root, "fig_eta_se_floor.png"), width = 1500, height = 1050, res = 200)
par(mar = c(4.6, 4.8, 3.2, 1.4))
ymax <- 0.32
plot(a$lp_sd, a$eta_meanSE, type = "b", pch = 19, col = "#1b6ca8", lwd = 2,
     xlab = "linear-predictor SD  (= xsd x 7.42;  literal design = 7.42)",
     ylab = "eta standard error",
     main = "Simulation 1 (eta = 0.4): eta SE floors at ~0.11 at any covariate scale",
     ylim = c(0, ymax), xlim = c(0, 7.6), cex.main = 0.98)
lines(a$lp_sd, a$eta_empSD, type = "b", pch = 1, col = "#1b6ca8", lwd = 1.5, lty = 2)
abline(h = 0.113, col = "#1b6ca8", lty = 3)
text(6.0, 0.128, "CR floor ~0.11", col = "#1b6ca8", cex = 0.85)
# paper reference lines
abline(h = 0.044, col = "#c0392b", lwd = 2); text(6.2, 0.052, "paper ADVI SEE 0.044", col = "#c0392b", cex = 0.82)
abline(h = 0.021, col = "#8e44ad", lwd = 2); text(6.2, 0.010, "paper NUTS SEE 0.021", col = "#8e44ad", cex = 0.82)
# local NUTS empirical SD (n=8) at literal design
points(7.42, 0.243, pch = 17, col = "#e67e22", cex = 1.3); text(6.3, 0.243, "local NUTS SD 0.24", col = "#e67e22", cex = 0.82)
legend("topleft", bty = "n", cex = 0.85,
       legend = c("mean posterior SD (avgSEM ~ CRLB)", "empirical SD across reps (avgSEE)"),
       col = "#1b6ca8", pch = c(19, 1), lty = c(1, 2), lwd = c(2, 1.5))
dev.off()

## ---- Figure 2: full-likelihood MCMLE == pseudolikelihood (no efficiency gain) ----
map <- list(c("011", 0.11, 0.4), c("044", 0.44, 1.6), c("078", 0.78, 2.8), c("16", 1.6, 5.8))
rows <- list()
for (m in map) {
  f <- file.path(od, sprintf("mcmle_vs_pl_eta%s_n2.csv", m[1]))
  if (!file.exists(f)) next
  d <- read.csv(f)
  rows[[length(rows) + 1L]] <- data.frame(eta_avg = as.numeric(m[3]),
    CV_PL = mean(d$CV_PL), CV_MCMLE = mean(d$CV_MCMLE),
    ratio = mean(d$SE_ratio_PL_over_MCMLE))
}
r <- do.call(rbind, rows); r <- r[order(r$eta_avg), ]

png(file.path(root, "fig_mcmle_vs_pl.png"), width = 1500, height = 1050, res = 200)
par(mar = c(4.6, 4.8, 3.2, 1.4))
plot(r$eta_avg, r$CV_PL, type = "b", pch = 19, col = "#1b6ca8", lwd = 2,
     xlab = "eta (avg-scale equivalent; last point is near-critical)",
     ylab = "eta coefficient of variation  (SE / eta)",
     main = "Full-likelihood MCMLE gives no efficiency gain over pseudolikelihood",
     ylim = c(0, max(r$CV_PL) * 1.15), cex.main = 0.98)
lines(r$eta_avg, r$CV_MCMLE, type = "b", pch = 4, col = "#e67e22", lwd = 2, lty = 2)
# paper C1 eta=0.4 points
points(0.4, 0.044 / 0.4, pch = 15, col = "#c0392b", cex = 1.2)
points(0.4, 0.021 / 0.4, pch = 17, col = "#8e44ad", cex = 1.2)
text(0.9, 0.044/0.4, "paper ADVI", col = "#c0392b", cex = 0.8, pos = 4)
text(0.9, 0.021/0.4 - 0.02, "paper NUTS", col = "#8e44ad", cex = 0.8, pos = 4)
legend("topright", bty = "n", cex = 0.85,
       legend = c("pseudolikelihood (GLM)  CV", "full-likelihood MCMLE  CV",
                  sprintf("SE ratio PL/MCMLE = %.2f (all)", mean(r$ratio))),
       col = c("#1b6ca8", "#e67e22", NA), pch = c(19, 4, NA), lty = c(1, 2, NA), lwd = c(2, 2, NA))
dev.off()

cat("Wrote fig_eta_se_floor.png and fig_mcmle_vs_pl.png\n")
print(r, row.names = FALSE)
