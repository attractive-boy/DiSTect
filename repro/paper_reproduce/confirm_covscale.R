#!/usr/bin/env Rscript

# Confirmatory test of the "missing covariate scale" explanation.
#
# Two independent lines (CRLB sweep + 8-dataset NUTS) agree that the literal
# design (X ~ N(0,1), beta up to +-5, LP sd ~7.4) forces eta SE ~0.07-0.24, far
# looser than the paper's 0.021-0.044. The CRLB sweep predicts xsd ~0.25-0.5
# would bring eta SE into the paper's range. If the paper (implicitly)
# standardized the covariate / linear-predictor scale, reproducing at a milder
# xsd should recover Table C1's four metrics (avgBias/avgSEE/avgSEM/avgCI) for
# BOTH proposed and naive ADVI -- not just eta.
#
# This reuses the paper ADVI fitters unchanged; only the covariate SD in data
# generation is scaled. Everything else (beta truth, eta grid, prior, 2000-sweep
# Gibbs, fullrank ADVI) is identical to the main reproduction.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
source(file.path(root, "repro/paper_reproduce/paper_models.R"))
outdir <- file.path(root, "repro/paper_reproduce/output")

n_rep    <- as.integer(Sys.getenv("N_REP", "50"))
eta_vals <- as.numeric(strsplit(Sys.getenv("ETA_VALUES", "0.4"), ",")[[1]])
xsd_vals <- as.numeric(strsplit(Sys.getenv("XSD_VALUES", "1.0,0.5"), ",")[[1]])

# Sim1 dataset with covariate SD = xsd (beta truth unchanged).
sim1_xsd <- function(side = 30L, p = 20L, eta = 0.4, xsd = 1.0, seed = 1L) {
  set.seed(seed)
  coords <- grid_coords(side); n <- nrow(coords)
  X <- matrix(rnorm(n * p, sd = xsd), n, p)
  colnames(X) <- paste0("beta", seq_len(p))
  beta <- c(1, 2, 3, -4, -5, rep(0, p - 5L))
  y <- gibbs_autologistic(X, coords, beta, eta, sweeps = 2000L)
  list(X = X, coords = coords, y = y, beta = beta, eta = eta)
}

rows <- list()
cat(sprintf("== Confirm covariate-scale explanation | reps=%d eta=%s xsd=%s ==\n\n",
            n_rep, paste(eta_vals, collapse=","), paste(xsd_vals, collapse=",")))
t0 <- Sys.time()
for (xsd in xsd_vals) {
  for (eta in eta_vals) {
    for (rep_id in seq_len(n_rep)) {
      seed <- 800000L + as.integer(round(xsd*1000))*1000L + as.integer(round(eta*1000)) + rep_id
      d <- sim1_xsd(eta = eta, xsd = xsd, seed = seed)
      # proposed ADVI
      fit_p <- tryCatch(fit_paper_single(d$y, d$X, d$coords, method = "VI"), error=function(e) NULL)
      if (!is.null(fit_p)) {
        tp <- fit_summary_table(fit_p, d$beta, truth_eta = eta); tp$method <- "Proposed (ADVI)"
        tp$xsd <- xsd; tp$eta_setting <- eta; tp$rep <- rep_id; rows[[length(rows)+1L]] <- tp
      }
      # naive ADVI
      fit_n <- tryCatch(fit_paper_naive(d$y, d$X, method = "VI"), error=function(e) NULL)
      if (!is.null(fit_n)) {
        tn <- fit_summary_table(fit_n, d$beta, truth_eta = NULL); tn$method <- "Naive (ADVI)"
        tn$xsd <- xsd; tn$eta_setting <- eta; tn$rep <- rep_id; rows[[length(rows)+1L]] <- tn
      }
    }
    cat(sprintf("  xsd=%.2f eta=%.1f done (%.1f min)\n", xsd, eta,
                as.numeric(difftime(Sys.time(), t0, units="mins"))))
  }
}

raw <- do.call(rbind, rows)
raw_path <- file.path(outdir, sprintf("confirm_covscale_raw_n%s.csv", n_rep))
write.csv(raw, raw_path, row.names = FALSE)

# aggregate to Table-C1 style metrics
agg_fun <- function(d) data.frame(
  xsd=d$xsd[1], eta_setting=d$eta_setting[1], method=d$method[1], param=d$param[1],
  avgBias=round(mean(d$abs_bias),3), avgSEE=round(sd(d$mean),3),
  avgSEM=round(mean(d$sd),3), avgCI=round(mean(d$covered),3), n=nrow(d))
g <- split(raw, list(raw$xsd, raw$eta_setting, raw$method, raw$param), drop=TRUE)
agg <- do.call(rbind, lapply(g, agg_fun)); rownames(agg) <- NULL
agg_path <- file.path(outdir, sprintf("confirm_covscale_aggregate_n%s.csv", n_rep))
write.csv(agg, agg_path, row.names = FALSE)

cat("\nWrote", raw_path, "\nWrote", agg_path, "\n\n")
cat("Paper Table C1 (eta=0.4) targets:\n")
cat("  Proposed ADVI: b1 bias .332 SEE .055 ; eta bias .087 SEE .044 CI .95\n")
cat("  Naive ADVI:    b1 bias .401 SEE .019 ; b5 bias 3.296\n\n")
show <- agg[agg$param %in% c("beta1","beta5","eta"),
            c("xsd","method","param","avgBias","avgSEE","avgSEM","avgCI")]
show <- show[order(show$xsd, show$method, show$param),]
print(show, row.names = FALSE)
