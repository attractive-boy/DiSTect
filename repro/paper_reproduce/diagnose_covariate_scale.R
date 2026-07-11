#!/usr/bin/env Rscript

# Covariate-scale (near-separation) diagnostic.
#
# Question: is the systematic eta over-estimation and the ~5x precision gap vs
# the paper driven by the strong-covariate / near-separation regime implied by a
# literal reading of Simulation 1 (X ~ N(0,1), beta = (1,2,3,-4,-5))?
#
# Xb has SD = sqrt(1+4+9+16+25) = 7.42 at xsd = 1, so ~50% of spots are near
# covariate-deterministic. We keep the generative + fitting model fixed (paper
# plug-in conditional logistic / MPLE) and only shrink the covariate SD `xsd`,
# which scales the linear-predictor magnitude by xsd. If eta bias and eta SD
# collapse toward the paper targets (~0.09, ~0.05) as xsd decreases, separation
# is the dominant cause and the paper used a milder LP scale than a literal read.
#
# GLM (MPLE) is used for speed and to isolate separation from any ADVI/prior
# artifact; earlier diagnostics already showed GLM reproduces the eta bias.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_file), "../.."))
source(file.path(root, "repro/paper_reproduce/paper_models.R"))

outdir <- file.path(root, "repro/paper_reproduce/output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_rep <- as.integer(Sys.getenv("N_REP", "50"))
eta_values <- as.numeric(strsplit(Sys.getenv("ETA_VALUES", "0.4,1.6,2.8"), ",")[[1]])
xsd_values <- as.numeric(strsplit(Sys.getenv("XSD_VALUES", "1.0,0.75,0.5,0.35,0.25"), ",")[[1]])

beta_true <- c(1, 2, 3, -4, -5)
lp_sd_unit <- sqrt(sum(beta_true^2)) # 7.416: SD of Xb per unit covariate SD

# Generate a Simulation-1-style dataset with covariate SD = xsd (beta truth fixed).
gen_dataset_xsd <- function(side = 30L, p = 20L, eta = 0.4, xsd = 1.0, seed = 1L) {
  set.seed(seed)
  coords <- grid_coords(side)
  n <- nrow(coords)
  X <- matrix(rnorm(n * p, sd = xsd), n, p)
  colnames(X) <- paste0("beta", seq_len(p))
  beta <- c(beta_true, rep(0, p - 5L))
  y <- gibbs_autologistic(X, coords, beta, eta, sweeps = 2000L)
  list(X = X, coords = coords, y = y, beta = beta, eta = eta)
}

# Plug-in conditional logistic (MPLE), neighbor-average convention.
fit_glm_eta <- function(d) {
  cvec <- neighbor_sum_count(d$coords, d$y)$avg
  dat <- data.frame(y = d$y, d$X, spatial = cvec)
  fit <- suppressWarnings(glm(y ~ . - 1, data = dat, family = binomial()))
  co <- coef(summary(fit))
  est <- coef(fit)
  eta_se <- if ("spatial" %in% rownames(co)) unname(co["spatial", "Std. Error"]) else NA_real_
  list(
    eta = unname(est["spatial"]),
    eta_se = eta_se,
    beta = unname(est[paste0("beta", 1:5)])
  )
}

rows <- list()
cat("== Covariate-scale (near-separation) diagnostic ==\n")
cat(sprintf("replicates=%d | eta=%s | xsd=%s\n\n",
            n_rep, paste(eta_values, collapse = ","), paste(xsd_values, collapse = ",")))

t0 <- Sys.time()
for (eta in eta_values) {
  for (xsd in xsd_values) {
    for (rep_id in seq_len(n_rep)) {
      seed <- 500000L + as.integer(round(eta * 1000)) * 1000L +
        as.integer(round(xsd * 1000)) + rep_id
      d <- gen_dataset_xsd(eta = eta, xsd = xsd, seed = seed)
      f <- fit_glm_eta(d)
      # separation flag: any |Xb| very large fraction, and glm fitted probs at 0/1
      xb <- as.numeric(d$X %*% d$beta)
      rows[[length(rows) + 1L]] <- data.frame(
        eta_true = eta,
        xsd = xsd,
        lp_sd = xsd * lp_sd_unit,
        rep = rep_id,
        disease_rate = mean(d$y),
        frac_extreme = mean(abs(xb) > 5),
        eta_hat = f$eta,
        eta_se = f$eta_se,
        eta_signed = f$eta - eta,
        eta_abs = abs(f$eta - eta),
        beta_infl = mean(f$beta / beta_true),  # mean signal-beta inflation ratio
        stringsAsFactors = FALSE
      )
    }
  }
  cat(sprintf("  eta=%.1f done (%.1f min elapsed)\n",
              eta, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

raw <- do.call(rbind, rows)

agg_fun <- function(d) {
  data.frame(
    eta_true = d$eta_true[1],
    xsd = d$xsd[1],
    lp_sd = round(d$lp_sd[1], 2),
    disease_rate = round(mean(d$disease_rate), 3),
    eta_mean = round(mean(d$eta_hat), 3),
    eta_signed = round(mean(d$eta_signed), 3),
    eta_absBias = round(mean(d$eta_abs), 3),
    eta_empSD = round(sd(d$eta_hat), 3),      # avgSEM analog (spread across reps)
    eta_meanSE = round(mean(d$eta_se), 3),    # avgSEE analog (GLM asymptotic SE)
    beta_infl = round(mean(d$beta_infl), 3),
    n = nrow(d),
    stringsAsFactors = FALSE
  )
}
groups <- split(raw, list(raw$eta_true, raw$xsd), drop = TRUE)
agg <- do.call(rbind, lapply(groups, agg_fun))
agg <- agg[order(agg$eta_true, -agg$xsd), ]
rownames(agg) <- NULL

raw_path <- file.path(outdir, sprintf("covariate_scale_raw_n%s.csv", n_rep))
agg_path <- file.path(outdir, sprintf("covariate_scale_aggregate_n%s.csv", n_rep))
write.csv(raw, raw_path, row.names = FALSE)
write.csv(agg, agg_path, row.names = FALSE)

cat("\nWrote", raw_path, "\n")
cat("Wrote", agg_path, "\n\n")
cat("Paper targets: eta absBias ~ 0.087/0.083/0.095 ; eta SD ~ 0.044-0.070\n\n")
print(agg, row.names = FALSE)
