#!/usr/bin/env Rscript

# Last confirmatory shot: meanfield ADVI (not fullrank) at reduced covariate
# scale xsd=0.5. The fullrank ADVI is pathological here (Pareto k up to 65) and
# gives eta SEE ~0.18, ~3.5x above the CRLB (0.05) and the paper (0.044). Test
# whether the rstan-default meanfield ADVI, at the scale where CRLB matches the
# paper, achieves SEE ~0.05 -- which would pin the reproduction gap to
# (mild scale) + (meanfield, not fullrank).

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- normalizePath(dirname(script_file))
source(file.path(root, "paper_models.R"))
outdir <- file.path(root, "output")

n_rep    <- as.integer(Sys.getenv("N_REP", "50"))
eta_vals <- as.numeric(strsplit(Sys.getenv("ETA_VALUES", "0.4"), ",")[[1]])
xsd_vals <- as.numeric(strsplit(Sys.getenv("XSD_VALUES", "0.5,1.0"), ",")[[1]])
algo     <- Sys.getenv("ALGO", "meanfield")

fit_single_mf <- function(y, X, coords) {
  c_avg <- neighbor_sum_count(coords, y)$avg
  data <- list(N = nrow(X), P = ncol(X), x = X, y = as.integer(y), c_avg = c_avg)
  model <- get_stan_model("paper_single", paper_single_stan)
  vb(model, data = data, algorithm = algo, seed = 128L, iter = 10000,
     tol_rel_obj = 1e-5, refresh = 0)
}
fit_naive_mf <- function(y, X) {
  data <- list(N = nrow(X), P = ncol(X), x = X, y = as.integer(y))
  model <- get_stan_model("paper_naive", paper_naive_stan)
  vb(model, data = data, algorithm = algo, seed = 128L, iter = 10000,
     tol_rel_obj = 1e-5, refresh = 0)
}
sim1_xsd <- function(side = 30L, p = 20L, eta = 0.4, xsd = 1.0, seed = 1L) {
  set.seed(seed)
  coords <- grid_coords(side); n <- nrow(coords)
  X <- matrix(rnorm(n * p, sd = xsd), n, p); colnames(X) <- paste0("beta", seq_len(p))
  beta <- c(1, 2, 3, -4, -5, rep(0, p - 5L))
  list(X = X, coords = coords, y = gibbs_autologistic(X, coords, beta, eta, sweeps = 2000L),
       beta = beta, eta = eta)
}

rows <- list()
cat(sprintf("== %s ADVI | reps=%d eta=%s xsd=%s ==\n\n", algo, n_rep,
            paste(eta_vals, collapse=","), paste(xsd_vals, collapse=",")))
t0 <- Sys.time()
for (xsd in xsd_vals) for (eta in eta_vals) {
  for (rep_id in seq_len(n_rep)) {
    seed <- 810000L + as.integer(round(xsd*1000))*1000L + as.integer(round(eta*1000)) + rep_id
    d <- sim1_xsd(eta = eta, xsd = xsd, seed = seed)
    fp <- tryCatch(fit_single_mf(d$y, d$X, d$coords), error=function(e) NULL)
    if (!is.null(fp)) { t <- fit_summary_table(fp, d$beta, truth_eta=eta); t$method<-"Proposed"; t$xsd<-xsd; rows[[length(rows)+1L]]<-t }
    fn <- tryCatch(fit_naive_mf(d$y, d$X), error=function(e) NULL)
    if (!is.null(fn)) { t <- fit_summary_table(fn, d$beta, truth_eta=NULL); t$method<-"Naive"; t$xsd<-xsd; rows[[length(rows)+1L]]<-t }
  }
  cat(sprintf("  xsd=%.2f eta=%.1f done (%.1f min)\n", xsd, eta, as.numeric(difftime(Sys.time(),t0,units="mins"))))
}
raw <- do.call(rbind, rows)
write.csv(raw, file.path(outdir, sprintf("confirm_meanfield_raw_n%s.csv", n_rep)), row.names=FALSE)
agg_fun <- function(d) data.frame(xsd=d$xsd[1], method=d$method[1], param=d$param[1],
  avgBias=round(mean(d$abs_bias),3), avgSEE=round(sd(d$mean),3), avgSEM=round(mean(d$sd),3), avgCI=round(mean(d$covered),3), n=nrow(d))
g <- split(raw, list(raw$xsd, raw$method, raw$param), drop=TRUE)
agg <- do.call(rbind, lapply(g, agg_fun)); rownames(agg)<-NULL
write.csv(agg, file.path(outdir, sprintf("confirm_meanfield_aggregate_n%s.csv", n_rep)), row.names=FALSE)
cat("\nPaper Table C1 (eta=0.4) Proposed ADVI: b1 bias .332 SEE .055; eta bias .087 SEE .044 CI .95\n")
cat("CRLB at xsd=0.5: SE(eta)~0.050 ; fullrank ADVI gave eta SEE 0.178\n\n")
show <- agg[agg$param %in% c("beta1","beta5","eta"), c("xsd","method","param","avgBias","avgSEE","avgSEM","avgCI")]
print(show[order(show$xsd, show$method, show$param),], row.names=FALSE)
