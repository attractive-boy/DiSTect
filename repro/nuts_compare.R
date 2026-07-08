# NUTS vs VI comparison on the single-slice toy data (DiSTect btaf530)
suppressMessages({ library(rstan); library(ggplot2); library(dplyr) })
options(mc.cores = parallel::detectCores())
set.seed(128)
seed <- 128   # dsgd_single's NUTS branch references a global `seed` (package bug workaround)

root <- "/Users/licheng/Documents/DiSTect"
source(file.path(root, "R", "dsgd.R"))
load(file.path(root, "data", "toy.rda"))
outdir <- file.path(root, "repro"); dir.create(outdir, showWarnings = FALSE)

x <- as.matrix(toy[, c("gene1","gene2","gene3","x","y")])
y <- toy$disease

cat("\n========== VI (fullrank ADVI) ==========\n")
t_vi <- system.time(fit_vi <- dsgd_single(y, x, method = "VI"))
s_vi <- rstan::summary(fit_vi)$summary
cat(sprintf("VI wall time: %.1f s\n", t_vi["elapsed"]))

cat("\n========== NUTS (HMC, iter=800, warmup=400, 2 chains) ==========\n")
t_nuts <- system.time(
  fit_nuts <- dsgd_single(y, x, method = "NUTS", niter = 800, nwarmup = 400, nchain = 2)
)
s_nuts <- rstan::summary(fit_nuts)$summary
cat(sprintf("NUTS wall time: %.1f s\n", t_nuts["elapsed"]))

pars <- c("beta[1]","beta[2]","beta[3]","eta","w")
cmp <- data.frame(
  param      = pars,
  VI_mean    = round(s_vi[pars,  "mean"], 3),
  NUTS_mean  = round(s_nuts[pars,"mean"], 3),
  VI_sd      = round(s_vi[pars,  "sd"],   3),
  NUTS_sd    = round(s_nuts[pars,"sd"],   3),
  NUTS_Rhat  = round(s_nuts[pars,"Rhat"], 3),
  NUTS_neff  = round(s_nuts[pars,"n_eff"], 0)
)
cat("\n---- Point estimate / uncertainty comparison ----\n")
print(cmp, row.names = FALSE)

cat(sprintf("\nSpeed-up: NUTS is %.1fx slower than VI (%.1fs vs %.1fs)\n",
            t_nuts["elapsed"]/t_vi["elapsed"], t_nuts["elapsed"], t_vi["elapsed"]))

saveRDS(list(cmp = cmp, t_vi = t_vi["elapsed"], t_nuts = t_nuts["elapsed"]),
        file.path(outdir, "nuts_vs_vi.rds"))
cat("\n========== NUTS-VI DONE ==========\n")
