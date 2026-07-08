#!/usr/bin/env Rscript
# Sim 4 -- PSEUDOLIKELIHOOD BIAS: fitting with observed neighbor labels attenuates
# the spatial autocorrelation eta; method/eta_fulllik.R corrects it by Monte-Carlo
# moment matching on the full autologistic. Shows eta_PL << eta_true, eta_debiased ~ eta_true.
root <- "/Users/licheng/Documents/DiSTect"
source(file.path(root, "repro/method/neighbors.R"))
source(file.path(root, "repro/method/fit_polyagamma.R"))
source(file.path(root, "repro/method/eta_fulllik.R"))
set.seed(3)

gen <- function(eta_true, side=40, P=8, b0=-2) {
  N<-side*side; gx<-((seq_len(N)-1)%%side)+1; gy<-((seq_len(N)-1)%/%side)+1
  coords<-data.frame(x=gx,y=gy); X<-matrix(rnorm(N*P),N,P); colnames(X)<-paste0("g",1:P)
  beta<-c(1.2,-1.0,rep(0,P-2)); y<-rbinom(N,1,0.3)
  for(s in 1:40) y<-rbinom(N,1,plogis(b0+X%*%beta+eta_true*neighbor_sum(coords,y)))
  list(y=y,X=X,coords=coords,beta=beta)
}

cat("== Sim 4: pseudolikelihood eta attenuation + full-likelihood debiasing ==\n")
out <- data.frame()
for (eta_true in c(0.5, 1.0, 1.5, 2.0)) {
  d <- gen(eta_true)
  fit <- fit_pgcavi_single(d$y, d$X, coords=d$coords)   # pseudolik eta
  db  <- debias_eta(d$y, d$X, d$coords, beta = fit$beta, offset = fit$intercept) # corrected
  out <- rbind(out, data.frame(eta_true=eta_true,
                               eta_pseudolik=round(fit$eta,3),
                               eta_debiased=round(db$eta,3)))
  cat(sprintf("  eta_true=%.1f | pseudolik=%.3f | debiased=%.3f\n",
              eta_true, fit$eta, db$eta))
}
write.csv(out, file.path(root,"repro/sim/sim4_pseudolik_bias.csv"), row.names=FALSE)
cat("\nExpect: eta_pseudolik systematically < eta_true; eta_debiased closer to the diagonal.\n")
