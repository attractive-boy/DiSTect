#!/usr/bin/env Rscript
# Sim 3 -- COUNT REGIME: why log(X+1) is the wrong covariate at single-cell / low
# sequencing depth, and how the NB size-factor normalization (method/likelihood_nb.R)
# fixes it. Two reproducible diagnostics:
#   (3a) VARIANCE STABILIZATION: log(X+1) has a strong mean-variance dependence
#        (heteroscedasticity that biases high-dim DE); NB Pearson residuals are flat.
#   (3b) DEPTH INVARIANCE: NB removes the per-cell depth confound; log does not.
# These are the mechanistic justifications; the predictive payoff is expected on real
# low-depth data (Xenium/MERFISH) and is validated there in the real-data pipeline.
root <- "/Users/licheng/Documents/DiSTect"
source(file.path(root, "repro/method/likelihood_nb.R"))
set.seed(23)

# genes across a wide TRUE expression range, observed at low depth with per-cell size
sim_counts <- function(N=2000, P=200, depth=1) {
  base <- exp(rnorm(P, -0.5, 1.2))                 # per-gene expression level (wide range)
  size_i <- exp(rnorm(N, 0, 0.6))                  # per-cell depth confound
  mu <- depth * outer(size_i, base)
  counts <- matrix(rnbinom(N*P, mu=pmax(mu,1e-4), size=1.5), N, P)
  list(counts=counts, base=base, size_i=size_i)
}

## (3a) mean-variance dependence: |cor(gene mean, gene variance)| across covariate
cat("== Sim 3a: variance stabilization (lower |cor(mean,var)| = better) ==\n")
tab <- data.frame()
for (depth in c(0.5, 1, 2, 5)) {
  d <- sim_counts(depth=depth)
  mv <- function(M){ m<-colMeans(M); v<-apply(M,2,var); abs(cor(m, v)) }
  logc <- mv(log1p(d$counts)); nbc <- mv(nb_pearson(d$counts))
  tab <- rbind(tab, data.frame(depth=depth, log_meanvar_cor=round(logc,3),
                               nb_meanvar_cor=round(nbc,3)))
}
print(tab, row.names = FALSE)

## (3b) depth invariance: correlation of covariate with the per-cell depth confound
##      (a good covariate should be ~uncorrelated with technical depth)
cat("\n== Sim 3b: residual depth confound (lower |cor(cell score, depth)| = better) ==\n")
tab2 <- data.frame()
for (depth in c(0.5, 1, 2, 5)) {
  d <- sim_counts(depth=depth)
  score <- function(M) rowMeans(M)                 # per-cell aggregate covariate
  logd <- abs(cor(score(log1p(d$counts)), d$size_i))
  nbd  <- abs(cor(score(nb_pearson(d$counts)),  d$size_i))
  tab2 <- rbind(tab2, data.frame(depth=depth, log_depth_cor=round(logd,3),
                                 nb_depth_cor=round(nbd,3)))
}
print(tab2, row.names = FALSE)
write.csv(cbind(tab, tab2[,-1]), file.path(root,"repro/sim/sim3_count_regime.csv"), row.names=FALSE)
cat("\nExpect: NB has far lower mean-variance coupling AND far lower depth confound than log(X+1),\n")
cat("        most pronounced at low depth -- the single-cell regime.\n")
