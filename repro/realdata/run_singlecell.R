#!/usr/bin/env Rscript
# Single-cell / imaging-based spatial pipeline (Xenium / MERFISH / STARmap+).
# Demonstrates the method at single-cell SCALE and low DEPTH -- the regime that
# motivates both the O(N) inference and the NB count layer, and that the original
# O(N^2) log-Gaussian model cannot handle.
#
# Data: if a real dataset is provided (see SINGLECELL_DATA.md) point DATA_RDS to an
# .rds list(counts=cells x genes, coords=cells x 2 [x,y], disease=0/1). Otherwise a
# SYNTHETIC single-cell dataset is generated so the pipeline runs out-of-the-box.
root <- "/Users/licheng/Documents/DiSTect"
for (f in c("neighbors.R","likelihood_nb.R","fit_polyagamma.R","selection_fdr.R"))
  source(file.path(root,"repro/method",f))
set.seed(2024)
DATA_RDS <- Sys.getenv("DATA_RDS", "")

if (nzchar(DATA_RDS) && file.exists(DATA_RDS)) {
  cat("== Loading real single-cell dataset:", DATA_RDS, "==\n")
  d <- readRDS(DATA_RDS)
  # single-cell coords are continuous -> snap to a lattice at the neighborhood scale
  # so the rook-neighbor engine applies (radius ~ typical cell spacing).
  gr <- as.integer(as.factor(round(d$coords[,1]/20))); gc_ <- as.integer(as.factor(round(d$coords[,2]/20)))
  coords <- data.frame(x=gr, y=gc_)
  counts <- d$counts; y <- d$disease
} else {
  cat("== No DATA_RDS given -> SYNTHETIC single-cell dataset (Xenium-like) ==\n")
  side <- 250; N <- side*side                 # 62,500 cells
  gx<-((seq_len(N)-1)%%side)+1; gy<-((seq_len(N)-1)%/%side)+1
  coords<-data.frame(x=gx,y=gy)
  P<-120                                       # targeted panel
  z<-matrix(rnorm(N*P),N,P); colnames(z)<-paste0("gene",1:P)
  beta<-c(runif(6,1.0,1.8)*sample(c(-1,1),6,TRUE), rep(0,P-6))
  y<-rbinom(N,1,0.25)
  for(s in 1:15) y<-rbinom(N,1,plogis(-1.5 + scale(z)%*%beta + 1.0*neighbor_sum(coords,y)))
  depth<-0.8; size_i<-exp(rnorm(N,0,0.6))      # low depth + capture variation
  counts<-matrix(rnbinom(N*P, mu=pmax(depth*exp(0.5*z)*size_i,1e-3), size=1.5),N,P)
  colnames(counts)<-colnames(z)
}

cat(sprintf("  %d cells x %d genes | %.0f%% disease | mean count %.2f, %.0f%% zeros\n",
    nrow(counts), ncol(counts), 100*mean(y), mean(counts), 100*mean(counts==0)))

## NB covariates + PG-CAVI at single-cell scale
X <- scale(nb_pearson(counts)); X[is.na(X)]<-0
t <- system.time(fit <- fit_pgcavi_single(y, X, coords=coords))["elapsed"]
gt <- fit$table[!fit$table$param %in% c("intercept","eta"), ]  # gene rows (works for real symbols)
sel <- select_fdr_z(gt$std_effect, level=0.05)
cat(sprintf("\n== PG-CAVI at single-cell scale: %d cells fit in %.1fs, %d iters ==\n",
            nrow(counts), t, fit$iter))
cat(sprintf("  eta=%.3f | %d genes selected @FDR 0.05\n", fit$eta, sum(sel$selected)))
cat("  top disease-associated genes:",
    paste(head(gt$param[order(-gt$std_effect)],8), collapse=","), "\n")
cat("\nDONE. (The original O(N^2) log-Gaussian model is infeasible at this N.)\n")
