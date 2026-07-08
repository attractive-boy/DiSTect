#!/usr/bin/env Rscript
# Full discovery report for the real Xenium breast run (cluster-based tumor labels).
# Refits PG-CAVI once, characterizes convergence, and dumps BOTH tails of the gene table
# (tumor-enriched vs tumor-depleted) + in-sample predictive fit. Saves table + fit for reuse.
root <- "/Users/licheng/Documents/DiSTect"
for (f in c("neighbors.R","likelihood_nb.R","fit_polyagamma.R","selection_fdr.R"))
  source(file.path(root,"repro/method",f))
set.seed(2024)
RDS <- Sys.getenv("DATA_RDS","/tmp/xenium_breast/xenium_breast_rep1_clusterlab.rds")
d   <- readRDS(RDS); y <- d$disease
gr  <- as.integer(as.factor(round(d$coords[,1]/20))); gcy <- as.integer(as.factor(round(d$coords[,2]/20)))
coords <- data.frame(x=gr, y=gcy)
c_vec  <- neighbor_sum(coords, y)
X <- scale(nb_pearson(d$counts)); X[is.na(X)] <- 0
MI <- as.integer(Sys.getenv("MAX_ITER","800"))

t  <- system.time(fit <- fit_pgcavi_single(y, X, c_vec=c_vec, max_iter=MI))["elapsed"]
tr <- fit$elbo; rel <- abs(diff(tail(tr,2)))/abs(tail(tr,1))
cat(sprintf("N=%d P=%d | fit %.0fs, %d iters | eta=%.3f | ELBO rel-change=%.2e (%s)\n",
    nrow(X), ncol(X), t, fit$iter, fit$eta, rel, ifelse(fit$iter<MI,"CONVERGED","hit cap")))

## in-sample predictive fit
Z  <- cbind(1, X, c_vec); p <- plogis(as.numeric(Z %*% fit$raw$mean))
r  <- rank(p); n1 <- as.numeric(sum(y)); n0 <- as.numeric(sum(y==0))   # as.numeric: avoid int overflow
auc <- (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0)
cat(sprintf("in-sample: acc=%.3f  Brier=%.3f  AUC=%.3f\n", mean((p>0.5)==y), mean((p-y)^2), auc))

## gene table, both tails
gt <- fit$table[!fit$table$param %in% c("intercept","eta"), ]
gt$q <- p.adjust(2*pnorm(-gt$std_effect), "BH")
gt <- gt[order(-gt$std_effect), ]
cat(sprintf("genes @FDR0.05: %d/%d\n", sum(gt$q<0.05), nrow(gt)))
cat("\n== TUMOR-ENRICHED (mean>0), top 10 ==\n"); print(head(gt[gt$mean>0, c("param","mean","std_effect","q")],10), row.names=FALSE)
cat("\n== TUMOR-DEPLETED (mean<0), top 10 ==\n"); print(head(gt[gt$mean<0, c("param","mean","std_effect","q")],10), row.names=FALSE)
write.csv(gt, "/tmp/xenium_breast/xenium_discovery_table.csv", row.names=FALSE)
saveRDS(fit, "/tmp/xenium_breast/xenium_fit.rds")
cat("\nwrote xenium_discovery_table.csv + xenium_fit.rds\n")
