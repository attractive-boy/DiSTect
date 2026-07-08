#!/usr/bin/env Rscript
# Real-data pipeline: the FULL upgraded method on HER2+ breast (her2st, Andersson 2021).
# PG-CAVI inference (O(N)) + NB count covariates + FDR-controlled selection +
# full-likelihood eta debiasing + calibrated leave-one-patient-out prediction.
# Runs ALL 8 patients in seconds -- the regime the original 25h ADVI could not reach.
root <- "/Users/licheng/Documents/DiSTect"; her2 <- "/tmp/her2st/data"
for (f in c("neighbors.R","likelihood_nb.R","fit_polyagamma.R","selection_fdr.R",
            "eta_fulllik.R","predict_calibrated.R"))
  source(file.path(root, "repro/method", f))
set.seed(128)
if (!dir.exists(her2)) stop("Run repro/realdata/fetch_her2st.sh first.")

SECTIONS <- c(A="A1",B="B1",C="C1",D="D1",E="E1",F="F1",G="G2",H="H1")
N_HVG <- 50

load_section <- function(sec) {
  cnt <- read.delim(gzfile(file.path(her2,"ST-cnts",paste0(sec,".tsv.gz"))),
                    row.names=1, check.names=FALSE)
  lab <- read.delim(file.path(her2,"ST-pat","lbl",paste0(sec,"_labeled_coordinates.tsv")),
                    check.names=FALSE)
  lab$key <- paste0(round(lab$x),"x",round(lab$y)); lab<-lab[!duplicated(lab$key),]
  rownames(lab)<-lab$key
  common<-intersect(rownames(cnt),lab$key); cnt<-cnt[common,]; lab<-lab[common,]
  keep<-!grepl("undetermined",lab$label,ignore.case=TRUE); cnt<-cnt[keep,]; lab<-lab[keep,]
  co<-do.call(rbind,strsplit(rownames(cnt),"x"))
  list(counts=as.matrix(cnt),
       coords=data.frame(x=as.numeric(co[,1]),y=as.numeric(co[,2])),
       disease=as.integer(grepl("cancer",lab$label,ignore.case=TRUE)))
}

cat("== Load 8 patients ==\n")
secs <- lapply(SECTIONS, load_section)
for (i in seq_along(secs)) cat(sprintf("  %s: %d spots, %.0f%% cancer\n",
  names(SECTIONS)[i], nrow(secs[[i]]$counts), 100*mean(secs[[i]]$disease)))

# common gene filter (total count >= 300 in every section), NB covariates, top HVG
gene_ok <- Reduce(intersect, lapply(secs, function(s) colnames(s$counts)[colSums(s$counts)>=300]))
Xlist <- lapply(secs, function(s) nb_pearson(s$counts[, gene_ok, drop=FALSE]))
Xall  <- do.call(rbind, Xlist)
y     <- unlist(lapply(secs,`[[`,"disease"))
label <- rep(seq_along(secs), sapply(secs, function(s) nrow(s$counts)))
coord <- do.call(rbind, lapply(secs,`[[`,"coords"))
hvg   <- names(sort(apply(Xall,2,var),decreasing=TRUE))[1:N_HVG]
X     <- scale(Xall[,hvg]); X[is.na(X)] <- 0
cat(sprintf("\n%d spots, %d common genes -> top %d HVG (NB), %d patients, %.0f%% cancer\n",
            length(y), length(gene_ok), N_HVG, length(secs), 100*mean(y)))

## --- 1) disease-gene discovery (all 8 patients jointly) ---
cat("\n== Disease-gene discovery (PG-CAVI, all 8 patients) ==\n")
t <- system.time(fit <- fit_pgcavi_single(y, X, coords=coord, label=label))["elapsed"]
gt  <- fit$table[fit$table$param %in% hvg, ]
sel <- select_fdr_z(gt$std_effect, level=0.10)
gt$qvalue <- sel$qvalue; gt <- gt[order(-gt$std_effect), ]
cat(sprintf("  fit %.2fs (vs original ADVI ~25h) | eta=%.3f | %d genes selected @FDR 0.10\n",
            t, fit$eta, sum(sel$selected)))
cat("  top 15 disease-associated genes:\n")
print(head(transform(gt, mean=round(mean,2), std_effect=round(std_effect,1),
                     qvalue=signif(qvalue,2))[,c("param","mean","std_effect","qvalue")],15),
      row.names=FALSE)

## --- 2) full-likelihood eta debiasing ---
db <- debias_eta(y, X, coord, beta=fit$beta, offset=fit$intercept, label=label)
cat(sprintf("\n== eta: pseudolik=%.3f  full-likelihood debiased=%.3f ==\n", fit$eta, db$eta))

## --- 3) calibrated leave-one-patient-out prediction (spatial vs non-spatial) ---
cat("\n== Leave-one-patient-out prediction (calibrated) ==\n")
res <- data.frame()
for (p in seq_along(secs)) {
  tr<-which(label!=p); te<-which(label==p)
  f <- fit_pgcavi_single(y[tr], X[tr,], coords=coord[tr,], label=label[tr])
  A <- build_adjacency(coord[te,])
  th<-f$raw$mean; S<-f$raw$Sigma
  # spatial posterior-predictive; theta layout = [intercept, genes, eta]
  Xte <- cbind(1, X[te,])                          # include intercept column
  mth <- th[c(1, 2:(N_HVG+1))]; Sth <- S[c(1,2:(N_HVG+1)), c(1,2:(N_HVG+1))]
  eta_hat <- th[length(th)]
  pp <- predict_posterior(c(mth, eta_hat),
                          Matrix::bdiag(Sth, S[length(th),length(th)]) |> as.matrix(),
                          Xte, A, S=60, sweep=30)
  m <- calibration(pp, y[te])
  res <- rbind(res, cbind(patient=names(SECTIONS)[p], round(m,3)))
}
print(res, row.names=FALSE)
cat(sprintf("\nmean accuracy=%.3f  Brier=%.3f  ECE=%.3f\n",
            mean(as.numeric(res$accuracy)), mean(as.numeric(res$brier)),
            mean(as.numeric(res$ece))))
cat("\nDONE: full upgraded pipeline on all 8 HER2+ patients.\n")
