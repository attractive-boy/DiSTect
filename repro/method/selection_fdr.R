# FDR-controlled disease-gene selection (Direction B).
#
# DiSTect selects genes by a hard |mean/sd| > 1.96 threshold with NO multiplicity
# control. This module provides two principled selectors:
#  (1) frequentist-style: treat standardized effect as a z-score -> BH q-values;
#  (2) Bayesian: from spike-and-slab inclusion probabilities gamma_j, control the
#      expected Bayesian FDR (Newton et al. 2004).

# BH q-values from standardized effects (|mean/sd| ~ |z|)
select_fdr_z <- function(std_effect, level = 0.10) {
  p <- 2 * pnorm(-abs(std_effect))
  q <- p.adjust(p, method = "BH")
  data.frame(std_effect = std_effect, p = p, qvalue = q,
             selected = q < level)
}

# Bayesian FDR from inclusion probabilities: sort by gamma desc, include the
# largest set whose average local-fdr (1 - gamma) stays below `level`.
select_fdr_bayes <- function(gamma, level = 0.10) {
  ord  <- order(gamma, decreasing = TRUE)
  lfdr <- 1 - gamma[ord]
  cummean <- cumsum(lfdr) / seq_along(lfdr)
  k    <- sum(cummean <= level)
  sel  <- logical(length(gamma)); sel[ord[seq_len(k)]] <- TRUE
  data.frame(gamma = gamma, local_fdr = 1 - gamma, selected = sel)
}
