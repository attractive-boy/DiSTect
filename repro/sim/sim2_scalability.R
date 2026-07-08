#!/usr/bin/env Rscript
# Sim 2 -- SCALABILITY: PG-CAVI is O(N) and the fastest of the three engines.
# Adds PG-CAVI timings across N and merges with repro/m0/scaling.csv (dense-ADVI vs
# sparse-ADVI) into a 3-way comparison + figure.
root <- "/Users/licheng/Documents/DiSTect"
source(file.path(root, "repro/method/neighbors.R"))
source(file.path(root, "repro/method/fit_polyagamma.R"))
suppressMessages(library(ggplot2))
set.seed(7)

gen_fit_time <- function(N, P = 20) {
  side <- ceiling(sqrt(N))
  gx <- ((seq_len(N)-1) %% side)+1; gy <- ((seq_len(N)-1) %/% side)+1
  coords <- data.frame(x = gx, y = gy)
  X <- matrix(rnorm(N*P), N, P)
  beta <- c(1.5,-1.2, rep(0,P-2))
  y <- rbinom(N,1,plogis(X%*%beta + as.integer((gx%/%5+gy%/%5)%%2)*1.5 - 1))
  cc <- neighbor_sum(coords, y)
  system.time(fit_pgcavi_single(y, X, c_vec = cc))["elapsed"]
}

Ns <- c(500,1000,2000,3000,10000,30000,100000,300000)
cat("== Sim 2: PG-CAVI scaling ==\n")
pg <- data.frame(method="PG-CAVI", N=Ns,
                 fit_seconds = sapply(Ns, function(n){ t<-gen_fit_time(n); cat(sprintf("  N=%d: %.2fs\n",n,t)); t }))

scsv <- file.path(root, "repro/m0/scaling.csv")
out  <- pg
if (file.exists(scsv)) {
  adv <- read.csv(scsv)[, c("method","N","fit_seconds")]
  adv$method <- ifelse(adv$method=="dense","dense-ADVI O(N^2)","sparse-ADVI O(N)")
  out <- rbind(adv, pg)
}
out <- out[!is.na(out$fit_seconds), ]
write.csv(out, file.path(root,"repro/sim/sim2_scaling.csv"), row.names = FALSE)

p <- ggplot(out, aes(N, fit_seconds, color = method)) +
  geom_line(linewidth=1) + geom_point(size=2.3) +
  scale_x_log10() + scale_y_log10() + annotation_logticks(sides="bl") +
  labs(x="number of spots N (log)", y="fit time, s (log)",
       title="Inference scaling: dense-ADVI vs sparse-ADVI vs PG-CAVI", color=NULL) +
  theme_bw() + theme(legend.position = c(0.24,0.82))
ggsave(file.path(root,"repro/sim/fig_sim2_scaling.png"), p, width=6.5, height=4.5, dpi=120)
cat("\nwrote repro/sim/sim2_scaling.csv + fig_sim2_scaling.png\n")
print(out, row.names = FALSE)
