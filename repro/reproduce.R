# Reproduction of DiSTect (Zhao et al., Bioinformatics 2025, btaf530)
# End-to-end demonstration of the paper's method on the packaged toy datasets.
suppressMessages({
  library(rstan)
  library(rjags)
  library(coda)
  library(sp)
  library(ggplot2)
  library(dplyr)
  library(GGally)
  library(network)
})
options(mc.cores = parallel::detectCores())
set.seed(128)

root <- "/Users/licheng/Documents/DiSTect"
# Source package functions (predict() here masks stats::predict on purpose)
for (f in c("dsgd.R", "prediction.R", "missing_imputation.R", "coef_plot.R", "network_plot.R"))
  source(file.path(root, "R", f))
load(file.path(root, "data", "toy.rda"))
load(file.path(root, "data", "toy2.rda"))
outdir <- file.path(root, "repro"); dir.create(outdir, showWarnings = FALSE)

sep <- function(t) cat("\n\n========== ", t, " ==========\n")

## ---- 1. Single-slice disease-associated gene detection (Sim 1 / main model) ----
sep("1. SINGLE-SLICE MODEL (toy: 300 spots, 3 genes)")
x  <- as.matrix(toy[, c("gene1","gene2","gene3","x","y")])
y  <- toy$disease
model.single <- dsgd(list_y = y, matrix_x = x)          # VI (fullrank)
s1 <- rstan::summary(model.single)$summary
print(round(s1[c("beta[1]","beta[2]","beta[3]","eta","w"),
              c("mean","sd","2.5%","97.5%")], 4))
cat("\nStandardized effect size |mean/sd| (gene ranking metric from paper):\n")
b <- s1[grep("^beta\\[", rownames(s1)), ]
z <- setNames(abs(b[,"mean"]/b[,"sd"]), c("gene1","gene2","gene3"))
print(round(sort(z, decreasing = TRUE), 3))
try(ggsave(file.path(outdir, "fig_coef_single.png"),
       plot_coef(model.single, x, n = 3), width = 5, height = 4, dpi = 120))

## ---- 2. Gene-gene interaction network (Downstream analysis ii) ----
sep("2. INTERACTION MODEL + NETWORK")
model.int <- dsgd(list_y = toy$disease,
                  matrix_x = toy[, c("gene1","gene2","gene3","x","y")],
                  interaction = c("gene1","gene2","gene3"))
si <- rstan::summary(model.int)$summary
bi <- si[grep("^beta\\[", rownames(si)), ]
lbl <- c("gene1","gene2","gene3","gene1*gene2","gene1*gene3","gene2*gene3")
zint <- data.frame(term = lbl, mean = round(bi[,"mean"],4), sd = round(bi[,"sd"],4),
                   std_effect = round(abs(bi[,"mean"]/bi[,"sd"]),3))
print(zint, row.names = FALSE)
cat("\nInteraction terms retained (|std effect| > 1.96, paper threshold):\n")
print(zint[4:6, ][abs(zint[4:6,"std_effect"]) > 1.96, ])
try(ggsave(file.path(outdir, "fig_network.png"),
       plot_network(model.int, c("gene1","gene2","gene3")), width = 5, height = 5, dpi = 120))

## ---- 3. Multiple-slice hierarchical model (Sim 2 / model 3) ----
sep("3. MULTIPLE-SLICE MODEL (toy2: 200 spots, label/random effect)")
x2 <- as.matrix(toy2[, c("gene1","gene2","gene3","x","y")])
y2 <- toy2$disease
label <- toy2$label
model.multi <- dsgd(list_y = y2, matrix_x = x2, label_list = label)
s3 <- rstan::summary(model.multi)$summary
keep <- intersect(c("beta[1]","beta[2]","beta[3]","eta","rho","sigma_square","w"),
                  rownames(s3))
print(round(s3[keep, c("mean","sd","2.5%","97.5%")], 4))

## ---- 4. Spatial disease prediction (Downstream analysis iii) ----
sep("4. PREDICTION on held-out spots (Gibbs, 100 sweeps)")
set.seed(1)
idx   <- sample(nrow(toy), 60)                 # held-out test spots
test  <- toy[idx, c("gene1","gene2","gene3","x","y")]
truth <- toy$disease[idx]
pred  <- predict(model.single, test, sweep = 100)
acc   <- mean(pred == truth)
cat(sprintf("Test spots: %d | Prediction accuracy: %.3f\n", length(truth), acc))
print(table(truth = truth, predicted = pred))

## ---- 5. Missing disease-status imputation (Downstream analysis iv) ----
sep("5. MISSING IMPUTATION (JAGS/Gibbs, neighbour-based)")
set.seed(7)
y_miss <- toy$disease
miss   <- sample(nrow(toy), 30)                # knock out 30 labels
true_v <- y_miss[miss]
y_miss[miss] <- NA
coords <- as.matrix(toy[, c("x","y")])
imp    <- missing_imputation(y_miss, coords)
imp_p  <- imp[miss]                            # posterior P(disease)
imp_hard <- as.integer(imp_p > 0.5)
cat(sprintf("Missing spots: %d | Imputation accuracy (>0.5): %.3f\n",
            length(true_v), mean(imp_hard == true_v)))
print(data.frame(spot = miss, true = true_v, post_prob = round(imp_p,3),
                 imputed = imp_hard)[1:10, ], row.names = FALSE)

sep("DONE")
cat("Figures written to", outdir, "\n")
