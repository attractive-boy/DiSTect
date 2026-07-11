#!/usr/bin/env Rscript

# Reverse-search candidate definitions for the published avgSEE/avgSEM columns
# using existing 200-replicate outputs. No model fitting is performed.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "repro/paper_reproduce/reverse_search_uncertainty_formulas.R"
root <- normalizePath(file.path(dirname(script_file), "../.."))
outdir <- file.path(root, "repro/paper_reproduce/output")

required <- file.path(outdir, c(
  "sim1_raw_n200.csv", "sim1_advi_key_comparison_n200.csv",
  "sim2_raw_n200.csv", "sim2_advi_key_comparison_n200.csv",
  "sim3_raw_n200.csv", "sim3_advi_key_comparison_n200.csv"
))
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing required existing outputs: ", paste(missing, collapse = ", "))

summarize_group <- function(d) {
  n <- nrow(d)
  err <- d$mean - d$truth
  c(
    empirical_sd = sd(d$mean),
    empirical_se_mean = sd(d$mean) / sqrt(n),
    rmse = sqrt(mean(err^2)),
    rmse_se_mean = sqrt(mean(err^2)) / sqrt(n),
    mean_abs_error = mean(abs(err)),
    mean_abs_error_se_mean = mean(abs(err)) / sqrt(n),
    mean_posterior_sd = mean(d$sd),
    mean_posterior_sd_se_mean = mean(d$sd) / sqrt(n),
    rms_posterior_sd = sqrt(mean(d$sd^2)),
    rms_posterior_sd_se_mean = sqrt(mean(d$sd^2)) / sqrt(n),
    mean_posterior_variance = mean(d$sd^2),
    sqrt_mean_posterior_variance_over_n = sqrt(mean(d$sd^2) / n)
  )
}

long_candidates <- function(raw, keys) {
  split_key <- interaction(raw[keys], drop = TRUE, lex.order = TRUE)
  parts <- split(raw, split_key)
  do.call(rbind, lapply(parts, function(d) {
    id <- d[1, keys, drop = FALSE]
    rownames(id) <- NULL
    vals <- summarize_group(d)
    rbind(
      cbind(id[rep(1, 6), , drop = FALSE], metric = "SEE",
            candidate = unname(names(vals)[1:6]), value = unname(vals[1:6]), row.names = NULL),
      cbind(id[rep(1, 6), , drop = FALSE], metric = "SEM",
            candidate = unname(names(vals)[7:12]), value = unname(vals[7:12]), row.names = NULL)
    )
  }))
}

# Standardize targets and raw group identifiers.
t1 <- read.csv(file.path(outdir, "sim1_advi_key_comparison_n200.csv"))
t1$table <- "C1-C3"
t1$group <- paste(t1$method, t1$eta_setting, t1$param, sep = "|")
t1 <- rbind(
  data.frame(table = t1$table, group = t1$group, method = t1$method,
             param = t1$param, metric = "SEE", target = t1$target_avgSEE),
  data.frame(table = t1$table, group = t1$group, method = t1$method,
             param = t1$param, metric = "SEM", target = t1$target_avgSEM)
)
r1 <- read.csv(file.path(outdir, "sim1_raw_n200.csv"))
c1 <- long_candidates(r1, c("method", "eta_setting", "param"))
c1$table <- "C1-C3"
c1$group <- paste(c1$method, c1$eta_setting, c1$param, sep = "|")

t2 <- read.csv(file.path(outdir, "sim2_advi_key_comparison_n200.csv"))
t2$table <- "C4"
t2$group <- paste(t2$method, t2$sigma2, t2$rho, t2$param, sep = "|")
t2 <- rbind(
  data.frame(table = t2$table, group = t2$group, method = t2$method,
             param = t2$param, metric = "SEE", target = t2$target_avgSEE),
  data.frame(table = t2$table, group = t2$group, method = t2$method,
             param = t2$param, metric = "SEM", target = t2$target_avgSEM)
)
r2 <- read.csv(file.path(outdir, "sim2_raw_n200.csv"))
c2 <- long_candidates(r2, c("method", "sigma2", "rho", "param"))
c2$table <- "C4"
c2$group <- paste(c2$method, c2$sigma2, c2$rho, c2$param, sep = "|")

t3 <- read.csv(file.path(outdir, "sim3_advi_key_comparison_n200.csv"))
t3$table <- "C5-C6"
t3$group <- paste(t3$method, t3$setting, t3$param, sep = "|")
t3 <- rbind(
  data.frame(table = t3$table, group = t3$group, method = t3$method,
             param = t3$param, metric = "SEE", target = t3$target_avgSEE),
  data.frame(table = t3$table, group = t3$group, method = t3$method,
             param = t3$param, metric = "SEM", target = t3$target_avgSEM)
)
r3 <- read.csv(file.path(outdir, "sim3_raw_n200.csv"))
c3 <- long_candidates(r3, c("method", "setting", "param"))
c3$table <- "C5-C6"
c3$group <- paste(c3$method, c3$setting, c3$param, sep = "|")

targets <- rbind(t1, t2, t3)
candidates <- rbind(
  c1[, c("table", "group", "metric", "candidate", "value")],
  c2[, c("table", "group", "metric", "candidate", "value")],
  c3[, c("table", "group", "metric", "candidate", "value")]
)
audit <- merge(targets, candidates, by = c("table", "group", "metric"))
audit$ratio <- audit$value / audit$target
audit$abs_log_error <- abs(log(audit$ratio))
audit$parameter_class <- ifelse(audit$param == "eta", "eta", "beta")

score_one <- function(d, scope) {
  data.frame(
    scope = scope,
    metric = d$metric[1],
    candidate = d$candidate[1],
    n_cells = nrow(d),
    median_ratio = median(d$ratio),
    median_abs_log_error = median(d$abs_log_error),
    rmse_log_error = sqrt(mean(log(d$ratio)^2)),
    within_25pct = mean(d$ratio >= 0.8 & d$ratio <= 1.25),
    within_2x = mean(d$ratio >= 0.5 & d$ratio <= 2),
    stringsAsFactors = FALSE
  )
}

score_scope <- function(d, scope) {
  groups <- split(d, list(d$metric, d$candidate), drop = TRUE)
  do.call(rbind, lapply(groups, score_one, scope = scope))
}

scores <- rbind(
  score_scope(audit, "all"),
  do.call(rbind, lapply(split(audit, audit$table), function(d) score_scope(d, d$table[1]))),
  do.call(rbind, lapply(split(audit, audit$parameter_class),
                       function(d) score_scope(d, paste0("parameter:", d$parameter_class[1]))))
)
scores <- scores[order(scores$scope, scores$metric, scores$rmse_log_error), ]

best <- do.call(rbind, lapply(split(audit, interaction(audit$table, audit$group, audit$metric,
                                                       drop = TRUE)), function(d) {
  d[which.min(d$abs_log_error), ]
}))
best <- best[order(best$table, best$group, best$metric), ]

num_score <- vapply(scores, is.numeric, logical(1))
scores[num_score] <- lapply(scores[num_score], function(x) round(x, 5))
num_best <- vapply(best, is.numeric, logical(1))
best[num_best] <- lapply(best[num_best], function(x) round(x, 5))

score_path <- file.path(outdir, "uncertainty_formula_search_scores.csv")
best_path <- file.path(outdir, "uncertainty_formula_search_best.csv")
write.csv(scores, score_path, row.names = FALSE)
write.csv(best, best_path, row.names = FALSE)

cat("Wrote", score_path, "\n")
cat("Wrote", best_path, "\n\n")
cat("Top candidates by scope and metric:\n")
top <- do.call(rbind, lapply(split(scores, list(scores$scope, scores$metric), drop = TRUE),
                            function(d) head(d, 3)))
print(top[, c("scope", "metric", "candidate", "n_cells", "median_ratio",
              "rmse_log_error", "within_2x")], row.names = FALSE)
