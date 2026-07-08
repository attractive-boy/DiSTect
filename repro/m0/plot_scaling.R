# Plot dense-vs-sparse scaling curves from scaling.csv
suppressMessages(library(ggplot2))
root <- "/Users/licheng/Documents/DiSTect"; d <- file.path(root, "repro/m0")
df <- read.csv(file.path(d, "scaling.csv"))
df <- df[!is.na(df$fit_seconds), ]
df$method <- factor(df$method, levels = c("dense","sparse"),
                    labels = c("dense O(N^2)", "sparse O(N)"))
cols <- c("dense O(N^2)" = "#C0392B", "sparse O(N)" = "#1F77B4")

mk <- function(y, ylab, file) {
  p <- ggplot(df, aes(N, .data[[y]], color = method)) +
    geom_line(linewidth = 1) + geom_point(size = 2.4) +
    scale_x_log10() + scale_y_log10() +
    scale_color_manual(values = cols) +
    annotation_logticks(sides = "bl") +
    labs(x = "number of spots N (log)", y = ylab, color = NULL,
         title = paste0("DiSTect fit ", ylab, ": dense vs sparse")) +
    theme_bw() + theme(legend.position = c(0.18, 0.85))
  ggsave(file.path(d, file), p, width = 6, height = 4.5, dpi = 120)
}
mk("fit_seconds", "fit time (s)",   "fig_scaling_time.png")
mk("peak_rss_mb", "peak RSS (MB)",  "fig_scaling_mem.png")

cat("== scaling.csv ==\n"); print(read.csv(file.path(d,"scaling.csv")), row.names = FALSE)
