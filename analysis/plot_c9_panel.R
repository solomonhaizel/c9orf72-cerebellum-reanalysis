#!/usr/bin/env Rscript
# Fig 1b: total C9ORF72, control vs C9-ALS (reanalysis of GSE67196).
# R port of plot_c9_panel.py. Reads data/total_c9_normalized.tsv, writes
# Fig1b_C9ORF72.{pdf,svg,png}. Mean +/- SEM, two-sided Mann-Whitney U.
suppressMessages(library(ggplot2))

args    <- commandArgs(trailingOnly = TRUE)
infile  <- if (length(args) >= 1) args[1] else "data/total_c9_normalized.tsv"
outbase <- if (length(args) >= 2) args[2] else "Fig1b_C9ORF72"

d <- read.delim(infile, stringsAsFactors = FALSE)
d <- d[d$group %in% c("control", "c9ALS"), ]
d$group <- factor(d$group, levels = c("control", "c9ALS"), labels = c("Control", "C9-ALS"))
d$x <- as.integer(d$group)

# p-value comes from the SciPy stats step (analysis/c9_stats.py); the figure
# does not recompute the test. Run c9_stats.py first to produce results/stats.tsv.
st   <- read.delim("results/stats.tsv", stringsAsFactors = FALSE)
pval <- as.numeric(st$p_value[1])
ptxt <- if (pval < 0.001) sprintf("P = %.4f", pval) else sprintf("P = %.3f", pval)

agg <- do.call(rbind, lapply(split(d, d$group), function(g) data.frame(
  group = g$group[1], x = g$x[1], mean = mean(g$norm_count),
  sem = sd(g$norm_count) / sqrt(nrow(g)))))

fill_cols <- c("Control" = "#8A8A8A", "C9-ALS" = "#1B6CA8")
edge_cols <- c("Control" = "#4A4A4A", "C9-ALS" = "#0E4D7A")

set.seed(7)
d$xj <- d$x + runif(nrow(d), -0.10, 0.10)
top <- max(d$norm_count); ybr <- top * 1.06; tick <- top * 0.035

g <- ggplot() +
  geom_segment(data = agg, aes(x = x - 0.23, xend = x + 0.23, y = mean, yend = mean,
               colour = group), linewidth = 0.9) +
  geom_errorbar(data = agg, aes(x = x, ymin = mean - sem, ymax = mean + sem,
               colour = group), width = 0.10, linewidth = 0.6) +
  geom_point(data = d, aes(xj, norm_count, fill = group, colour = group),
             shape = 21, size = 2.4, stroke = 0.5) +
  annotate("segment", x = 1, xend = 1, y = ybr, yend = ybr + tick, colour = "#222222", linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = ybr, yend = ybr + tick, colour = "#222222", linewidth = 0.4) +
  annotate("segment", x = 1, xend = 2, y = ybr + tick, yend = ybr + tick, colour = "#222222", linewidth = 0.4) +
  annotate("text", x = 1.5, y = ybr + top * 0.05, label = ptxt, fontface = "italic", size = 3.1, vjust = 0) +
  scale_fill_manual(values = fill_cols) +
  scale_colour_manual(values = edge_cols) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Control", "C9-ALS"), limits = c(0.4, 2.6)) +
  scale_y_continuous(limits = c(0, top * 1.28), expand = c(0, 0)) +
  labs(x = NULL, y = expression(atop("Total " * italic("C9ORF72"), "(normalized counts)")),
       title = "b.") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13, hjust = 0),
        axis.text = element_text(colour = "black", size = 9),
        axis.title.y = element_text(size = 9),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        axis.line = element_line(colour = "black", linewidth = 0.5))

ggsave(paste0(outbase, ".pdf"), g, width = 2.4, height = 2.8)
ggsave(paste0(outbase, ".png"), g, width = 2.4, height = 2.8, dpi = 600)
ggsave(paste0(outbase, ".svg"), g, width = 2.4, height = 2.8)

cat(sprintf("control n=%d mean=%.1f\n", sum(d$group=="Control"), mean(d$norm_count[d$group=="Control"])))
cat(sprintf("c9ALS   n=%d mean=%.1f\n", sum(d$group=="C9-ALS"), mean(d$norm_count[d$group=="C9-ALS"])))
cat(sprintf("fold=%.3f  p=%.4f\n", mean(d$norm_count[d$group=="C9-ALS"])/mean(d$norm_count[d$group=="Control"]), pval))
