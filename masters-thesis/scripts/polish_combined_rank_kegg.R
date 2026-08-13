library(ggplot2)

source_file <- "/Users/debangana/Thesis/Parkinson/annotation/top1000gProfiler.csv"
output_file <- "masters-thesis/assets/figures/combined-rank-kegg.png"

kegg_df <- read.csv(source_file, check.names = FALSE)
kegg_df <- subset(kegg_df, source == "KEGG")
kegg_df <- kegg_df[order(kegg_df$adjusted_p_value), ]
kegg_df$pathway <- gsub(" - ", " -\n", kegg_df$term_name)
kegg_df$pathway <- gsub("Metabolism of xenobiotics by\ncytochrome", "Metabolism of xenobiotics by\ncytochrome", kegg_df$pathway)
kegg_df$pathway <- factor(kegg_df$pathway, levels = rev(kegg_df$pathway))

plot <- ggplot(kegg_df, aes(x = intersection_size, y = pathway, fill = adjusted_p_value)) +
  geom_col(width = 0.84) +
  scale_fill_gradient(
    low = "#d9414c",
    high = "#2c7fb8",
    name = "p.adjust",
    breaks = c(0.01, 0.02, 0.03, 0.04)
  ) +
  scale_x_continuous(
    breaks = c(0, 2.5, 5, 7.5, 10, 12.5),
    limits = c(0, 12.1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "Top KEGG Pathway Enrichment",
    x = "Count",
    y = NULL
  ) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 28, hjust = 0.5, margin = margin(b = 18)),
    axis.text.y = element_text(face = "bold", size = 17, color = "#111111", lineheight = 0.95),
    axis.text.x = element_text(face = "bold", size = 17, color = "#111111"),
    axis.title.x = element_text(face = "bold", size = 20, color = "#111111", margin = margin(t = 10)),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 14),
    panel.grid.major = element_line(color = "#dcdcdc", linewidth = 0.55),
    panel.grid.minor = element_line(color = "#ededed", linewidth = 0.35),
    plot.margin = margin(18, 24, 18, 22),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(output_file, plot, width = 10.8, height = 8.0, dpi = 220, bg = "white")
