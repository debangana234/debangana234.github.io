library(ggplot2)

input_file <- "/Users/debangana/Thesis/Parkinson/annotation/top_cpgsLateControlPostConv.csv"
output_file <- "masters-thesis/assets/figures/dmp-volcano-polished.png"

dmp_table <- read.csv(input_file, row.names = 1, check.names = FALSE)
dmp_table$CpG_ID <- rownames(dmp_table)
dmp_table$neg_log10_p <- -log10(pmax(dmp_table$P.Value, .Machine$double.xmin))
dmp_table$direction <- "Not significant"
dmp_table$direction[dmp_table$adj.P.Val < 0.05 & dmp_table$logFC > 0] <- "Hypermethylated"
dmp_table$direction[dmp_table$adj.P.Val < 0.05 & dmp_table$logFC < 0] <- "Hypomethylated"

label_ids <- c(
  head(dmp_table$CpG_ID[dmp_table$direction == "Hypermethylated"], 4),
  head(dmp_table$CpG_ID[dmp_table$direction == "Hypomethylated"], 4)
)
label_table <- dmp_table[dmp_table$CpG_ID %in% label_ids, ]

volcano_plot <- ggplot(dmp_table, aes(x = logFC, y = neg_log10_p)) +
  geom_point(
    data = subset(dmp_table, direction == "Not significant"),
    color = "#9aa5b1",
    alpha = 0.26,
    size = 0.58
  ) +
  geom_point(
    data = subset(dmp_table, direction == "Hypomethylated"),
    aes(color = direction),
    alpha = 0.9,
    size = 1.05
  ) +
  geom_point(
    data = subset(dmp_table, direction == "Hypermethylated"),
    aes(color = direction),
    alpha = 0.9,
    size = 1.05
  ) +
  geom_hline(yintercept = -log10(0.05), color = "#6b7280", linetype = "dashed", linewidth = 0.35) +
  geom_vline(xintercept = c(-0.5, 0.5), color = "#6b7280", linetype = "dashed", linewidth = 0.35) +
  geom_text(
    data = label_table,
    aes(label = CpG_ID, color = direction),
    size = 2.45,
    fontface = "bold",
    vjust = -0.65,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Hypomethylated" = "#2563eb",
      "Hypermethylated" = "#dc2626"
    ),
    breaks = c("Hypermethylated", "Hypomethylated")
  ) +
  coord_cartesian(xlim = c(-2.7, 2.7), ylim = c(0, 12.8), clip = "off") +
  labs(
    x = "Log2 fold change",
    y = "-log10(P-value)",
    color = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "#e7eaee", linewidth = 0.45),
    panel.grid.minor = element_line(color = "#f1f3f5", linewidth = 0.3),
    axis.title = element_text(face = "bold", color = "#1f252c"),
    axis.text = element_text(color = "#4b5563"),
    legend.position = "right",
    legend.text = element_text(face = "bold", color = "#1f252c", size = 8),
    legend.key.height = unit(0.18, "in"),
    legend.key.width = unit(0.18, "in"),
    plot.margin = margin(12, 18, 12, 12)
  )

ggsave(output_file, volcano_plot, width = 5.2, height = 3.8, dpi = 420, bg = "white")
