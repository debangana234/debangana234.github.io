library(ggplot2)
library(grid)

base_dir <- "/Users/debangana/Thesis/Parkinson"
figure_dir <- "masters-thesis/assets/figures"

stage_colors <- c(
  early_control = "#69B7C7",
  late_control = "#2B5CFF",
  PreConverter = "#F2A7B8",
  PostConverter = "#E84A3A"
)

stage_labels <- c(
  early_control = "Early control",
  late_control = "Late control",
  PreConverter = "Pre-converter",
  PostConverter = "Post-converter"
)

theme_web_plot <- function(base_size = 18) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      text = element_text(face = "bold", color = "black"),
      plot.title = element_blank(),
      axis.title = element_text(face = "bold", size = rel(1.0), color = "black"),
      axis.text = element_text(face = "bold", color = "grey20"),
      panel.grid.major = element_line(color = "grey86", linewidth = 0.65),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(face = "bold", color = "black"),
      plot.margin = margin(8, 18, 8, 12),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

## PCA of top 500 combined-rank CpGs
pca_scores <- read.csv(
  file.path(base_dir, "GCB_figures", "top500_combined_rank_PCA", "PCA_500_combined_rank_sites_TREND_scores.csv"),
  check.names = FALSE
)
pca_scores$stage <- factor(pca_scores$stage, levels = names(stage_colors))

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = stage)) +
  geom_hline(yintercept = 0, linewidth = 0.75, color = "grey82") +
  geom_vline(xintercept = 0, linewidth = 0.75, color = "grey82") +
  stat_ellipse(aes(fill = stage), geom = "polygon", type = "norm", level = 0.68,
               alpha = 0.11, linewidth = 0, show.legend = FALSE) +
  geom_point(size = 3.9, alpha = 0.92) +
  scale_color_manual(values = stage_colors, labels = stage_labels, drop = FALSE) +
  scale_fill_manual(values = stage_colors, labels = stage_labels, drop = FALSE) +
  labs(
    x = "PC1 (16.13% variance explained)",
    y = "PC2 (4.72% variance explained)",
    color = "Stage"
  ) +
  theme_web_plot(18) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 17)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4.6, alpha = 1)))

ggsave(file.path(figure_dir, "combined-rank-pca.png"), pca_plot,
       width = 8.8, height = 6.2, dpi = 260, bg = "white")

## GO enrichment bubble plot
go_data <- read.csv(
  file.path(base_dir, "GCB_figures", "poster_bold_plots", "poster_bold_GO_bubbleplot_data.csv"),
  check.names = FALSE
)
go_data$term_label <- factor(go_data$term_label, levels = rev(go_data$term_label))

go_plot <- ggplot(go_data, aes(x = negative_log10_of_adjusted_p_value, y = term_label)) +
  geom_point(aes(size = intersection_size, color = negative_log10_of_adjusted_p_value), alpha = 0.95) +
  scale_color_gradient(low = "#2c7fb8", high = "#d0183f", name = expression(-log[10]("adjusted p"))) +
  scale_size_continuous(name = "Gene count", range = c(4.2, 15), breaks = c(10, 20, 30, 40, 50)) +
  labs(x = expression(-log[10]("adjusted p-value")), y = NULL) +
  theme_web_plot(18) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 15, lineheight = 0.88),
    axis.text.x = element_text(size = 15),
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

ggsave(file.path(figure_dir, "combined-rank-go-enrichment.png"), go_plot,
       width = 8.8, height = 6.2, dpi = 260, bg = "white")

## KEGG pathway bar plot
kegg_source <- file.path(base_dir, "annotation", "top1000gProfiler.csv")
kegg_df <- read.csv(kegg_source, check.names = FALSE)
kegg_df <- subset(kegg_df, source == "KEGG")
kegg_df <- kegg_df[order(kegg_df$adjusted_p_value), ]
kegg_df$pathway <- gsub(" - ", " -\n", kegg_df$term_name)
kegg_df$pathway <- factor(kegg_df$pathway, levels = rev(kegg_df$pathway))

kegg_plot <- ggplot(kegg_df, aes(x = intersection_size, y = pathway, fill = adjusted_p_value)) +
  geom_col(width = 0.82) +
  scale_fill_gradient(low = "#d9414c", high = "#2c7fb8",
                      name = "p.adjust", breaks = c(0.01, 0.02, 0.03, 0.04)) +
  scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10, 12.5), limits = c(0, 12.1),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Count", y = NULL) +
  theme_web_plot(18) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 15, lineheight = 0.9),
    axis.text.x = element_text(size = 15),
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

ggsave(file.path(figure_dir, "combined-rank-kegg.png"), kegg_plot,
       width = 8.8, height = 6.2, dpi = 260, bg = "white")

## Selected model confusion matrix
confusion_counts <- read.csv(
  file.path(base_dir, "GCB_figures", "highly_variable_CpG_modeling", "conference_outputs", "selected_200CpG_celltype_confusion_matrix_counts.csv"),
  check.names = FALSE
)
confusion_counts$Prediction <- factor(confusion_counts$Prediction, levels = c("Control", "Parkinson"))
confusion_counts$Reference <- factor(confusion_counts$Reference, levels = c("Control", "Parkinson"))

confusion_plot <- ggplot(confusion_counts, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white", linewidth = 1.1) +
  geom_text(aes(label = Freq), size = 6.2, fontface = "bold", color = "#151515") +
  scale_fill_gradient(low = "#d8eef3", high = "#e76f66", name = "Count") +
  labs(x = "Actual label", y = "Predicted label") +
  coord_equal() +
  theme_minimal(base_size = 15, base_family = "Helvetica") +
  theme(
    axis.title = element_text(face = "bold", size = 15, color = "#151515"),
    axis.text = element_text(face = "bold", size = 13, color = "#333333"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12, color = "#151515"),
    legend.text = element_text(face = "bold", size = 11, color = "#333333"),
    panel.grid = element_blank(),
    plot.margin = margin(12, 14, 12, 12),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(figure_dir, "predictive-modeling-confusion-matrix.png"), confusion_plot,
       width = 6.6, height = 5.8, dpi = 260, bg = "white")
