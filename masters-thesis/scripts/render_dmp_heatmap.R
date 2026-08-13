library(pheatmap)
library(grid)

beta_file <- "/Users/debangana/Thesis/Parkinson/GCB_figures/PTPRN2_DMP_GCB2026/TREND_significant_DMP_CpG_heatmap_beta_values.csv"
annotation_file <- "/Users/debangana/Thesis/Parkinson/GCB_figures/PTPRN2_DMP_GCB2026/TREND_significant_DMP_CpG_heatmap_sample_annotation.csv"
output_file <- "masters-thesis/assets/figures/dmp-heatmap-latepost-tight.png"

stage_order <- c("early_control", "late_control", "PreConverter", "PostConverter")
stage_labels <- c(
  "early_control" = "Early control",
  "late_control" = "Late control",
  "PreConverter" = "Pre-converter",
  "PostConverter" = "Post-converter"
)

stage_colors <- c(
  "Early control" = "#BFDDE8",
  "Late control" = "#3DC9C9",
  "Pre-converter" = "#EAB6C3",
  "Post-converter" = "#EC7F75"
)

beta_values <- read.csv(beta_file, row.names = 1, check.names = FALSE)
sample_annotation <- read.csv(annotation_file, row.names = 1, check.names = FALSE)
sample_annotation$stage <- factor(sample_annotation$stage, levels = stage_order)
sample_annotation <- sample_annotation[order(sample_annotation$stage, sample_annotation$age), , drop = FALSE]
sample_annotation$Stage <- factor(stage_labels[as.character(sample_annotation$stage)], levels = unname(stage_labels))
sample_annotation$Age <- sample_annotation$age

annotation_display <- sample_annotation[, c("Age", "Stage"), drop = FALSE]
late_post_samples <- rownames(annotation_display)[annotation_display$Stage %in% c("Late control", "Post-converter")]
late_post_matrix <- as.matrix(beta_values[, late_post_samples[late_post_samples %in% colnames(beta_values)], drop = FALSE])
annotation_late_post <- droplevels(annotation_display[colnames(late_post_matrix), , drop = FALSE])

heatmap_colors <- colorRampPalette(c("#2563EB", "#F8FAFC", "#DC2626"))(101)
heatmap_breaks <- seq(0, 1, length.out = length(heatmap_colors) + 1)

age_colors <- colorRampPalette(c("#E5F3F3", "#87CFC5", "#1C7C70"))(100)
annotation_colors_late_post <- list(
  Age = age_colors,
  Stage = stage_colors[c("Late control", "Post-converter")]
)

bold_heatmap_text <- function(grob) {
  if (inherits(grob, "text")) {
    grob$gp$font <- NULL
    grob$gp$fontface <- "bold"
    grob$gp$col <- "#1f252c"
  }
  if (!is.null(grob$children)) {
    grob$children <- lapply(grob$children, bold_heatmap_text)
  }
  if (!is.null(grob$grobs)) {
    grob$grobs <- lapply(grob$grobs, bold_heatmap_text)
  }
  grob
}

heatmap_plot <- pheatmap(
  late_post_matrix,
  color = heatmap_colors,
  breaks = heatmap_breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_late_post,
  annotation_colors = annotation_colors_late_post,
  show_colnames = FALSE,
  show_rownames = FALSE,
  border_color = NA,
  fontsize = 12,
  fontsize_row = 5,
  fontsize_col = 5,
  annotation_names_col = FALSE,
  legend = TRUE,
  annotation_legend = TRUE,
  main = "91 DMPs between late control (cyan) and PostConverter (pink)",
  silent = TRUE
)

heatmap_plot$gtable <- bold_heatmap_text(heatmap_plot$gtable)

png(output_file, width = 10.4, height = 6.8, units = "in", res = 260, bg = "white")
grid.draw(heatmap_plot$gtable)
dev.off()
