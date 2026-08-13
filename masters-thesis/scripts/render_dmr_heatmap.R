library(pheatmap)
library(grid)

beta_file <- "/Users/debangana/Thesis/Parkinson/GCB_figures/DMR_GCB2026/DMR_average_beta_values_GCB2026.csv"
annotation_file <- "/Users/debangana/Thesis/Parkinson/GCB_figures/DMR_GCB2026/DMR_heatmap_sample_annotation_GCB2026.csv"
output_file <- "masters-thesis/assets/figures/dmr-heatmap-beta-polished.png"

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

dmr_values <- read.csv(beta_file, row.names = 1, check.names = FALSE)
sample_annotation <- read.csv(annotation_file, row.names = 1, check.names = FALSE)
sample_annotation$stage <- factor(sample_annotation$stage, levels = stage_order)
sample_annotation <- sample_annotation[
  sample_annotation$stage %in% c("late_control", "PostConverter"),
  ,
  drop = FALSE
]
sample_annotation <- sample_annotation[order(sample_annotation$stage, sample_annotation$age), , drop = FALSE]
sample_annotation$Stage <- factor(stage_labels[as.character(sample_annotation$stage)], levels = unname(stage_labels))
sample_annotation$Age <- sample_annotation$age

annotation_display <- sample_annotation[, c("Stage", "Age"), drop = FALSE]
ordered_samples <- rownames(annotation_display)
dmr_matrix <- as.matrix(dmr_values[, ordered_samples[ordered_samples %in% colnames(dmr_values)], drop = FALSE])
annotation_display <- annotation_display[colnames(dmr_matrix), , drop = FALSE]

short_region_labels <- sub("^(DMR_[0-9]+ \\| [^|]+).*", "\\1", rownames(dmr_matrix))
rownames(dmr_matrix) <- short_region_labels

heatmap_colors <- colorRampPalette(c("#2563EB", "#F8FAFC", "#DC2626"))(101)
heatmap_breaks <- seq(0, 1, length.out = length(heatmap_colors) + 1)
age_colors <- colorRampPalette(c("#E5F3F3", "#87CFC5", "#1C7C70"))(100)

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
  dmr_matrix,
  color = heatmap_colors,
  breaks = heatmap_breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_display,
  annotation_colors = list(
    Stage = stage_colors[c("Late control", "Post-converter")],
    Age = age_colors
  ),
  show_colnames = FALSE,
  show_rownames = TRUE,
  border_color = NA,
  fontsize = 13,
  fontsize_row = 13,
  annotation_names_col = TRUE,
  legend = TRUE,
  annotation_legend = TRUE,
  main = "DMR-level methylation: late control vs post-converter",
  silent = TRUE
)

heatmap_plot$gtable <- bold_heatmap_text(heatmap_plot$gtable)

png(output_file, width = 10.8, height = 8.4, units = "in", res = 260, bg = "white")
grid.draw(heatmap_plot$gtable)
dev.off()
