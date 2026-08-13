library(ggplot2)

beta_file <- "/Users/debangana/Thesis/Parkinson/GCB_figures/DMR_GCB2026/extracted_DMR_CpGs_beta_values.csv"
region_file <- "/Users/debangana/Thesis/Parkinson/GCB_figures/DMR_GCB2026/DMR_region_table_GCB2026.csv"
annotation_file <- "/Users/debangana/Thesis/Parkinson/GCB_figures/DMR_GCB2026/DMR_heatmap_sample_annotation_GCB2026.csv"
output_dir <- "masters-thesis/assets/figures"
heatmap_colors <- grDevices::colorRampPalette(c("#2563EB", "#F8FAFC", "#DC2626"))(101)

targets <- data.frame(
  DMR_ID = c("DMR_3", "DMR_1"),
  gene = c("LDHC", "ACAP1"),
  direction = c("Hypomethylated region", "Hypermethylated region"),
  output = c("dmr-region-ldhc-thesis.png", "dmr-region-acap1-thesis.png"),
  stringsAsFactors = FALSE
)

beta_values <- read.csv(beta_file, row.names = 1, check.names = FALSE)
region_table <- read.csv(region_file, check.names = FALSE)
sample_annotation <- read.csv(annotation_file, row.names = 1, check.names = FALSE)

sample_annotation <- sample_annotation[
  sample_annotation$stage %in% c("late_control", "PostConverter"),
  ,
  drop = FALSE
]
sample_annotation$stage <- factor(sample_annotation$stage, levels = c("late_control", "PostConverter"))
sample_annotation <- sample_annotation[order(sample_annotation$stage, sample_annotation$age), , drop = FALSE]

plot_dmr_heatmap <- function(region_matrix, target, region) {
  row_means <- rowMeans(region_matrix, na.rm = TRUE)
  cpg_order <- names(sort(row_means, decreasing = TRUE))
  region_matrix <- region_matrix[cpg_order, rownames(sample_annotation), drop = FALSE]

  plot_data <- as.data.frame(as.table(region_matrix), stringsAsFactors = FALSE)
  colnames(plot_data) <- c("CpG", "Sample", "Beta")
  plot_data$Sample <- factor(plot_data$Sample, levels = rownames(sample_annotation))
  plot_data$CpG <- factor(plot_data$CpG, levels = rev(cpg_order))

  late_count <- sum(sample_annotation$stage == "late_control")
  stage_labels <- data.frame(
    x = c(late_count / 2 + 0.5, late_count + (nrow(sample_annotation) - late_count) / 2 + 0.5),
    y = length(cpg_order) + 1.08,
    label = c("LateControl", "PostConverter")
  )

  plot_data$x <- as.numeric(plot_data$Sample)
  plot_data$y <- as.numeric(plot_data$CpG)

  ggplot(plot_data, aes(x = x, y = y, fill = Beta)) +
    geom_tile(width = 0.98, height = 0.98) +
    annotate(
      "rect",
      xmin = 0.5,
      xmax = late_count + 0.5,
      ymin = length(cpg_order) + 0.62,
      ymax = length(cpg_order) + 0.86,
      fill = "#25cbd3"
    ) +
    annotate(
      "rect",
      xmin = late_count + 0.5,
      xmax = nrow(sample_annotation) + 0.5,
      ymin = length(cpg_order) + 0.62,
      ymax = length(cpg_order) + 0.86,
      fill = "#ef8c9a"
    ) +
    geom_text(
      data = stage_labels,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4.2,
      color = "#1f252c"
    ) +
    scale_fill_gradientn(
      colors = heatmap_colors,
      values = seq(0, 1, length.out = length(heatmap_colors)),
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = c("0", "0.25", "0.50", "0.75", "1"),
      name = "Beta"
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, nrow(sample_annotation) + 0.5)) +
    scale_y_continuous(
      expand = c(0.01, 0),
      breaks = seq_along(cpg_order),
      labels = rev(cpg_order),
      limits = c(0.5, length(cpg_order) + 1.25),
      position = "right"
    ) +
    labs(
      title = paste0(target$DMR_ID, ": ", target$gene, " | ", target$direction),
      subtitle = paste0(region$seqnames, ":", region$start, "-", region$end),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "#1f252c"),
      plot.subtitle = element_text(face = "bold", size = 11, hjust = 0.5, color = "#4b5563"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(face = "bold", size = 10, color = "#1f252c"),
      panel.grid = element_blank(),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(face = "bold", size = 8),
      plot.margin = margin(12, 12, 12, 12),
      plot.background = element_rect(fill = "#ffffff", color = NA),
      panel.background = element_rect(fill = "#ffffff", color = NA)
    )
}

for (target_index in seq_len(nrow(targets))) {
  target <- targets[target_index, ]
  region <- region_table[region_table$DMR_ID == target$DMR_ID, ]
  cpg_ids <- trimws(unlist(strsplit(region$cpg.probes, ",")))
  cpg_ids <- cpg_ids[cpg_ids %in% rownames(beta_values)]
  region_matrix <- as.matrix(beta_values[cpg_ids, rownames(sample_annotation), drop = FALSE])

  plot <- plot_dmr_heatmap(region_matrix, target, region)
  ggsave(
    filename = file.path(output_dir, target$output),
    plot = plot,
    width = 8.8,
    height = max(3.7, 2.6 + 0.34 * length(cpg_ids)),
    dpi = 220,
    bg = "white"
  )
}
