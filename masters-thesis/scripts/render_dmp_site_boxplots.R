library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)

beta_file <- "/Users/debangana/Thesis/Parkinson/GCB_figures/PTPRN2_DMP_GCB2026/extracted_PTPRN2_and_significant_TREND_CpGs.csv"
annotation_file <- "/Users/debangana/Thesis/Parkinson/GCB_figures/PTPRN2_DMP_GCB2026/TREND_significant_DMP_CpG_heatmap_sample_annotation.csv"
output_dir <- "masters-thesis/assets/figures"

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

selected_sites <- tibble::tribble(
  ~CpG_ID, ~Gene, ~Feature, ~Position, ~Direction, ~Output,
  "cg19000186", "CNGA1", "5'UTR", "chr4:47956011", "Hypomethylated", "dmp-site-cg19000186.png",
  "cg11821245", "LDHC", "TSS200", "chr11:18433683", "Hypomethylated", "dmp-site-cg11821245.png",
  "cg09775972", "OR4K1", "1stExon", "chr14:20404281", "Hypermethylated", "dmp-site-cg09775972.png",
  "cg17316096", "PTPRN2", "Body", "chr7:158158218", "Hypermethylated", "dmp-site-cg17316096.png"
)

beta_values <- read.csv(beta_file, row.names = 1, check.names = FALSE)
sample_annotation <- read.csv(annotation_file, row.names = 1, check.names = FALSE)
sample_annotation$stage <- factor(sample_annotation$stage, levels = stage_order)
sample_annotation$Stage <- factor(stage_labels[as.character(sample_annotation$stage)], levels = unname(stage_labels))
sample_annotation$Sample_Name <- rownames(sample_annotation)

plot_data <- beta_values[selected_sites$CpG_ID, , drop = FALSE] |>
  rownames_to_column("CpG_ID") |>
  pivot_longer(cols = -CpG_ID, names_to = "Sample_Name", values_to = "Beta") |>
  left_join(sample_annotation[, c("Sample_Name", "Stage")], by = "Sample_Name") |>
  filter(!is.na(Stage)) |>
  left_join(selected_sites, by = "CpG_ID")

stage_counts <- sample_annotation |>
  filter(!is.na(Stage)) |>
  count(Stage, name = "n")

make_site_plot <- function(site_row) {
  site_data <- plot_data |>
    filter(CpG_ID == site_row$CpG_ID)
  y_range <- range(site_data$Beta, na.rm = TRUE)
  y_pad <- max(0.025, diff(y_range) * 0.12)
  y_top <- y_range[2] + y_pad
  y_bottom <- max(0, y_range[1] - y_pad)
  count_data <- stage_counts |>
    mutate(y = y_top)

  site_title <- paste0(
    site_row$CpG_ID,
    " | Gene: ",
    site_row$Gene,
    " (",
    site_row$Feature,
    "), ",
    site_row$Position
  )

  ggplot(site_data, aes(x = Stage, y = Beta, fill = Stage)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.7, color = "#2f343a", alpha = 0.86) +
    geom_jitter(width = 0.08, size = 1.45, alpha = 0.9, color = "#0B3CFF") +
    geom_text(
      data = count_data,
      aes(x = Stage, y = y, label = paste0("n: ", n)),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4.1
    ) +
    scale_fill_manual(values = stage_colors) +
    coord_cartesian(ylim = c(y_bottom, y_top + y_pad * 0.45)) +
    labs(
      title = site_title,
      subtitle = site_row$Direction,
      x = NULL,
      y = "Methylation beta value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#e7eaee", linewidth = 0.42),
      panel.grid.minor = element_line(color = "#f2f4f6", linewidth = 0.28),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 10.5, color = "#1f252c", hjust = 0.5),
      plot.subtitle = element_text(face = "bold", size = 10.5, color = ifelse(site_row$Direction == "Hypermethylated", "#dc2626", "#2563eb"), hjust = 0.5),
      axis.title.y = element_text(face = "bold", color = "#1f252c", size = 10.5),
      axis.text.x = element_text(face = "bold", color = "#343a40", size = 8.5, angle = 20, hjust = 1),
      axis.text.y = element_text(face = "bold", color = "#343a40", size = 9),
      plot.margin = margin(12, 12, 12, 12)
    )
}

for (row_index in seq_len(nrow(selected_sites))) {
  site_row <- selected_sites[row_index, ]
  ggsave(
    file.path(output_dir, site_row$Output),
    make_site_plot(site_row),
    width = 6.0,
    height = 4.2,
    dpi = 420,
    bg = "white"
  )
}
