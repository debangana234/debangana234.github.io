library(grid)
library(png)
library(ggplot2)
library(dplyr)

figure_dir <- "masters-thesis/assets/figures"

stage_colors <- c(
  early_control = "#9ED9E8",
  late_control = "#12D7D7",
  PreConverter = "#F7B8C8",
  PostConverter = "#F26F63",
  control = "#12D7D7",
  Parkinson = "#F26F63"
)

site_theme <- theme_minimal(base_size = 17) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 15, color = "#111820"),
    axis.text = element_text(face = "bold", size = 12, color = "#3c3f43"),
    strip.text = element_text(face = "bold", size = 15, color = "#111820"),
    panel.grid.major = element_line(color = "#e7e7e7", linewidth = 0.45),
    panel.grid.minor = element_line(color = "#f1f1f1", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(14, 18, 14, 18),
    panel.spacing = unit(1.2, "lines")
  )

make_public_plot <- function() {
  public_file <- "/Users/debangana/Thesis/Parkinson/annotation/GSE111629AnnoEpiAge.csv"
  public_df <- read.csv(public_file, check.names = FALSE)
  names(public_df)[names(public_df) == ""] <- "row_id"
  public_df <- public_df |>
    filter(!is.na(`age.ch1`), !is.na(Horvath), !is.na(`disease.state.ch1`)) |>
    mutate(
      stage = factor(`disease.state.ch1`, levels = c("control", "Parkinson")),
      age = as.numeric(`age.ch1`)
    )

  r2_df <- public_df |>
    group_by(stage) |>
    summarise(r2 = summary(lm(Horvath ~ age, data = cur_data()))$r.squared, .groups = "drop") |>
    mutate(
      label = paste0("R² = ", sprintf("%.2f", r2)),
      x = c(40, 58),
      y = c(78, 78)
    )

  public_plot <- ggplot(public_df, aes(x = age, y = Horvath, color = stage)) +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "#7A140E", fill = "#C9C9C9", linewidth = 1.05) +
    geom_text(
      data = r2_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 5.2,
      fontface = "bold",
      color = "#111820"
    ) +
    facet_wrap(~ stage, ncol = 2) +
    scale_color_manual(values = stage_colors) +
    labs(
      title = "Horvath Age vs. Chronological age (public dataset)",
      x = "age",
      y = "Horvath"
    ) +
    coord_cartesian(ylim = c(28, 82), clip = "off") +
    site_theme

  ggsave(
    file.path(figure_dir, "horvath-clock-performance-horvath-dataset.png"),
    public_plot,
    width = 8.4,
    height = 6.0,
    units = "in",
    dpi = 320,
    bg = "white"
  )
}

make_trend_plot_from_original <- function() {
  trend_data_file <- "masters-thesis/assets/data/trend_horvath_age_estimates.csv"
  if (file.exists(trend_data_file)) {
    trend_df <- read.csv(trend_data_file, check.names = FALSE) |>
      filter(!is.na(age), !is.na(Horvath), !is.na(stage)) |>
      mutate(
        stage = factor(stage, levels = c("early_control", "late_control", "PreConverter", "PostConverter"))
      )

    r2_df <- trend_df |>
      group_by(stage) |>
      summarise(r2 = summary(lm(Horvath ~ age, data = cur_data()))$r.squared, .groups = "drop") |>
      mutate(
        label = paste0("R² = ", sprintf("%.2f", r2)),
        x = c(51, 55, 51, 68),
        y = c(80, 80, 80, 80)
      )

    trend_plot <- ggplot(trend_df, aes(x = age, y = Horvath, color = stage)) +
      geom_point(size = 2.5, alpha = 0.86) +
      geom_smooth(method = "lm", se = TRUE, color = "#7A140E", fill = "#C9C9C9", linewidth = 1.1) +
      geom_text(
        data = r2_df,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = 0,
        size = 5.4,
        fontface = "bold",
        color = "#111820"
      ) +
      facet_wrap(~ stage, ncol = 2) +
      scale_color_manual(values = stage_colors) +
      labs(
        title = "Horvath Age vs Chronological Age (TREND dataset)",
        x = "age",
        y = "Horvath"
      ) +
      coord_cartesian(ylim = c(42, 83), clip = "off") +
      site_theme +
      theme(
        plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 12)),
        axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 16)
      )

    ggsave(
      file.path(figure_dir, "horvath-clock-performance-trend.png"),
      trend_plot,
      width = 8.4,
      height = 6.0,
      units = "in",
      dpi = 320,
      bg = "white"
    )
    return(invisible(NULL))
  }

  source_file <- "/Users/debangana/Thesis/Parkinson/Rplots/EpigeneticAge/Horvath_vs_ChoronoAgeTrend.png"
  output_file <- file.path(figure_dir, "horvath-clock-performance-trend.png")
  img <- readPNG(source_file)

  # Remove the original right-side legend and use the plot panels as the data layer.
  cropped <- img[, 1:1135, , drop = FALSE]

  png(output_file, width = dim(cropped)[2], height = dim(cropped)[1], bg = "white")
  grid.raster(cropped, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

  grid.rect(x = 0.5, y = 0.982, width = 1, height = 0.05, gp = gpar(col = NA, fill = "white"))
  grid.rect(x = 0.47, y = 0.94, width = 0.78, height = 0.035, gp = gpar(col = NA, fill = "white"))
  grid.rect(x = 0.47, y = 0.492, width = 0.78, height = 0.055, gp = gpar(col = NA, fill = "white"))

  grid.rect(x = 0.20, y = 0.875, width = 0.24, height = 0.075, gp = gpar(col = NA, fill = "white"))
  grid.rect(x = 0.64, y = 0.875, width = 0.30, height = 0.075, gp = gpar(col = NA, fill = "white"))
  grid.rect(x = 0.20, y = 0.43, width = 0.24, height = 0.075, gp = gpar(col = NA, fill = "white"))
  grid.rect(x = 0.64, y = 0.43, width = 0.30, height = 0.075, gp = gpar(col = NA, fill = "white"))

  grid.text(
    "Horvath Age vs Chronological Age (TREND dataset)",
    x = unit(0.075, "npc"),
    y = unit(0.989, "npc"),
    just = c("left", "top"),
    gp = gpar(fontsize = 22, fontface = "bold", col = "#111820")
  )

  facet_labels <- data.frame(
    label = c("early_control", "late_control", "PreConverter", "PostConverter"),
    x = c(0.285, 0.71, 0.285, 0.71),
    y = c(0.94, 0.94, 0.495, 0.495)
  )
  for (i in seq_len(nrow(facet_labels))) {
    grid.text(
      facet_labels$label[i],
      x = unit(facet_labels$x[i], "npc"),
      y = unit(facet_labels$y[i], "npc"),
      just = c("center", "center"),
      gp = gpar(fontsize = 20, fontface = "bold", col = "#111820")
    )
  }

  r2_labels <- data.frame(
    label = c("R² = 0.54", "R² = 0.40", "R² = 0.48", "R² = 0.28"),
    x = c(0.205, 0.615, 0.205, 0.615),
    y = c(0.875, 0.875, 0.43, 0.43)
  )
  for (i in seq_len(nrow(r2_labels))) {
    grid.text(
      r2_labels$label[i],
      x = unit(r2_labels$x[i], "npc"),
      y = unit(r2_labels$y[i], "npc"),
      just = c("center", "center"),
      gp = gpar(fontsize = 34, fontface = "bold", col = "#111820")
    )
  }

  dev.off()
}

make_public_plot()
make_trend_plot_from_original()
