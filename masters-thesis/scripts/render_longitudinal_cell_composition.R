library(ggplot2)
library(dplyr)
library(tidyr)

metadata_file <- "/Users/debangana/Thesis/Parkinson/annotation/metadataUpdated.csv"
celltype_file <- "/Users/debangana/Thesis/Parkinson/annotation/cellcountMinfiEpic.csv"
output_file <- "masters-thesis/assets/figures/longitudinal-cell-composition.png"
summary_file <- "masters-thesis/assets/figures/longitudinal-cell-composition-slopes.csv"

celltype_order <- c("Neu", "NK", "Mono")

celltype_labels <- c(
  Neu = "Neutrophils",
  NK = "NK cells",
  Mono = "Monocytes"
)

stage_colours <- c(
  PreConverter = "#d95f8d",
  PostConverter = "#d95b43"
)

stage_labels <- c(
  PreConverter = "Pre-converter",
  PostConverter = "Post-converter"
)

metadata <- read.csv(metadata_file, check.names = FALSE)
celltype <- read.csv(celltype_file, check.names = FALSE)
celltype$Sample_Name <- metadata$Sample_Name

metadata <- metadata[!grepl("^7004_", metadata$Sample_Name), ]
celltype <- celltype[!grepl("^7004_", celltype$Sample_Name), ]

celltype_metadata <- merge(metadata, celltype, by = "Sample_Name", all = TRUE) %>%
  select(Sample_Name, TREND.ID, group, stage, visit.date, diag, age, Neu, NK, Mono) %>%
  filter(group == "converter", stage %in% c("PreConverter", "PostConverter")) %>%
  mutate(
    TREND.ID = factor(TREND.ID),
    stage = factor(stage, levels = c("PreConverter", "PostConverter")),
    visit.date = as.Date(visit.date),
    diag = as.Date(diag),
    years_relative_to_diag = as.numeric(difftime(visit.date, diag, units = "days")) / 365.25
  ) %>%
  filter(!is.na(years_relative_to_diag))

long_data <- celltype_metadata %>%
  pivot_longer(
    cols = all_of(celltype_order),
    names_to = "CellType",
    values_to = "Proportion"
  ) %>%
  mutate(CellType = factor(CellType, levels = celltype_order))

slope_results <- data.frame()

for (this_celltype in celltype_order) {
  model_data <- long_data %>% filter(CellType == this_celltype)
  fit <- lm(Proportion ~ years_relative_to_diag, data = model_data)
  slope <- unname(coef(fit)[["years_relative_to_diag"]])

  slope_results <- rbind(
    slope_results,
    data.frame(
      CellType = this_celltype,
      slope_per_year_relative_to_diagnosis = slope,
      slope_per_5_years_relative_to_diagnosis = slope * 5,
      stringsAsFactors = FALSE
    )
  )
}

write.csv(slope_results, summary_file, row.names = FALSE)

slope_labels <- slope_results %>%
  mutate(
    CellType = factor(CellType, levels = celltype_order),
    label = paste0("slope/year = ", sprintf("%+.4f", slope_per_year_relative_to_diagnosis))
  )

label_positions <- long_data %>%
  group_by(CellType) %>%
  summarise(
    x = min(years_relative_to_diag, na.rm = TRUE),
    y = max(Proportion, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x = x + 0.04 * diff(range(long_data$years_relative_to_diag, na.rm = TRUE)),
    y = y - 0.06 * diff(range(long_data$Proportion, na.rm = TRUE)),
    CellType = factor(CellType, levels = celltype_order)
  ) %>%
  left_join(slope_labels, by = "CellType")

plot <- ggplot(long_data, aes(x = years_relative_to_diag, y = Proportion)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#f2dcff", alpha = 0.45) +
  geom_vline(xintercept = 0, color = "#1f252c", linewidth = 0.75) +
  geom_line(
    aes(group = TREND.ID, color = TREND.ID),
    linewidth = 0.68,
    alpha = 0.58
  ) +
  geom_point(
    aes(fill = stage, color = TREND.ID),
    shape = 21,
    size = 2.7,
    stroke = 0.45,
    alpha = 0.98
  ) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = FALSE,
    color = "#9f0000",
    linewidth = 1.35,
    linetype = "longdash"
  ) +
  geom_text(
    data = label_positions,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.3,
    fontface = "bold",
    color = "#8a1f1f"
  ) +
  facet_wrap(~CellType, scales = "free_y", ncol = 3, labeller = labeller(CellType = celltype_labels)) +
  scale_fill_manual(values = stage_colours, labels = stage_labels) +
  scale_x_continuous(breaks = seq(-12, 6, 3)) +
  labs(
    title = "Cell Composition Relative to Diagnosis in Converters",
    subtitle = "Converters aligned to diagnosis date.",
    x = "Years relative to diagnosis",
    y = "Cell type proportion",
    fill = "Stage"
  ) +
  guides(color = "none") +
  theme_minimal(base_family = "Arial", base_size = 14) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0, color = "#1f252c"),
    plot.subtitle = element_text(size = 9.5, hjust = 0, color = "#54616f", margin = margin(b = 6)),
    axis.title = element_text(size = 14, face = "bold", color = "#1f252c"),
    axis.text = element_text(size = 11, color = "#3d4752"),
    strip.text = element_text(size = 15, face = "bold", color = "#1f252c"),
    panel.grid.major = element_line(color = "#e5e7eb", linewidth = 0.8),
    panel.grid.minor = element_line(color = "#f1f3f5", linewidth = 0.45),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(12, 18, 8, 12)
  )

ggsave(output_file, plot, width = 13.2, height = 6.3, dpi = 320, bg = "white")
