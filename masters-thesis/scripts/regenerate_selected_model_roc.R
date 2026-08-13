library(ggplot2)
library(pROC)

prediction_file <- "/Users/debangana/Thesis/Parkinson/GCB_figures/highly_variable_CpG_modeling/conference_outputs/selected_200CpG_celltype_predictions.csv"
output_file <- "masters-thesis/assets/figures/predictive-modeling-roc-curve.png"

predictions <- read.csv(prediction_file, check.names = FALSE)
predictions$actual <- factor(predictions$actual, levels = c("Control", "Parkinson"))

roc_curve <- roc(
  response = predictions$actual,
  predictor = predictions$probability_Parkinson,
  levels = c("Control", "Parkinson"),
  direction = "<",
  quiet = TRUE
)

roc_df <- data.frame(
  false_positive_rate = 1 - roc_curve$specificities,
  true_positive_rate = roc_curve$sensitivities
)

auc_value <- as.numeric(auc(roc_curve))

roc_plot <- ggplot(roc_df, aes(x = false_positive_rate, y = true_positive_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#8b8f94", linewidth = 0.75) +
  geom_path(color = "#d9414c", linewidth = 1.2, lineend = "round") +
  annotate("label",
           x = 0.62,
           y = 0.18,
           label = paste0("AUC = ", round(auc_value, 3)),
           label.size = 0,
           fill = "white",
           color = "#111111",
           fontface = "bold",
           size = 4.6) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = TRUE) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_minimal(base_size = 15, base_family = "Helvetica") +
  theme(
    axis.title = element_text(face = "bold", size = 15, color = "#151515"),
    axis.text = element_text(face = "bold", size = 13, color = "#333333"),
    panel.grid.major = element_line(color = "#e0e0e0", linewidth = 0.55),
    panel.grid.minor = element_line(color = "#eeeeee", linewidth = 0.35),
    plot.margin = margin(12, 14, 12, 12),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(output_file, roc_plot, width = 6.6, height = 5.8, dpi = 260, bg = "white")
