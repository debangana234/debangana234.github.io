library(dplyr)
library(ggplot2)
library(lmerTest)
library(emmeans)

metadata_file <- "/Users/debangana/Thesis/Parkinson/annotation/metadataUpdated.csv"
celltype_file <- "/Users/debangana/Thesis/Parkinson/annotation/cellcountMinfiEpic.csv"
output_file <- "masters-thesis/assets/figures/cell-composition-model-tests.txt"
plot_file <- "masters-thesis/assets/figures/cell-composition-mixed-model-results.png"

metadata <- read.csv(metadata_file, check.names = FALSE)
celltype <- read.csv(celltype_file, check.names = FALSE)
celltype$Sample_Name <- metadata$Sample_Name

metadata <- metadata[!grepl("^7004_", metadata$Sample_Name), ]
celltype <- celltype[!grepl("^7004_", celltype$Sample_Name), ]

model_data <- merge(metadata, celltype, by = "Sample_Name", all = TRUE) %>%
  filter(stage %in% c("early_control", "late_control", "PreConverter", "PostConverter")) %>%
  mutate(
    stage = factor(stage, levels = c("early_control", "late_control", "PreConverter", "PostConverter")),
    sex = factor(sex),
    TREND.ID = factor(TREND.ID),
    Lymphocyte = CD8T + CD4T + NK + Bcell,
    NL_ratio = Neu / Lymphocyte
  ) %>%
  select(Sample_Name, TREND.ID, stage, group, sex, age, Neu, Lymphocyte, NL_ratio) %>%
  filter(complete.cases(.))

response_labels <- c(
  Neu = "Neutrophils",
  Lymphocyte = "Total lymphocytes",
  NL_ratio = "Neutrophil-to-lymphocyte ratio"
)

stage_labels <- c(
  early_control = "Early\ncontrol",
  late_control = "Late\ncontrol",
  PreConverter = "Pre-\nconverter",
  PostConverter = "Post-\nconverter"
)

stage_colours <- c(
  early_control = "#87bcca",
  late_control = "#2f6f8f",
  PreConverter = "#d95f8d",
  PostConverter = "#d95b43"
)

estimated_means <- data.frame()
stage_tests <- data.frame()

fit_and_capture <- function(response) {
  formula_text <- paste0(response, " ~ stage + age + sex + (1 | TREND.ID)")
  model <- lmer(as.formula(formula_text), data = model_data, REML = FALSE)
  reduced <- lmer(as.formula(paste0(response, " ~ age + sex + (1 | TREND.ID)")), data = model_data, REML = FALSE)
  stage_lrt <- anova(reduced, model)
  stage_p <- stage_lrt$`Pr(>Chisq)`[2]
  stage_anova <- anova(model, type = 3)
  this_label <- response_labels[[response]]
  emmeans_df <- as.data.frame(emmeans(model, ~ stage)) %>%
    mutate(
      response = response,
      response_label = this_label,
      stage_p = stage_p
    )

  estimated_means <<- rbind(estimated_means, emmeans_df)
  stage_tests <<- rbind(
    stage_tests,
    data.frame(
      response = response,
      response_label = this_label,
      likelihood_ratio_p = stage_p,
      type_III_p = stage_anova["stage", "Pr(>F)"],
      stringsAsFactors = FALSE
    )
  )

  capture.output({
    cat("\n========================================\n")
    cat(response, "\n")
    cat("Model:", formula_text, "\n")
    cat("----------------------------------------\n")
    cat("Likelihood-ratio test for categorical stage effect:\n")
    print(stage_lrt)
    cat("\nFixed effects:\n")
    print(coef(summary(model)))
    cat("\nType III tests:\n")
    print(stage_anova)
    cat("\nPairwise stage contrasts, Tukey adjusted:\n")
    print(pairs(emmeans(model, ~ stage), adjust = "tukey"))
  })
}

sample_summary <- capture.output({
  cat("Sample counts by stage:\n")
  print(table(model_data$stage))
  cat("\nNumber of unique TREND participants:", length(unique(model_data$TREND.ID)), "\n")
  cat("\nParticipants with repeated measurements:", sum(table(model_data$TREND.ID) > 1), "\n")
})

results <- c(
  sample_summary,
  fit_and_capture("Neu"),
  fit_and_capture("Lymphocyte"),
  fit_and_capture("NL_ratio")
)

estimated_means$response_label <- factor(estimated_means$response_label, levels = response_labels)
estimated_means$stage <- factor(estimated_means$stage, levels = levels(model_data$stage))

stage_label_data <- stage_tests %>%
  mutate(
    response_label = factor(response_label, levels = response_labels),
    label = paste0("stage effect p=", signif(likelihood_ratio_p, 2))
  )

label_positions <- estimated_means %>%
  group_by(response_label) %>%
  summarise(
    x = 1,
    y = max(upper.CL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(stage_label_data, by = "response_label")

plot <- ggplot(estimated_means, aes(x = stage, y = emmean, group = response_label)) +
  geom_line(color = "#1f252c", linewidth = 0.75, alpha = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = stage), width = 0.11, linewidth = 0.85) +
  geom_point(aes(fill = stage), shape = 21, size = 3.8, stroke = 0.65, color = "white") +
  geom_text(
    data = label_positions,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -0.55,
    size = 3.4,
    fontface = "bold",
    color = "#1f252c"
  ) +
  facet_wrap(~response_label, scales = "free_y", ncol = 3) +
  scale_x_discrete(labels = stage_labels, drop = FALSE) +
  scale_color_manual(values = stage_colours, guide = "none") +
  scale_fill_manual(values = stage_colours, labels = gsub("\n", " ", stage_labels)) +
  labs(
    title = "Mixed-model Stage Effects",
    subtitle = "Estimated marginal means adjusted for age, sex, and repeated TREND measurements.",
    x = NULL,
    y = "Model-estimated value",
    fill = "Stage"
  ) +
  theme_minimal(base_family = "Arial", base_size = 14) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0, color = "#1f252c"),
    plot.subtitle = element_text(size = 10.5, hjust = 0, color = "#54616f", margin = margin(b = 8)),
    axis.title.y = element_text(size = 14, face = "bold", color = "#1f252c"),
    axis.text = element_text(size = 11, color = "#3d4752"),
    axis.text.x = element_text(face = "bold"),
    strip.text = element_text(size = 14, face = "bold", color = "#1f252c"),
    panel.grid.major = element_line(color = "#e5e7eb", linewidth = 0.8),
    panel.grid.minor = element_line(color = "#f1f3f5", linewidth = 0.45),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(12, 18, 8, 12)
  )

ggsave(plot_file, plot, width = 12.8, height = 5.8, dpi = 320, bg = "white")

writeLines(results, output_file)
cat(paste(results, collapse = "\n"))
