library(caret)
library(e1071)
library(pROC)

set.seed(42)

base_dir <- "/Users/debangana/Thesis/Parkinson"
model_data_file <- file.path(base_dir, "annotation", "model_data.csv")
ranked_cpg_file <- file.path(
  base_dir,
  "GCB_figures",
  "highly_variable_CpG_modeling",
  "model_metrics",
  "all_CpGs_ranked_by_within_dataset_variance.csv"
)
output_file <- "masters-thesis/assets/data/predictive-model-family-auroc.csv"

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

model_data <- read.csv(model_data_file, check.names = FALSE)
model_data <- model_data[, !colnames(model_data) %in% c("", "X")]
model_data$disease_state[model_data$disease_state == "control"] <- "Control"
model_data$disease_state <- factor(model_data$disease_state, levels = c("Control", "Parkinson"))

ranked_cpgs <- read.csv(ranked_cpg_file, check.names = FALSE)$CpG
site_counts <- c(100, 200, 300, 400)

split_group <- interaction(model_data$dataset, model_data$disease_state)
train_index <- createDataPartition(split_group, p = 0.8, list = FALSE)
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

calculate_auc <- function(actual, score) {
  roc_curve <- roc(
    response = factor(actual, levels = c("Control", "Parkinson")),
    predictor = score,
    levels = c("Control", "Parkinson"),
    quiet = TRUE
  )
  as.numeric(auc(roc_curve))
}

results <- data.frame()

for (site_count in site_counts) {
  selected_cpgs <- ranked_cpgs[seq_len(site_count)]

  x_train <- train_data[, selected_cpgs, drop = FALSE]
  x_test <- test_data[, selected_cpgs, drop = FALSE]
  y_train <- train_data$disease_state
  y_test <- test_data$disease_state

  preprocess_steps <- preProcess(x_train, method = c("medianImpute", "center", "scale"))
  x_train_scaled <- predict(preprocess_steps, x_train)
  x_test_scaled <- predict(preprocess_steps, x_test)

  y_train_numeric <- ifelse(y_train == "Parkinson", 1, 0)
  logistic_model <- suppressWarnings(glm(
    y_train_numeric ~ .,
    data = data.frame(y_train_numeric = y_train_numeric, x_train_scaled),
    family = "binomial",
    control = glm.control(maxit = 50)
  ))
  logistic_probability <- suppressWarnings(as.numeric(predict(logistic_model, newdata = x_test_scaled, type = "response")))
  if (anyNA(logistic_probability)) {
    logistic_probability[is.na(logistic_probability)] <- mean(logistic_probability, na.rm = TRUE)
  }
  results <- rbind(
    results,
    data.frame(
      Model = "Logistic Regression",
      Sites_Used = site_count,
      AUROC = calculate_auc(y_test, logistic_probability)
    )
  )

  svm_linear <- svm(x_train_scaled, y_train, kernel = "linear", decision.values = TRUE)
  svm_linear_prediction <- predict(svm_linear, x_test_scaled, decision.values = TRUE)
  svm_linear_score <- as.numeric(attr(svm_linear_prediction, "decision.values"))
  results <- rbind(
    results,
    data.frame(
      Model = "SVM (Linear Kernel)",
      Sites_Used = site_count,
      AUROC = calculate_auc(y_test, svm_linear_score)
    )
  )

  svm_radial <- svm(x_train_scaled, y_train, kernel = "radial", decision.values = TRUE)
  svm_radial_prediction <- predict(svm_radial, x_test_scaled, decision.values = TRUE)
  svm_radial_score <- as.numeric(attr(svm_radial_prediction, "decision.values"))
  results <- rbind(
    results,
    data.frame(
      Model = "SVM (Radial Kernel)",
      Sites_Used = site_count,
      AUROC = calculate_auc(y_test, svm_radial_score)
    )
  )
}

results$AUROC <- round(results$AUROC, 4)
write.csv(results, output_file, row.names = FALSE)
print(results)
