library(methylclockData)

metadata_file <- "/Users/debangana/Thesis/Parkinson/annotation/metadataUpdated.csv"
beta_file <- "/Users/debangana/Thesis/Parkinson/annotation/beta_values.csv"
output_file <- "masters-thesis/assets/data/trend_horvath_age_estimates.csv"

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

metadata <- read.csv(metadata_file, check.names = FALSE)
metadata <- metadata[!grepl("^7004_", metadata$Sample_Name), ]

coef_horvath <- get_coefHorvath()
horvath_cpgs <- setdiff(coef_horvath$CpGmarker, "(Intercept)")

header <- strsplit(readLines(beta_file, n = 1), ",", fixed = TRUE)[[1]]
header <- gsub('^"|"$', "", header)
sample_columns <- header[-1]
missing_samples <- setdiff(metadata$Sample_Name, sample_columns)
if (length(missing_samples) > 0) {
  stop("Missing samples in beta matrix: ", paste(missing_samples, collapse = ", "))
}

selected_lines <- c(readLines(beta_file, n = 1))
con <- file(beta_file, open = "r")
on.exit(close(con), add = TRUE)
invisible(readLines(con, n = 1))
repeat {
  lines <- readLines(con, n = 5000)
  if (length(lines) == 0) break
  ids <- sub('^"?([^",]+)"?,.*$', "\\1", lines)
  selected_lines <- c(selected_lines, lines[ids %in% horvath_cpgs])
  if (sum(ids %in% horvath_cpgs) > 0 && length(intersect(sub('^"?([^",]+)"?,.*$', "\\1", selected_lines[-1]), horvath_cpgs)) == length(horvath_cpgs)) {
    break
  }
}

beta_df <- read.csv(text = paste(selected_lines, collapse = "\n"), row.names = 1, check.names = FALSE)
observed_cpgs <- intersect(horvath_cpgs, rownames(beta_df))
if (length(observed_cpgs) / length(horvath_cpgs) <= 0.8) {
  stop("Too few Horvath CpGs found in beta matrix.")
}

beta_df <- beta_df[observed_cpgs, metadata$Sample_Name, drop = FALSE]
coef_subset <- coef_horvath[match(observed_cpgs, coef_horvath$CpGmarker), ]
linear_age <- as.numeric(t(as.matrix(beta_df)) %*% coef_subset$CoefficientTraining) +
  coef_horvath$CoefficientTraining[coef_horvath$CpGmarker == "(Intercept)"]
horvath_age <- ifelse(linear_age < 0, 21 * exp(linear_age) - 1, 21 * linear_age + 20)

trend_age <- data.frame(
  Sample_Name = metadata$Sample_Name,
  TREND.ID = metadata$TREND.ID,
  age = metadata$age,
  sex = metadata$sex,
  stage = metadata$stage,
  group = metadata$group,
  Horvath = horvath_age,
  stringsAsFactors = FALSE
)

write.csv(trend_age, output_file, row.names = FALSE)
