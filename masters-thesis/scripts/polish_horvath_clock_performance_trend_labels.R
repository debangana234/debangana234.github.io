library(grid)
library(png)

source_file <- "/Users/debangana/Thesis/Parkinson/Rplots/EpigeneticAge/Horvath_vs_ChoronoAgeTrend.png"
output_file <- "masters-thesis/assets/figures/horvath-clock-performance-trend.png"

img <- readPNG(source_file)

png(output_file, width = dim(img)[2], height = dim(img)[1], bg = "white")
grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# Cover original title, facet-strip labels, x-axis label, and legend labels.
grid.rect(x = 0.5, y = 0.982, width = 1, height = 0.055, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.43, y = 0.938, width = 0.68, height = 0.032, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.43, y = 0.492, width = 0.68, height = 0.058, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.45, y = 0.022, width = 0.14, height = 0.055, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.9, y = 0.455, width = 0.2, height = 0.26, gp = gpar(col = NA, fill = "white"))

grid.text(
  "Horvath Age vs Chronological Age (TREND dataset)",
  x = unit(0.063, "npc"),
  y = unit(0.988, "npc"),
  just = c("left", "top"),
  gp = gpar(fontsize = 28, fontface = "bold", col = "#111820")
)

facet_labels <- data.frame(
  label = c("early_control", "late_control", "PreConverter", "PostConverter"),
  x = c(0.25, 0.61, 0.25, 0.61),
  y = c(0.938, 0.938, 0.495, 0.495)
)

for (i in seq_len(nrow(facet_labels))) {
  grid.text(
    facet_labels$label[i],
    x = unit(facet_labels$x[i], "npc"),
    y = unit(facet_labels$y[i], "npc"),
    just = c("center", "center"),
    gp = gpar(fontsize = 25, fontface = "bold", col = "#111820")
  )
}

grid.text(
  "age",
  x = unit(0.45, "npc"),
  y = unit(0.025, "npc"),
  just = c("center", "bottom"),
  gp = gpar(fontsize = 24, fontface = "bold", col = "#111820")
)

grid.text(
  "stage",
  x = unit(0.872, "npc"),
  y = unit(0.535, "npc"),
  just = c("center", "center"),
  gp = gpar(fontsize = 27, fontface = "bold", col = "#111820")
)

legend_labels <- data.frame(
  label = c("early_control", "late_control", "PreConverter", "PostConverter"),
  y = c(0.493, 0.456, 0.419, 0.382)
)

for (i in seq_len(nrow(legend_labels))) {
  grid.text(
    legend_labels$label[i],
    x = unit(0.872, "npc"),
    y = unit(legend_labels$y[i], "npc"),
    just = c("left", "center"),
    gp = gpar(fontsize = 22, fontface = "bold", col = "#111820")
  )
}

dev.off()
