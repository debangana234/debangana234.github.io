library(grid)
library(png)

source_file <- "/Users/debangana/Thesis/Parkinson/Rplots/EpigeneticAge/EAAhorvathTrend.png"
output_file <- "masters-thesis/assets/figures/epigenetic-age-acceleration-stage-boxplot.png"

img <- readPNG(source_file)

png(output_file, width = dim(img)[2], height = dim(img)[1], bg = "white")
grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# Cover only the original text regions, leaving the plotted data untouched.
grid.rect(x = 0.5, y = 0.975, width = 1, height = 0.05, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.5, y = 0.03, width = 1, height = 0.075, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.045, y = 0.5, width = 0.09, height = 0.52, gp = gpar(col = NA, fill = "white"))

grid.text(
  "Epigenetic Age Acceleration by Stage",
  x = unit(0.07, "npc"),
  y = unit(0.983, "npc"),
  just = c("left", "top"),
  gp = gpar(fontsize = 40, fontface = "bold", col = "#111820")
)

grid.text(
  "Stage",
  x = unit(0.535, "npc"),
  y = unit(0.012, "npc"),
  just = c("center", "bottom"),
  gp = gpar(fontsize = 27, fontface = "bold", col = "#111820")
)

stage_labels <- c("early_control", "late_control", "PreConverter", "PostConverter")
stage_x <- c(0.2, 0.42, 0.64, 0.86)
for (i in seq_along(stage_labels)) {
  grid.text(
    stage_labels[i],
    x = unit(stage_x[i], "npc"),
    y = unit(0.053, "npc"),
    just = c("center", "bottom"),
    gp = gpar(fontsize = 29, fontface = "bold", col = "#3b414a")
  )
}

grid.text(
  "EAA (Horvath)",
  x = unit(0.018, "npc"),
  y = unit(0.5, "npc"),
  rot = 90,
  just = c("center", "bottom"),
  gp = gpar(fontsize = 30, fontface = "bold", col = "#111820")
)

dev.off()
