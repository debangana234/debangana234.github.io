library(grid)
library(png)

source_file <- "/Users/debangana/Thesis/Parkinson/Rplots/EpigeneticAge/EpivsChronoHorvathGSE111629.png"
output_file <- "masters-thesis/assets/figures/horvath-clock-performance-horvath-dataset.png"

img <- readPNG(source_file)

png(output_file, width = dim(img)[2], height = dim(img)[1], bg = "white")
grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# Cover original facet labels, axis labels, and legend labels before redrawing them in bold.
grid.rect(x = 0.5, y = 0.982, width = 1, height = 0.07, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.295, y = 0.92, width = 0.14, height = 0.035, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.58, y = 0.92, width = 0.16, height = 0.035, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.485, y = 0.032, width = 0.22, height = 0.07, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.035, y = 0.53, width = 0.07, height = 0.42, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.9, y = 0.495, width = 0.22, height = 0.22, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.245, y = 0.935, width = 0.22, height = 0.08, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.66, y = 0.935, width = 0.24, height = 0.08, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.245, y = 0.875, width = 0.18, height = 0.06, gp = gpar(col = NA, fill = "white"))
grid.rect(x = 0.66, y = 0.875, width = 0.19, height = 0.06, gp = gpar(col = NA, fill = "white"))

grid.text(
  "Horvath Age vs. Chronological age (public dataset)",
  x = unit(0.055, "npc"),
  y = unit(0.988, "npc"),
  just = c("left", "top"),
  gp = gpar(fontsize = 30, fontface = "bold", col = "#111820")
)

facet_labels <- data.frame(
  label = c("control", "Parkinson"),
  x = c(0.295, 0.58),
  y = c(0.92, 0.92)
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

r_squared_labels <- data.frame(
  label = c("R² = 0.76", "R² = 0.49"),
  x = c(0.235, 0.645),
  y = c(0.875, 0.875)
)

for (i in seq_len(nrow(r_squared_labels))) {
  grid.text(
    r_squared_labels$label[i],
    x = unit(r_squared_labels$x[i], "npc"),
    y = unit(r_squared_labels$y[i], "npc"),
    just = c("center", "center"),
    gp = gpar(fontsize = 27, fontface = "bold", col = "#111820")
  )
}

grid.text(
  "age",
  x = unit(0.485, "npc"),
  y = unit(0.045, "npc"),
  just = c("center", "bottom"),
  gp = gpar(fontsize = 24, fontface = "bold", col = "#111820")
)

grid.text(
  "Horvath",
  x = unit(0.033, "npc"),
  y = unit(0.53, "npc"),
  rot = 90,
  just = c("center", "bottom"),
  gp = gpar(fontsize = 24, fontface = "bold", col = "#111820")
)

grid.text(
  "stage",
  x = unit(0.86, "npc"),
  y = unit(0.515, "npc"),
  just = c("center", "center"),
  gp = gpar(fontsize = 27, fontface = "bold", col = "#111820")
)

legend_labels <- data.frame(
  label = c("control", "Parkinson"),
  y = c(0.475, 0.438)
)

for (i in seq_len(nrow(legend_labels))) {
  grid.text(
    legend_labels$label[i],
    x = unit(0.86, "npc"),
    y = unit(legend_labels$y[i], "npc"),
    just = c("left", "center"),
    gp = gpar(fontsize = 24, fontface = "bold", col = "#111820")
  )
}

dev.off()
