here::i_am("Code/Figure4.R")

# ============================================================
# Calibration functions for Figure 4
#
# Required input:
#   Output/Intermediate/calibration_parameters.rds
#
# Outputs:
#   Output/Figure4/Figure4(a).png
#   Output/Figure4/Figure4(b).png
#   Output/Figure4/Figure4(c).png
# ============================================================

library(ggplot2)

calibration_file <- here::here(
  "Output",
  "Intermediate",
  "calibration_parameters.rds"
)

if (!file.exists(calibration_file)) {
  stop(
    "Required calibration file was not found: ",
    calibration_file,
    "\nRun Code/prepare_calibration.R first."
  )
}

calibration_results <- readRDS(calibration_file)
vario.fit2 <- calibration_results$vario_fit2

required_plots <- c("plot1", "plot2", "plot3")

if (!all(required_plots %in% names(vario.fit2))) {
  stop(
    "The calibration file does not contain all Figure 4 plots. ",
    "Run Code/prepare_calibration.R again."
  )
}

figure4_output_dir <- here::here("Output", "Figure4")

dir.create(
  figure4_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Styling helper
# ------------------------------------------------------------
style_figure4_plot <- function(
    p,
    base_size = 22) {

  p +
    theme_gray(
      base_size = base_size
    ) +
    theme(
      text = element_text(
        size = base_size
      ),
      axis.title = element_text(
        size = base_size
      ),
      axis.text = element_text(
        size = base_size * 0.8,
        colour = "black"
      ),
      legend.position = "none",
      panel.grid.major = element_line(
        colour = "white",
        linewidth = 0.7
      ),
      panel.grid.minor = element_line(
        colour = "white",
        linewidth = 0.35
      )
    )
}

# ------------------------------------------------------------
# Restyle plots
# ------------------------------------------------------------

# Figure 4(a): zoom x-axis to 1.4
p4a <- style_figure4_plot(
  vario.fit2$plot1 +
    coord_cartesian(
      xlim = c(0, 1.4)
    ) +
    scale_x_continuous(
      breaks = c(
        0,
        0.5,
        1.0
      )
    ),
  base_size = 22
)

p4b <- style_figure4_plot(
  vario.fit2$plot2,
  base_size = 22
)

p4c <- style_figure4_plot(
  vario.fit2$plot3,
  base_size = 22
)

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
ggsave(
  filename = file.path(figure4_output_dir, "Figure4(a).png"),
  plot = p4a,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(figure4_output_dir, "Figure4(b).png"),
  plot = p4b,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(figure4_output_dir, "Figure4(c).png"),
  plot = p4c,
  width = 5,
  height = 6,
  units = "in",
  dpi = 300
)

message(
  "All Figure 4 panels were saved to: ",
  figure4_output_dir
)
