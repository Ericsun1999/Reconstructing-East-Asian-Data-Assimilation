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

figure4_output_dir <- here::here(
  "Output",
  "Figure4"
)

dir.create(
  figure4_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = file.path(
    figure4_output_dir,
    "Figure4(a).png"
  ),
  plot = vario.fit2$plot1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(
    figure4_output_dir,
    "Figure4(b).png"
  ),
  plot = vario.fit2$plot2,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(
    figure4_output_dir,
    "Figure4(c).png"
  ),
  plot = vario.fit2$plot3,
  width = 5,
  height = 6,
  units = "in",
  dpi = 300
)

message(
  "All Figure 4 panels were saved to: ",
  figure4_output_dir
)
