here::i_am("Code/Figure4.R")

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

var.fit2 <- calibration_results$var_fit2
vario.fit2 <- calibration_results$vario_fit2
