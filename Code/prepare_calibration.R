here::i_am("Code/prepare_calibration.R")

# ============================================================
# Shared parameter estimation used by Figures 4, 5, and 7(e)
#
# Output:
#   Output/Intermediate/calibration_parameters.rds
# ============================================================

library(readxl)
library(dplyr)
library(sp)
library(spacetime)
library(zoo)
library(gstat)
library(mgcv)
library(ggplot2)
library(mvtnorm)
library(MASS)

# ------------------------------------------------------------
# 1. Load REACHES data
# ------------------------------------------------------------

temperature <- read_excel(
  here::here("Data", "temperature index value.v1.xlsx"),
  col_types = c(
    "skip", "skip", "numeric", "numeric",
    "skip", "skip", "skip", "skip",
    "skip", "numeric", "numeric",
    "skip", "skip"
  )
)

colnames(temperature) <- c(
  "level",
  "year",
  "long",
  "lat"
)

temp2 <- temperature %>%
  group_by(year, long, lat)

# ------------------------------------------------------------
# 2. Shared functions
# ------------------------------------------------------------

par_est_initial <- function(temp2) {

  # 把 Figure4.R 裡原本完整的函數內容放在這裡，
  # 暫時不要改動其計算方式。

}

findgam <- function(y, z, y_target) {

  # 把 Figure4.R 裡原本完整的函數內容放在這裡。

}

par_estimation_calibrate <- function(
    psill1,
    psill,
    range,
    n1 = 500,
    n2 = 200,
    n3 = 40,
    plot = FALSE) {

  # 把整理後的完整函數內容放在這裡。
  #
  # 最後的 return list 應保留：
  # psill1, psill2, range, plot1, plot2, plot3

}

# ------------------------------------------------------------
# 3. Estimate shared parameters
# ------------------------------------------------------------

var.fit2 <- par_est_initial(temp2)

set.seed(10)

vario.fit2 <- par_estimation_calibrate(
  psill1 = var.fit2$psill[1],
  psill = var.fit2$psill,
  range = var.fit2$range,
  n1 = 500,
  n2 = 200,
  n3 = 40,
  plot = FALSE
)

# ------------------------------------------------------------
# 4. Save intermediate results
# ------------------------------------------------------------

intermediate_dir <- here::here(
  "Output",
  "Intermediate"
)

dir.create(
  intermediate_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

calibration_file <- file.path(
  intermediate_dir,
  "calibration_parameters.rds"
)

saveRDS(
  list(
    temp2 = temp2,
    var_fit2 = var.fit2,
    vario_fit2 = vario.fit2
  ),
  file = calibration_file
)

message(
  "Shared calibration parameters saved to: ",
  calibration_file
)
