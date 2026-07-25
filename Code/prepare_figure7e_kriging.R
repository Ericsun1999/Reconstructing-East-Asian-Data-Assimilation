here::i_am("Code/prepare_figure7e_kriging.R")

# ============================================================
# Generate the 12 x 14 kriging intermediate files used by
# Figure 7(e), using the same ordinal-aware kriging equations
# and calibrated parameters as Figure5.R.
#
# Required input:
#   Output/Intermediate/calibration_parameters.rds
#
# Outputs:
#   Data/reaches_kriging_grid12x14_mean.csv
#   Data/reaches_kriging_grid12x14_variance.csv
#
# Important:
#   - The second output stores prediction VARIANCES.
#   - Figure7e.R must take the square root before using these
#     values as standard deviations in the quantile mapping.
# ============================================================

library(here)
library(sp)
library(mvtnorm)

# ------------------------------------------------------------
# 1. Load the calibrated parameters and REACHES observations
# ------------------------------------------------------------

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

calibration_results <- readRDS(
  calibration_file
)

required_objects <- c(
  "temp2",
  "var_fit2",
  "vario_fit2"
)

if (!all(required_objects %in% names(calibration_results))) {
  stop(
    "The calibration RDS file does not contain all required objects: ",
    paste(required_objects, collapse = ", "),
    "."
  )
}

temp2 <- calibration_results$temp2
var.fit2 <- calibration_results$var_fit2
vario.fit2 <- calibration_results$vario_fit2

y2 <- var.fit2$y2

# ------------------------------------------------------------
# 2. Define the years and the 12 x 14 prediction grid
# ------------------------------------------------------------

year_start <- min(temp2$year, na.rm = TRUE)
year_end <- max(temp2$year, na.rm = TRUE)
year_all <- year_start:year_end

ncase <- vapply(
  year_all,
  function(current_year) {
    sum(temp2$year == current_year)
  },
  numeric(1)
)

year2 <- year_all[ncase >= 1]

if (length(year2) == 0L) {
  stop("No years with REACHES observations were found.")
}

loc <- expand.grid(
  long = seq(
    97.5,
    125,
    by = 2.5
  ),
  lat = seq(
    18,
    42.63158,
    by = 1.89473692308
  )
)

coordinates(loc) <- ~ long + lat
proj4string(loc) <- CRS(
  "+proj=longlat +datum=WGS84"
)

n_long <- length(unique(loc@coords[, 1]))
n_lat <- length(unique(loc@coords[, 2]))
n_prediction_locations <- nrow(loc@coords)

if (n_long != 12L ||
    n_lat != 14L ||
    n_prediction_locations != 168L) {
  stop(
    "Expected a 12 x 14 grid with 168 locations, but obtained ",
    n_long,
    " x ",
    n_lat,
    " with ",
    n_prediction_locations,
    " locations."
  )
}

# ------------------------------------------------------------
# 3. Use the same calibrated quantities as Figure5.R
#
# These definitions intentionally follow Figure5.R exactly.
# ------------------------------------------------------------

sigmay <- sqrt(
  vario.fit2$psill2
)

sigmae <- sqrt(
  vario.fit2$psill1
)

alpha <- vario.fit2$range / 100

cuts <- c(
  -Inf,
  -1.5,
  -0.5,
  0.5,
  Inf
)

vals <- c(
  -2,
  -1,
  0,
  1
)

solve_tol <- 1e-12

# ------------------------------------------------------------
# 4. Figure5.R kriging helper functions
# ------------------------------------------------------------

Ez_discrete <- function(
    sigma,
    cuts,
    vals) {

  stopifnot(
    length(vals) ==
      length(cuts) - 1
  )

  a <- head(cuts, -1) / sigma
  b <- tail(cuts, -1) / sigma

  probs <- pnorm(b) - pnorm(a)

  sum(vals * probs)
}

EZstar_h <- function(
    sigma,
    cuts,
    vals) {

  stopifnot(
    length(vals) ==
      length(cuts) - 1
  )

  a <- head(cuts, -1) / sigma
  b <- tail(cuts, -1) / sigma

  part <- sigma * (
    dnorm(a) - dnorm(b)
  )

  sum(vals * part)
}

cov_Z_pair <- function(
    rho,
    sigma,
    cuts,
    vals,
    Eh = NULL) {

  if (is.null(Eh)) {
    Eh <- Ez_discrete(
      sigma,
      cuts,
      vals
    )
  }

  K <- length(vals)

  a_std <- head(cuts, -1) / sigma
  b_std <- tail(cuts, -1) / sigma

  Sigma2 <- matrix(
    c(
      1,
      rho,
      rho,
      1
    ),
    nrow = 2,
    ncol = 2
  )

  Eh2 <- 0

  for (k in seq_len(K)) {
    for (l in seq_len(K)) {

      lower <- c(
        a_std[k],
        a_std[l]
      )

      upper <- c(
        b_std[k],
        b_std[l]
      )

      pij <- as.numeric(
        pmvnorm(
          lower = lower,
          upper = upper,
          mean = c(0, 0),
          sigma = Sigma2
        )
      )

      Eh2 <- Eh2 +
        vals[k] * vals[l] * pij
    }
  }

  Eh2 - Eh^2
}

cZY_vector <- function(
    s_coords,
    s0,
    sigma_Y2,
    sigma_E2,
    alpha,
    cuts = c(
      -Inf,
      -1.5,
      -0.5,
      0.5,
      Inf
    ),
    vals = c(
      -2,
      -1,
      0,
      1
    )) {

  sigma2 <- sigma_Y2 + sigma_E2
  sigma <- sqrt(sigma2)

  Ezstar_h <- EZstar_h(
    sigma,
    cuts,
    vals
  )

  s0_matrix <- matrix(
    s0,
    nrow = nrow(s_coords),
    ncol = ncol(s_coords),
    byrow = TRUE
  )

  dists <- sqrt(
    rowSums(
      (s_coords - s0_matrix)^2
    )
  )

  covYY <- sigma_Y2 *
    exp(-dists / alpha)

  (covYY / sigma2) * Ezstar_h
}

SigmaZ_matrix <- function(
    s_coords,
    sigma_Y2,
    sigma_E2,
    alpha,
    cuts = c(
      -Inf,
      -1.5,
      -0.5,
      0.5,
      Inf
    ),
    vals = c(
      -2,
      -1,
      0,
      1
    )) {

  n <- nrow(s_coords)
  sigma2 <- sigma_Y2 + sigma_E2
  sigma <- sqrt(sigma2)

  Eh <- Ez_discrete(
    sigma,
    cuts,
    vals
  )

  dmat <- as.matrix(
    dist(
      s_coords,
      method = "euclidean",
      upper = TRUE,
      diag = TRUE
    )
  )

  rho <- (
    sigma_Y2 *
      exp(-dmat / alpha)
  ) / sigma2

  diag(rho) <- 1

  Sig <- matrix(
    NA_real_,
    nrow = n,
    ncol = n
  )

  for (i in seq_len(n)) {
    for (j in i:n) {

      cij <- cov_Z_pair(
        rho = rho[i, j],
        sigma = sigma,
        cuts = cuts,
        vals = vals,
        Eh = Eh
      )

      Sig[i, j] <- cij

      if (j != i) {
        Sig[j, i] <- cij
      }
    }
  }

  Sig
}

# ------------------------------------------------------------
# 5. Predict one year
#
# This uses the same equations as Figure5.R. The only
# computational simplification is that SigmaZ is constructed
# once per year rather than once per prediction location.
# That does not change the statistical calculation.
# ------------------------------------------------------------

predict_one_year <- function(
    current_year,
    observation_data,
    prediction_coordinates,
    sigma_Y2,
    sigma_E2,
    alpha,
    cuts,
    vals,
    tol) {

  # ----------------------------------------------------------
  # 1. Extract observations for this year
  # ----------------------------------------------------------

  temp14 <- observation_data[
    observation_data$year == current_year,
  ]

  if (nrow(temp14) == 0L) {
    stop(
      "No REACHES observations were found for year ",
      current_year,
      "."
    )
  }

  locations <- temp14@coords
  z_obs <- as.numeric(temp14$level)

  if (length(z_obs) != nrow(locations)) {
    stop(
      "The observation values and coordinates have different ",
      "lengths for year ",
      current_year,
      "."
    )
  }

  # ----------------------------------------------------------
  # 2. Construct SigmaZ once for this year
  # ----------------------------------------------------------

  SigmaZ <- SigmaZ_matrix(
    s_coords = locations,
    sigma_Y2 = sigma_Y2,
    sigma_E2 = sigma_E2,
    alpha = alpha,
    cuts = cuts,
    vals = vals
  )

  SigmaZ <- (
    SigmaZ + t(SigmaZ)
  ) / 2

  # ----------------------------------------------------------
  # 3. Construct the right-hand side for the prediction means
  # ----------------------------------------------------------

  sigma2_total <- sigma_Y2 + sigma_E2
  sigma_total <- sqrt(sigma2_total)

  Ez <- Ez_discrete(
    sigma = sigma_total,
    cuts = cuts,
    vals = vals
  )

  rhs_mean <- z_obs - Ez

  # ----------------------------------------------------------
  # 4. Construct all 168 cross-covariance vectors
  #
  # CZY has:
  #   rows    = observation locations
  #   columns = prediction locations
  # ----------------------------------------------------------

  CZY <- vapply(
    seq_len(nrow(prediction_coordinates)),
    function(j) {

      s0 <- as.numeric(
        prediction_coordinates[j, ]
      )

      cZY_vector(
        s_coords = locations,
        s0 = s0,
        sigma_Y2 = sigma_Y2,
        sigma_E2 = sigma_E2,
        alpha = alpha,
        cuts = cuts,
        vals = vals
      )
    },
    numeric(nrow(locations))
  )

  if (!is.matrix(CZY)) {
    CZY <- matrix(
      CZY,
      nrow = nrow(locations)
    )
  }

  # ----------------------------------------------------------
  # 5. Solve all systems at once
  #
  # Column 1:
  #   SigmaZ^{-1} (z_obs - Ez)
  #
  # Columns 2--169:
  #   SigmaZ^{-1} c_zy(s0)
  # ----------------------------------------------------------

  all_rhs <- cbind(
    rhs_mean,
    CZY
  )

  all_solutions <- tryCatch(
    solve(
      SigmaZ,
      all_rhs,
      tol = tol
    ),
    error = function(e) {
      stop(
        "Failed to solve the kriging systems for year ",
        current_year,
        ". Reciprocal condition number = ",
        signif(rcond(SigmaZ), 6),
        ". Original error: ",
        conditionMessage(e)
      )
    }
  )

  w_mean <- all_solutions[, 1]

  W_variance <- all_solutions[
    ,
    -1,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # 6. Prediction means
  #
  # For each prediction location:
  #   c_zy' SigmaZ^{-1} (z_obs - Ez)
  # ----------------------------------------------------------

  prediction_mean <- drop(
    crossprod(
      CZY,
      w_mean
    )
  )

  # ----------------------------------------------------------
  # 7. Prediction variances
  #
  # For each prediction location:
  #   sigma_Y2 - c_zy' SigmaZ^{-1} c_zy
  # ----------------------------------------------------------

  covariance_reduction <- colSums(
    CZY * W_variance
  )

  prediction_variance <- sigma_Y2 -
    covariance_reduction

  # Avoid negligible negative values caused by floating-point
  # numerical error.
  prediction_variance <- pmax(
    prediction_variance,
    0
  )

  list(
    mean = prediction_mean,
    variance = prediction_variance
  )
}


# ------------------------------------------------------------
# 6. Generate all yearly predictions
# ------------------------------------------------------------

mean_matrix <- matrix(
  NA_real_,
  nrow = n_prediction_locations,
  ncol = length(year2)
)

variance_matrix <- matrix(
  NA_real_,
  nrow = n_prediction_locations,
  ncol = length(year2)
)

colnames(mean_matrix) <- as.character(year2)
colnames(variance_matrix) <- as.character(year2)

# pmvnorm can use randomized numerical integration.
# Fix the random sequence for reproducibility.
set.seed(10)

start_time <- Sys.time()

for (i in seq_along(year2)) {

  current_year <- year2[i]

  message(
    "Kriging year ",
    current_year,
    " (",
    i,
    "/",
    length(year2),
    ")"
  )

  yearly_result <- predict_one_year(
    current_year = current_year,
    observation_data = y2,
    prediction_coordinates = loc@coords,
    sigma_Y2 = sigmay,
    sigma_E2 = sigmae,
    alpha = alpha,
    cuts = cuts,
    vals = vals,
    tol = solve_tol
  )

  mean_matrix[, i] <- yearly_result$mean
  variance_matrix[, i] <- yearly_result$variance
}

elapsed_time <- difftime(
  Sys.time(),
  start_time,
  units = "mins"
)

# ------------------------------------------------------------
# 7. Save the two intermediate CSV files
# ------------------------------------------------------------

mean_output <- data.frame(
  long = loc@coords[, 1],
  lat = loc@coords[, 2],
  mean_matrix,
  check.names = FALSE
)

variance_output <- data.frame(
  long = loc@coords[, 1],
  lat = loc@coords[, 2],
  variance_matrix,
  check.names = FALSE
)

mean_output_file <- here::here(
  "Data",
  "reaches_kriging_grid12x14_mean.csv"
)

variance_output_file <- here::here(
  "Data",
  "reaches_kriging_grid12x14_variance.csv"
)

write.csv(
  mean_output,
  mean_output_file,
  row.names = FALSE
)

write.csv(
  variance_output,
  variance_output_file,
  row.names = FALSE
)

message(
  "12 x 14 kriging means saved to: ",
  mean_output_file
)

message(
  "12 x 14 prediction variances saved to: ",
  variance_output_file
)

message(
  "Elapsed time: ",
  round(
    as.numeric(elapsed_time),
    2
  ),
  " minutes."
)
