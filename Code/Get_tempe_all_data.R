here::i_am("Code/Get_tempe_all_data.R")

# ============================================================
# Generate REACHES kriging results for the complete 53 x 49
# grid and extract the three locations used in Figures 7--10.
#
# The kriging equations and calibrated parameter definitions
# are intentionally the same as those used in Figure5.R.
#
# Required input:
#   Output/Intermediate/calibration_parameters.rds
#
# Outputs:
#   Data/reaches_kriging_grid53x49_mean.csv
#   Data/reaches_kriging_grid53x49_variance.csv
#   Data/reaches_kriging_city3_mean.csv
#   Data/reaches_kriging_city3_sd.csv
#
# City-row order in the city3 files:
#   1. Hong Kong
#   2. Shanghai
#   3. Beijing
#
# Note:
#   The full-grid uncertainty file stores prediction VARIANCE.
#   The city3 uncertainty file stores prediction STANDARD
#   DEVIATION for compatibility with Figure7-8.R.
# ============================================================

library(sp)
library(mvtnorm)

# ------------------------------------------------------------
# 1. Load shared calibration results
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
    "The calibration file does not contain all required objects: ",
    paste(required_objects, collapse = ", "),
    "."
  )
}

temp2 <- calibration_results$temp2
var.fit2 <- calibration_results$var_fit2
vario.fit2 <- calibration_results$vario_fit2
y2 <- var.fit2$y2

# ------------------------------------------------------------
# 2. Define all available years and the complete 53 x 49 grid
# ------------------------------------------------------------

year2 <- sort(
  unique(
    as.integer(temp2$year)
  )
)

year2 <- year2[
  is.finite(year2)
]

if (length(year2) == 0L) {
  stop("No valid REACHES observation years were found.")
}

loc <- expand.grid(
  long = seq(
    98.25,
    124.25,
    by = 0.5
  ),
  lat = seq(
    18.25,
    42.25,
    by = 0.5
  )
)

coordinates(loc) <- ~ long + lat
proj4string(loc) <- CRS(
  "+proj=longlat +datum=WGS84"
)

n_long <- length(
  unique(loc@coords[, 1])
)

n_lat <- length(
  unique(loc@coords[, 2])
)

n_grid <- nrow(loc@coords)

if (n_long != 53L ||
    n_lat != 49L ||
    n_grid != 2597L) {
  stop(
    "Expected a 53 x 49 grid with 2597 locations, but obtained ",
    n_long,
    " x ",
    n_lat,
    " with ",
    n_grid,
    " locations."
  )
}

# ------------------------------------------------------------
# 3. Kriging functions copied from Figure5.R
# ------------------------------------------------------------

# These definitions preserve the original calculation used in
# Figure 5. The square-root transformation is intentionally kept
# because it is part of the original analysis code.
sigmay <- sqrt(vario.fit2$psill2)
sigmae <- sqrt(vario.fit2$psill1)
alpha <- vario.fit2$range / 100

cuts <- c(-Inf, -1.5, -0.5, 0.5, Inf)
vals <- c(-2, -1, 0, 1)

Ez_discrete <- function(sigma, cuts, vals) {
  stopifnot(length(vals) == length(cuts) - 1)

  a <- head(cuts, -1) / sigma
  b <- tail(cuts, -1) / sigma
  probs <- pnorm(b) - pnorm(a)

  sum(vals * probs)
}

EZstar_h <- function(sigma, cuts, vals) {
  stopifnot(length(vals) == length(cuts) - 1)

  a <- head(cuts, -1) / sigma
  b <- tail(cuts, -1) / sigma
  part <- sigma * (dnorm(a) - dnorm(b))

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
    c(1, rho, rho, 1),
    2,
    2
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
    cuts = c(-Inf, -1.5, -0.5, 0.5, Inf),
    vals = c(-2, -1, 0, 1)) {

  sigma2 <- sigma_Y2 + sigma_E2
  sigma <- sqrt(sigma2)

  Ezstar_h <- EZstar_h(
    sigma,
    cuts,
    vals
  )

  dists <- sqrt(
    rowSums(
      (
        s_coords -
          matrix(
            s0,
            nrow(s_coords),
            ncol(s_coords),
            byrow = TRUE
          )
      )^2
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
    cuts = c(-Inf, -1.5, -0.5, 0.5, Inf),
    vals = c(-2, -1, 0, 1)) {

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
    n,
    n
  )

  for (i in seq_len(n)) {
    for (j in i:n) {
      cij <- cov_Z_pair(
        rho[i, j],
        sigma,
        cuts,
        vals,
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

predict_grid_one_year <- function(
    current_year,
    observation_data,
    prediction_coordinates,
    sigma_Y2,
    sigma_E2,
    alpha,
    cuts = c(-Inf, -1.5, -0.5, 0.5, Inf),
    vals = c(-2, -1, 0, 1),
    tol = 1e-2) {

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

  # SigmaZ depends on the observation locations and model
  # parameters, not on the prediction location. Construct it
  # only once for this year.
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

  sigma2 <- sigma_Y2 + sigma_E2
  sigma <- sqrt(sigma2)

  Ez <- Ez_discrete(
    sigma = sigma,
    cuts = cuts,
    vals = vals
  )

  rhs_mean <- z_obs - Ez

  # Each column is the cross-covariance vector for one
  # prediction grid point.
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

  # Solve the mean system and all prediction-variance systems
  # together. This performs one matrix factorization per year.
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

  prediction_mean <- drop(
    crossprod(
      CZY,
      w_mean
    )
  )

  prediction_variance <- sigma_Y2 -
    colSums(
      CZY * W_variance
    )

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
# 4. Generate predictions for every grid location and year
# ------------------------------------------------------------

mean_matrix <- matrix(
  NA_real_,
  nrow = n_grid,
  ncol = length(year2)
)

variance_matrix <- matrix(
  NA_real_,
  nrow = n_grid,
  ncol = length(year2)
)

colnames(mean_matrix) <- as.character(year2)
colnames(variance_matrix) <- as.character(year2)

start_time <- Sys.time()

for (i in seq_along(year2)) {

  current_year <- year2[i]

  # pmvnorm() uses numerical integration. Give each year a
  # fixed random-number state so the result does not depend on
  # which other years were run before it.
  set.seed(
    500000L + current_year
  )

  message(
    "Kriging year ",
    current_year,
    " (",
    i,
    "/",
    length(year2),
    ")"
  )

  current_data <- y2[
    y2$year == current_year,
  ]

  duplicated_coordinates <- duplicated(
    as.data.frame(
      current_data@coords
    )
  )

  if (any(duplicated_coordinates)) {
    warning(
      "Year ",
      current_year,
      " contains ",
      sum(duplicated_coordinates),
      " duplicated observation coordinates. ",
      "The records are retained in this version."
    )
  }

  yearly_result <- predict_grid_one_year(
    current_year = current_year,
    observation_data = y2,
    prediction_coordinates = loc@coords,
    sigma_Y2 = sigmay,
    sigma_E2 = sigmae,
    alpha = alpha,
    cuts = cuts,
    vals = vals,
    tol = 1e-2
  )

  mean_matrix[, i] <- yearly_result$mean
  variance_matrix[, i] <- yearly_result$variance
}

elapsed_minutes <- as.numeric(
  difftime(
    Sys.time(),
    start_time,
    units = "mins"
  )
)

# ------------------------------------------------------------
# 5. Assemble and save complete-grid outputs
# ------------------------------------------------------------

grid_coordinates <- as.data.frame(
  loc@coords
)

grid_mean <- cbind(
  grid_coordinates,
  as.data.frame(
    mean_matrix,
    check.names = FALSE
  )
)

grid_variance <- cbind(
  grid_coordinates,
  as.data.frame(
    variance_matrix,
    check.names = FALSE
  )
)

names(grid_mean)[1:2] <- c(
  "long",
  "lat"
)

names(grid_variance)[1:2] <- c(
  "long",
  "lat"
)

grid_mean_file <- here::here(
  "Data",
  "reaches_kriging_grid53x49_mean.csv"
)

grid_variance_file <- here::here(
  "Data",
  "reaches_kriging_grid53x49_variance.csv"
)

write.csv(
  grid_mean,
  grid_mean_file,
  row.names = FALSE
)

write.csv(
  grid_variance,
  grid_variance_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Extract Hong Kong, Shanghai, and Beijing by coordinates
#
# These coordinates reproduce the original rows:
#   Hong Kong: row 456
#   Shanghai:  row 1425
#   Beijing:   row 2316
#
# Using coordinates is safer than depending on row numbers.
# ------------------------------------------------------------

city_locations <- data.frame(
  city = c(
    "HongKong",
    "Shanghai",
    "Beijing"
  ),
  long = c(
    113.75,
    121.25,
    116.25
  ),
  lat = c(
    22.25,
    31.25,
    39.75
  )
)

coordinate_key <- function(
    long,
    lat) {
  sprintf(
    "%.2f_%.2f",
    long,
    lat
  )
}

grid_keys <- coordinate_key(
  grid_mean$long,
  grid_mean$lat
)

city_keys <- coordinate_key(
  city_locations$long,
  city_locations$lat
)

city_indices <- match(
  city_keys,
  grid_keys
)

if (anyNA(city_indices)) {
  stop(
    "At least one city location could not be matched to the ",
    "53 x 49 prediction grid."
  )
}

city3_mean <- grid_mean[
  city_indices,
  ,
  drop = FALSE
]

city3_variance <- grid_variance[
  city_indices,
  ,
  drop = FALSE
]

city3_sd <- city3_variance

year_columns <- names(city3_sd)[
  -c(1, 2)
]

city3_sd[
  ,
  year_columns
] <- sqrt(
  pmax(
    as.matrix(
      city3_variance[
        ,
        year_columns,
        drop = FALSE
      ]
    ),
    0
  )
)

city3_mean_file <- here::here(
  "Data",
  "reaches_kriging_city3_mean.csv"
)

city3_sd_file <- here::here(
  "Data",
  "reaches_kriging_city3_sd.csv"
)

write.csv(
  city3_mean,
  city3_mean_file,
  row.names = FALSE
)

write.csv(
  city3_sd,
  city3_sd_file,
  row.names = FALSE
)

message(
  "Complete-grid kriging means saved to: ",
  grid_mean_file
)

message(
  "Complete-grid prediction variances saved to: ",
  grid_variance_file
)

message(
  "Three-city kriging means saved to: ",
  city3_mean_file
)

message(
  "Three-city prediction standard deviations saved to: ",
  city3_sd_file
)

message(
  "Elapsed time: ",
  round(elapsed_minutes, 2),
  " minutes."
)
