here::i_am("Code/Get_tempe_all_data.R")

# ============================================================
# Generate interval-censored REACHES kriging products for:
#   1. the complete 53 x 49 REACHES grid;
#   2. the 266 native LME grid locations used in Figure 6;
#   3. the three locations used in Figures 7--10.
#
# Required input:
#   Output/Intermediate/calibration_parameters.rds
#
# Outputs:
#   Output/Intermediate/REACHES/
#     reaches_kriging_grid53x49_mean.csv
#     reaches_kriging_grid53x49_variance.csv
#     reaches_kriging_lme_grid_mean.csv
#     reaches_kriging_lme_grid_variance.csv
#     reaches_kriging_city3_mean.csv
#     reaches_kriging_city3_sd.csv
#     reaches_kriging_metadata.csv
#
# City-row order:
#   1. Hong Kong
#   2. Shanghai
#   3. Beijing
#
# Notes:
#   - The complete-grid uncertainty file stores prediction
#     VARIANCE.
#   - The city uncertainty file stores prediction STANDARD
#     DEVIATION.
#   - Multiple documentary records at the same site-year are
#     retained as distinct noisy measurements of the same latent
#     spatial process value.
# ============================================================

library(sp)
library(mvtnorm)
library(readr)

# ------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------

checkpoint_every <- 25L
resume_from_checkpoint <- TRUE
rho_lookup_size <- 4001L

calibration_file <- here::here(
  "Output",
  "Intermediate",
  "calibration_parameters.rds"
)

lme_annual_file <- here::here(
  "Output",
  "Intermediate",
  "LME",
  "lme_annual_1368_1911.rds"
)

output_directory <- here::here(
  "Output",
  "Intermediate",
  "REACHES"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

checkpoint_file <- file.path(
  output_directory,
  "reaches_kriging_checkpoint.rds"
)

# ------------------------------------------------------------
# 2. Load calibration results
# ------------------------------------------------------------

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

if (!all(
  required_objects %in%
    names(calibration_results)
)) {
  stop(
    "The calibration file does not contain all required ",
    "objects: ",
    paste(
      required_objects,
      collapse = ", "
    ),
    "."
  )
}

temp2 <- calibration_results$temp2
var.fit2 <- calibration_results$var_fit2
vario.fit2 <- calibration_results$vario_fit2
y2 <- var.fit2$y2

# The calibrated psill values are variances. Do not take their
# square roots before passing them to covariance functions.
sigma_Y2 <- vario.fit2$psill2
sigma_E2 <- vario.fit2$psill1

# The pooled variogram range is on the kilometre scale because
# the observations use longitude/latitude coordinates. All
# prediction distances below are therefore great-circle
# distances in kilometres.
alpha_km <- vario.fit2$range

parameter_vector <- c(
  sigma_Y2 = sigma_Y2,
  sigma_E2 = sigma_E2,
  alpha_km = alpha_km
)

if (
  any(!is.finite(parameter_vector)) ||
    sigma_Y2 <= 0 ||
    sigma_E2 < 0 ||
    alpha_km <= 0
) {
  stop(
    "Invalid calibrated covariance parameters were read from: ",
    calibration_file
  )
}

# ------------------------------------------------------------
# 3. Define event years and the complete 53 x 49 grid
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
  stop(
    "No valid REACHES event years were found."
  )
}

if (length(year2) != 524L) {
  warning(
    "The manuscript describes 524 event years, but ",
    length(year2),
    " were found."
  )
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

n_grid <- nrow(
  loc@coords
)

if (
  n_long != 53L ||
    n_lat != 49L ||
    n_grid != 2597L
) {
  stop(
    "Expected a 53 x 49 grid with 2597 locations, ",
    "but obtained ",
    n_long,
    " x ",
    n_lat,
    " with ",
    n_grid,
    " locations."
  )
}

# ------------------------------------------------------------
# 4. Load the native LME grid used in Figure 6
# ------------------------------------------------------------

if (!file.exists(lme_annual_file)) {
  stop(
    "The annual LME archive was not found: ",
    lme_annual_file,
    "\nRun Code/DataPreparation/prepare_lme_annual.R first."
  )
}

lme_annual_archive <- readRDS(
  lme_annual_file
)

if (
  !"coordinates" %in%
    names(lme_annual_archive)
) {
  stop(
    "The annual LME archive does not contain coordinates."
  )
}

lme_coordinates <- lme_annual_archive$coordinates

required_lme_coordinate_columns <- c(
  "location_id",
  "lati",
  "long"
)

if (!all(
  required_lme_coordinate_columns %in%
    names(lme_coordinates)
)) {
  stop(
    "The annual LME coordinate table must contain: ",
    paste(
      required_lme_coordinate_columns,
      collapse = ", "
    ),
    "."
  )
}

lme_coordinates <- lme_coordinates[
  ,
  required_lme_coordinate_columns
]

lme_coordinates$location_id <- as.integer(
  lme_coordinates$location_id
)

lme_coordinates$lati <- as.numeric(
  lme_coordinates$lati
)

lme_coordinates$long <- as.numeric(
  lme_coordinates$long
)

if (
  anyNA(lme_coordinates) ||
    any(!is.finite(lme_coordinates$lati)) ||
    any(!is.finite(lme_coordinates$long))
) {
  stop(
    "The annual LME coordinate table contains invalid values."
  )
}

if (anyDuplicated(
  lme_coordinates[
    ,
    c(
      "long",
      "lati"
    )
  ]
)) {
  stop(
    "The annual LME coordinate table contains duplicated ",
    "locations."
  )
}

n_lme_grid <- nrow(
  lme_coordinates
)

if (n_lme_grid != 266L) {
  warning(
    "Expected 266 native LME locations, but found ",
    n_lme_grid,
    "."
  )
}

all_prediction_coordinates <- rbind(
  unname(
    loc@coords[
      ,
      c(
        "long",
        "lat"
      )
    ]
  ),
  as.matrix(
    lme_coordinates[
      ,
      c(
        "long",
        "lati"
      )
    ]
  )
)

grid_prediction_rows <- seq_len(
  n_grid
)

lme_prediction_rows <- n_grid +
  seq_len(
    n_lme_grid
  )

# ------------------------------------------------------------
# 5. Interval-censored Gaussian moment helpers
# ------------------------------------------------------------

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

Ez_discrete <- function(
    sigma,
    cuts,
    vals) {

  stopifnot(
    length(vals) ==
      length(cuts) - 1L
  )

  lower <- head(
    cuts,
    -1L
  ) /
    sigma

  upper <- tail(
    cuts,
    -1L
  ) /
    sigma

  probabilities <- pnorm(
    upper
  ) -
    pnorm(
      lower
    )

  sum(
    vals *
      probabilities
  )
}


Varz_discrete <- function(
    sigma,
    cuts,
    vals,
    mean_z = NULL) {

  if (is.null(mean_z)) {
    mean_z <- Ez_discrete(
      sigma,
      cuts,
      vals
    )
  }

  lower <- head(
    cuts,
    -1L
  ) /
    sigma

  upper <- tail(
    cuts,
    -1L
  ) /
    sigma

  probabilities <- pnorm(
    upper
  ) -
    pnorm(
      lower
    )

  sum(
    vals^2 *
      probabilities
  ) -
    mean_z^2
}


EZstar_h <- function(
    sigma,
    cuts,
    vals) {

  stopifnot(
    length(vals) ==
      length(cuts) - 1L
  )

  lower <- head(
    cuts,
    -1L
  ) /
    sigma

  upper <- tail(
    cuts,
    -1L
  ) /
    sigma

  truncated_first_moment <- sigma *
    (
      dnorm(
        lower
      ) -
        dnorm(
          upper
        )
    )

  sum(
    vals *
      truncated_first_moment
  )
}


ordinal_covariance <- function(
    rho,
    sigma) {

  if (
    !is.finite(rho) ||
      rho < -1 ||
      rho > 1
  ) {
    stop(
      "rho must lie in [-1, 1]."
    )
  }

  if (abs(rho) < 1e-14) {
    return(0)
  }

  # Z = -2 + I(W >= -1.5) + I(W >= -0.5) +
  #     I(W >= 0.5).
  # The constant does not affect covariance, so Cov(Z1, Z2)
  # is the sum of nine threshold-indicator covariances.

  standardized_thresholds <- c(
    -1.5,
    -0.5,
    0.5
  ) /
    sigma

  exceedance_probabilities <- 1 -
    pnorm(
      standardized_thresholds
    )

  correlation_matrix <- matrix(
    c(
      1,
      rho,
      rho,
      1
    ),
    nrow = 2,
    byrow = TRUE
  )

  covariance_sum <- 0

  for (
    threshold_1_index in
      seq_along(
        standardized_thresholds
      )
  ) {
    for (
      threshold_2_index in
        seq_along(
          standardized_thresholds
        )
    ) {

      threshold_1 <- standardized_thresholds[
        threshold_1_index
      ]

      threshold_2 <- standardized_thresholds[
        threshold_2_index
      ]

      lower_tail_probability <- as.numeric(
        mvtnorm::pmvnorm(
          lower = c(
            -Inf,
            -Inf
          ),
          upper = c(
            threshold_1,
            threshold_2
          ),
          mean = c(
            0,
            0
          ),
          corr = correlation_matrix,
          algorithm = mvtnorm::TVPACK()
        )
      )

      joint_exceedance_probability <-
        1 -
        pnorm(
          threshold_1
        ) -
        pnorm(
          threshold_2
        ) +
        lower_tail_probability

      covariance_sum <- covariance_sum +
        joint_exceedance_probability -
        exceedance_probabilities[
          threshold_1_index
        ] *
        exceedance_probabilities[
          threshold_2_index
        ]
    }
  }

  covariance_sum
}


build_covariance_lookup <- function(
    sigma_Y2,
    sigma_E2,
    sigma,
    number_of_points) {

  total_variance <- sigma_Y2 +
    sigma_E2

  maximum_rho <- sigma_Y2 /
    total_variance

  rho_grid <- seq(
    0,
    maximum_rho,
    length.out = number_of_points
  )

  covariance_grid <- vapply(
    rho_grid,
    ordinal_covariance,
    numeric(1),
    sigma = sigma
  )

  # The theoretical covariance is monotone in rho. cummax()
  # removes negligible numerical reversals before interpolation.
  covariance_grid <- cummax(
    covariance_grid
  )

  list(
    maximum_rho = maximum_rho,
    table = data.frame(
      rho = rho_grid,
      ordinal_covariance =
        covariance_grid
    ),
    function_object = approxfun(
      x = rho_grid,
      y = covariance_grid,
      method = "linear",
      rule = 2,
      ties = "ordered"
    )
  )
}

# ------------------------------------------------------------
# 6. Precompute fixed marginal quantities and covariance lookup
# ------------------------------------------------------------

total_variance <- sigma_Y2 +
  sigma_E2

sigma_total <- sqrt(
  total_variance
)

mean_z <- Ez_discrete(
  sigma = sigma_total,
  cuts = cuts,
  vals = vals
)

variance_z <- Varz_discrete(
  sigma = sigma_total,
  cuts = cuts,
  vals = vals,
  mean_z = mean_z
)

ezstar_h <- EZstar_h(
  sigma = sigma_total,
  cuts = cuts,
  vals = vals
)

covariance_lookup <- build_covariance_lookup(
  sigma_Y2 = sigma_Y2,
  sigma_E2 = sigma_E2,
  sigma = sigma_total,
  number_of_points = rho_lookup_size
)

readr::write_csv(
  covariance_lookup$table,
  file.path(
    output_directory,
    "ordinal_covariance_lookup.csv"
  )
)

# ------------------------------------------------------------
# 7. Fast and stable linear-system solver
# ------------------------------------------------------------

solve_symmetric_system <- function(
    covariance_matrix,
    right_hand_side) {

  jitter_values <- c(
    0,
    1e-12,
    1e-10,
    1e-8,
    1e-6
  )

  n <- nrow(
    covariance_matrix
  )

  for (jitter in jitter_values) {

    current_matrix <- covariance_matrix

    if (jitter > 0) {
      diag(
        current_matrix
      ) <- diag(
        current_matrix
      ) +
        jitter
    }

    cholesky_factor <- try(
      chol(
        current_matrix
      ),
      silent = TRUE
    )

    if (!inherits(
      cholesky_factor,
      "try-error"
    )) {

      solution <- backsolve(
        cholesky_factor,
        forwardsolve(
          t(
            cholesky_factor
          ),
          right_hand_side
        )
      )

      return(
        list(
          solution = solution,
          jitter = jitter
        )
      )
    }
  }

  stop(
    "Unable to obtain a positive-definite covariance ",
    "matrix even after adding numerical jitter. ",
    "Reciprocal condition number = ",
    signif(
      rcond(
        covariance_matrix
      ),
      6
    ),
    "."
  )
}

# ------------------------------------------------------------
# 8. Krige one event year
# ------------------------------------------------------------

predict_grid_one_year <- function(
    current_year,
    observation_data,
    prediction_coordinates,
    sigma_Y2,
    sigma_E2,
    alpha_km,
    total_variance,
    mean_z,
    variance_z,
    ezstar_h,
    covariance_lookup_function) {

  current_data <- observation_data[
    observation_data$year ==
      current_year,
  ]

  if (nrow(current_data) == 0L) {
    stop(
      "No REACHES observations were found for year ",
      current_year,
      "."
    )
  }

  observation_coordinates <- current_data@coords
  observed_indices <- as.numeric(
    current_data$level
  )

  if (
    length(observed_indices) !=
      nrow(observation_coordinates)
  ) {
    stop(
      "The observation values and coordinates have ",
      "different lengths for year ",
      current_year,
      "."
    )
  }

  # Great-circle observation-to-observation distances in km.
  observation_distance_matrix <- sp::spDists(
    observation_coordinates,
    longlat = TRUE
  )

  rho_matrix <- (
    sigma_Y2 *
      exp(
        -observation_distance_matrix /
          alpha_km
      )
  ) /
    total_variance

  SigmaZ <- matrix(
    covariance_lookup_function(
      as.vector(
        rho_matrix
      )
    ),
    nrow = nrow(
      rho_matrix
    ),
    ncol = ncol(
      rho_matrix
    )
  )

  # Distinct documentary records at an identical location retain
  # rho = sigma_Y2 / total_variance off the diagonal. Only the
  # self-covariance diagonal is Var(Z).
  diag(
    SigmaZ
  ) <- variance_z

  SigmaZ <- (
    SigmaZ +
      t(
        SigmaZ
      )
  ) /
    2

  # Vectorized observation-to-grid distances in km.
  observation_to_grid_distance <- sp::spDists(
    x = observation_coordinates,
    y = prediction_coordinates,
    longlat = TRUE
  )

  latent_cross_covariance <- sigma_Y2 *
    exp(
      -observation_to_grid_distance /
        alpha_km
    )

  CZY <- (
    latent_cross_covariance /
      total_variance
  ) *
    ezstar_h

  right_hand_side <- cbind(
    observed_indices -
      mean_z,
    CZY
  )

  solved_system <- solve_symmetric_system(
    covariance_matrix = SigmaZ,
    right_hand_side =
      right_hand_side
  )

  all_solutions <- solved_system$solution

  mean_weights <- all_solutions[
    ,
    1L
  ]

  variance_weights <- all_solutions[
    ,
    -1L,
    drop = FALSE
  ]

  prediction_mean <- drop(
    crossprod(
      CZY,
      mean_weights
    )
  )

  prediction_variance <- sigma_Y2 -
    colSums(
      CZY *
        variance_weights
    )

  prediction_variance <- pmax(
    prediction_variance,
    0
  )

  list(
    mean = prediction_mean,
    variance = prediction_variance,
    jitter = solved_system$jitter,
    number_of_observations =
      length(
        observed_indices
      )
  )
}

# ------------------------------------------------------------
# 9. Initialise or resume the full-grid calculation
# ------------------------------------------------------------

mean_matrix <- matrix(
  NA_real_,
  nrow = n_grid,
  ncol = length(year2),
  dimnames = list(
    NULL,
    as.character(
      year2
    )
  )
)

variance_matrix <- matrix(
  NA_real_,
  nrow = n_grid,
  ncol = length(year2),
  dimnames = list(
    NULL,
    as.character(
      year2
    )
  )
)

lme_mean_matrix <- matrix(
  NA_real_,
  nrow = n_lme_grid,
  ncol = length(year2),
  dimnames = list(
    as.character(
      lme_coordinates$location_id
    ),
    as.character(
      year2
    )
  )
)

lme_variance_matrix <- matrix(
  NA_real_,
  nrow = n_lme_grid,
  ncol = length(year2),
  dimnames = list(
    as.character(
      lme_coordinates$location_id
    ),
    as.character(
      year2
    )
  )
)

year_diagnostics <- data.frame(
  year = year2,
  number_of_observations =
    NA_integer_,
  numerical_jitter =
    NA_real_
)

completed_indices <- integer(0)

if (
  isTRUE(
    resume_from_checkpoint
  ) &&
    file.exists(
      checkpoint_file
    )
) {

  checkpoint <- readRDS(
    checkpoint_file
  )

  checkpoint_is_compatible <-
    identical(
      checkpoint$years,
      year2
    ) &&
    isTRUE(
      all.equal(
        checkpoint$parameters,
        parameter_vector,
        tolerance = 0
      )
    ) &&
    identical(
      dim(
        checkpoint$mean_matrix
      ),
      dim(
        mean_matrix
      )
    ) &&
    identical(
      dim(
        checkpoint$lme_mean_matrix
      ),
      dim(
        lme_mean_matrix
      )
    )

  if (checkpoint_is_compatible) {

    mean_matrix <- checkpoint$mean_matrix
    variance_matrix <-
      checkpoint$variance_matrix
    lme_mean_matrix <-
      checkpoint$lme_mean_matrix
    lme_variance_matrix <-
      checkpoint$lme_variance_matrix
    year_diagnostics <-
      checkpoint$year_diagnostics
    completed_indices <-
      checkpoint$completed_indices

    message(
      "Resuming from checkpoint with ",
      length(
        completed_indices
      ),
      " completed years."
    )

  } else {

    warning(
      "An incompatible checkpoint was found and ignored: ",
      checkpoint_file
    )
  }
}

# ------------------------------------------------------------
# 10. Generate predictions for every event year
# ------------------------------------------------------------

start_time <- Sys.time()

indices_to_run <- setdiff(
  seq_along(
    year2
  ),
  completed_indices
)

initial_completed_count <- length(
  completed_indices
)

for (
  iteration_index in
    seq_along(
      indices_to_run
    )
) {

  i <- indices_to_run[
    iteration_index
  ]

  current_year <- year2[
    i
  ]

  message(
    "Kriging year ",
    current_year,
    " (",
    initial_completed_count +
      iteration_index,
    "/",
    length(
      year2
    ),
    ")"
  )

  yearly_result <- predict_grid_one_year(
    current_year = current_year,
    observation_data = y2,
    prediction_coordinates =
      all_prediction_coordinates,
    sigma_Y2 = sigma_Y2,
    sigma_E2 = sigma_E2,
    alpha_km = alpha_km,
    total_variance =
      total_variance,
    mean_z = mean_z,
    variance_z = variance_z,
    ezstar_h = ezstar_h,
    covariance_lookup_function =
      covariance_lookup$function_object
  )

  mean_matrix[
    ,
    i
  ] <- yearly_result$mean[
    grid_prediction_rows
  ]

  variance_matrix[
    ,
    i
  ] <- yearly_result$variance[
    grid_prediction_rows
  ]

  lme_mean_matrix[
    ,
    i
  ] <- yearly_result$mean[
    lme_prediction_rows
  ]

  lme_variance_matrix[
    ,
    i
  ] <- yearly_result$variance[
    lme_prediction_rows
  ]

  year_diagnostics$number_of_observations[
    i
  ] <- yearly_result$number_of_observations

  year_diagnostics$numerical_jitter[
    i
  ] <- yearly_result$jitter

  completed_indices <- sort(
    unique(
      c(
        completed_indices,
        i
      )
    )
  )

  should_checkpoint <-
    length(
      completed_indices
    ) %% checkpoint_every == 0L ||
    length(
      completed_indices
    ) == length(
      year2
    )

  if (should_checkpoint) {
    saveRDS(
      list(
        years = year2,
        parameters = parameter_vector,
        mean_matrix = mean_matrix,
        variance_matrix =
          variance_matrix,
        lme_mean_matrix =
          lme_mean_matrix,
        lme_variance_matrix =
          lme_variance_matrix,
        year_diagnostics =
          year_diagnostics,
        completed_indices =
          completed_indices
      ),
      checkpoint_file
    )
  }
}

elapsed_minutes <- as.numeric(
  difftime(
    Sys.time(),
    start_time,
    units = "mins"
  )
)

# ------------------------------------------------------------
# 11. Assemble complete-grid outputs
# ------------------------------------------------------------

grid_coordinates <- as.data.frame(
  loc@coords
)

names(
  grid_coordinates
) <- c(
  "long",
  "lat"
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

grid_mean_file <- file.path(
  output_directory,
  "reaches_kriging_grid53x49_mean.csv.gz"
)

grid_variance_file <- file.path(
  output_directory,
  "reaches_kriging_grid53x49_variance.csv.gz"
)

readr::write_csv(
  grid_mean,
  grid_mean_file
)

readr::write_csv(
  grid_variance,
  grid_variance_file
)

# ------------------------------------------------------------
# 12. Save REACHES predictions at native LME grid locations
# ------------------------------------------------------------

lme_grid_mean <- cbind(
  lme_coordinates,
  as.data.frame(
    lme_mean_matrix,
    check.names = FALSE
  )
)

lme_grid_variance <- cbind(
  lme_coordinates,
  as.data.frame(
    lme_variance_matrix,
    check.names = FALSE
  )
)

lme_year_columns <- as.character(
  year2
)

names(
  lme_grid_mean
)[
  names(lme_grid_mean) %in%
    lme_year_columns
] <- paste0(
  "x",
  lme_year_columns
)

names(
  lme_grid_variance
)[
  names(lme_grid_variance) %in%
    lme_year_columns
] <- paste0(
  "x",
  lme_year_columns
)

lme_grid_mean_file <- file.path(
  output_directory,
  "reaches_kriging_lme_grid_mean.csv"
)

lme_grid_variance_file <- file.path(
  output_directory,
  "reaches_kriging_lme_grid_variance.csv"
)

readr::write_csv(
  lme_grid_mean,
  lme_grid_mean_file
)

readr::write_csv(
  lme_grid_variance,
  lme_grid_variance_file
)

# ------------------------------------------------------------
# 13. Extract Hong Kong, Shanghai, and Beijing
# ------------------------------------------------------------

city_locations <- data.frame(
  city = c(
    "HongKong",
    "Shanghai",
    "Beijing"
  ),
  long = c(
    114.167,
    121.433,
    116.283
  ),
  lat = c(
    22.333,
    31.167,
    39.933
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
    "At least one city location could not be matched to ",
    "the 53 x 49 prediction grid."
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

year_columns <- setdiff(
  names(
    city3_sd
  ),
  c(
    "long",
    "lat"
  )
)

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

city3_mean_file <- file.path(
  output_directory,
  "reaches_kriging_city3_mean.csv"
)

city3_sd_file <- file.path(
  output_directory,
  "reaches_kriging_city3_sd.csv"
)

readr::write_csv(
  city3_mean,
  city3_mean_file
)

readr::write_csv(
  city3_sd,
  city3_sd_file
)

# ------------------------------------------------------------
# 14. Save diagnostics and metadata
# ------------------------------------------------------------

metadata <- data.frame(
  number_of_event_years =
    length(
      year2
    ),
  number_of_grid_locations =
    n_grid,
  number_of_lme_grid_locations =
    n_lme_grid,
  process_variance =
    sigma_Y2,
  nugget_variance =
    sigma_E2,
  total_variance =
    total_variance,
  range_km =
    alpha_km,
  covariance_lookup_size =
    rho_lookup_size,
  elapsed_minutes =
    elapsed_minutes
)

metadata_file <- file.path(
  output_directory,
  "reaches_kriging_metadata.csv"
)

diagnostics_file <- file.path(
  output_directory,
  "reaches_kriging_year_diagnostics.csv"
)

readr::write_csv(
  metadata,
  metadata_file
)

readr::write_csv(
  year_diagnostics,
  diagnostics_file
)

if (file.exists(
  checkpoint_file
)) {
  file.remove(
    checkpoint_file
  )
}

message(
  "Complete-grid kriging means saved to: ",
  grid_mean_file
)

message(
  "Complete-grid prediction variances saved to: ",
  grid_variance_file
)

message(
  "LME-grid kriging means saved to: ",
  lme_grid_mean_file
)

message(
  "LME-grid prediction variances saved to: ",
  lme_grid_variance_file
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
  "Elapsed time for this run: ",
  round(
    elapsed_minutes,
    2
  ),
  " minutes."
)
