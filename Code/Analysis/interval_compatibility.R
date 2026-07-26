here::i_am("Code/Analysis/interval_compatibility.R")

# ============================================================
# Held-out interval-compatibility validation for Section 5.2
#
# For each repetition:
#   1. Keep event years with at least five distinct REACHES
#      sites.
#   2. Randomly hold out one site from each eligible year.
#   3. Re-estimate the pooled pilot variogram and apply the
#      interval-censoring calibration to the covariance
#      parameters using only the training data.
#   4. Predict the latent process at every held-out site using
#      the interval-censored best linear predictor.
#   5. Approximate the predictive distribution of
#         W = Y + epsilon
#      by N(yhat, MSPE + sigma_epsilon^2).
#   6. Check whether the 95%, 90%, 85%, and 80% predictive
#      intervals intersect the censoring interval associated
#      with the held-out ordinal category.
#
# Main outputs:
#   Output/Validation/interval_compatibility_summary.csv
#   Output/Validation/interval_compatibility_by_split.csv
#   Output/Validation/interval_compatibility_diagnostics.csv
#
# Detailed outputs:
#   Output/Intermediate/interval_compatibility/
# ============================================================

library(here)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(sp)
library(spacetime)
library(zoo)
library(gstat)
library(mgcv)
library(mvtnorm)

# ------------------------------------------------------------
# 1. Analysis settings
# ------------------------------------------------------------

reaches_file <- here::here(
  "Data",
  "temperature index value.v1.xlsx"
)

validation_output_dir <- here::here(
  "Output",
  "Validation"
)

intermediate_output_dir <- here::here(
  "Output",
  "Intermediate",
  "interval_compatibility"
)

dir.create(
  validation_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  intermediate_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

analysis_years <- 1368:1911

number_of_repetitions <- 30L
minimum_sites_per_year <- 5L

nominal_levels <- c(
  0.95,
  0.90,
  0.85,
  0.80
)

master_seed <- 10L

# Monte Carlo grid sizes used in the original calibration.
calibration_n1 <- 500L
calibration_n2 <- 200L
calibration_n3 <- 40L

# Numerical settings.
solve_tolerance <- 1e-12

ordinal_values <- c(
  -2,
  -1,
  0,
  1
)

ordinal_cuts <- c(
  -Inf,
  -1.5,
  -0.5,
  0.5,
  Inf
)

if (!file.exists(reaches_file)) {
  stop(
    "The REACHES input file was not found: ",
    reaches_file
  )
}

# ------------------------------------------------------------
# 2. Read and clean REACHES data
# ------------------------------------------------------------

temperature <- readxl::read_excel(
  reaches_file,
  col_types = c(
    "skip",
    "skip",
    "numeric",
    "numeric",
    "skip",
    "skip",
    "skip",
    "skip",
    "skip",
    "numeric",
    "numeric",
    "skip",
    "skip"
  )
)

names(temperature) <- c(
  "level",
  "year",
  "long",
  "lat"
)

temperature <- temperature %>%
  transmute(
    level = as.numeric(level),
    year = as.integer(year),
    long = as.numeric(long),
    lat = as.numeric(lat),
    source_order = row_number()
  ) %>%
  filter(
    year %in% analysis_years,
    level %in% ordinal_values,
    is.finite(long),
    is.finite(lat)
  )

# The historical code stated that repeated records at the same
# site and year should be represented by their modal category.
# This implementation performs that aggregation explicitly.
# When two or more categories tie for the mode, the category
# appearing first in the source file is retained.
mode_first <- function(
    values,
    source_order) {

  counts <- table(values)
  maximum_count <- max(counts)

  modal_values <- as.numeric(
    names(
      counts[
        counts == maximum_count
      ]
    )
  )

  if (length(modal_values) == 1L) {
    return(modal_values)
  }

  tied_rows <- which(
    values %in% modal_values
  )

  values[tied_rows[which.min(source_order[tied_rows])]]
}

duplicate_diagnostics <- temperature %>%
  count(
    year,
    long,
    lat,
    name = "number_of_records"
  ) %>%
  summarise(
    number_of_site_years = n(),
    duplicated_site_years = sum(
      number_of_records > 1L
    ),
    maximum_records_at_one_site_year = max(
      number_of_records
    )
  )

temperature_site_year <- temperature %>%
  group_by(
    year,
    long,
    lat
  ) %>%
  summarise(
    number_of_source_records = n(),
    number_of_distinct_levels = n_distinct(
      level
    ),
    level = mode_first(
      level,
      source_order
    ),
    .groups = "drop"
  ) %>%
  arrange(
    year,
    long,
    lat
  )

readr::write_csv(
  duplicate_diagnostics,
  file.path(
    intermediate_output_dir,
    "duplicate_site_year_diagnostics.csv"
  )
)

readr::write_csv(
  temperature_site_year,
  file.path(
    intermediate_output_dir,
    "reaches_site_year_mode.csv"
  )
)

# ------------------------------------------------------------
# 3. Hold-out split
# ------------------------------------------------------------

split_validation_data <- function(
    data,
    minimum_sites = 5L) {

  eligible_years <- data %>%
    count(
      year,
      name = "number_of_sites"
    ) %>%
    filter(
      number_of_sites >= minimum_sites
    ) %>%
    pull(
      year
    )

  if (length(eligible_years) == 0L) {
    stop(
      "No event years contain at least ",
      minimum_sites,
      " distinct sites."
    )
  }

  test_data <- data %>%
    filter(
      year %in% eligible_years
    ) %>%
    group_by(
      year
    ) %>%
    slice_sample(
      n = 1L
    ) %>%
    ungroup() %>%
    arrange(
      year
    )

  training_data <- anti_join(
    data,
    test_data %>%
      select(
        year,
        long,
        lat
      ),
    by = c(
      "year",
      "long",
      "lat"
    )
  )

  list(
    training_data = training_data,
    test_data = test_data,
    validation_years = test_data$year
  )
}

# ------------------------------------------------------------
# 4. Pilot variogram estimation
# ------------------------------------------------------------

estimate_pilot_variogram <- function(
    training_data) {

  y2 <- training_data %>%
    filter(
      year > 100L,
      year <= 2000L
    )

  coordinate_data <- unique(
    round(
      data.frame(
        long = y2$long,
        lat = y2$lat
      ),
      4
    )
  )

  sp::coordinates(
    coordinate_data
  ) <- ~ long + lat

  sp::proj4string(
    coordinate_data
  ) <- sp::CRS(
    "+proj=longlat +datum=WGS84"
  )

  years <- sort(
    unique(
      y2$year
    )
  )

  year_month_index <- zoo::as.yearmon(
    years
  )

  y2_rounded <- as.data.frame(
    round(
      as.matrix(
        y2 %>%
          select(
            level,
            year,
            long,
            lat
          )
      ),
      4
    )
  )

  index_matrix <- matrix(
    NA_integer_,
    nrow = nrow(y2_rounded),
    ncol = 2L
  )

  index_matrix[, 2] <- match(
    y2_rounded$year,
    years
  )

  index_matrix[, 1] <- match(
    y2_rounded$long * 10^7 +
      y2_rounded$lat,
    coordinate_data@coords[, 1] * 10^7 +
      coordinate_data@coords[, 2]
  )

  if (anyNA(index_matrix)) {
    stop(
      "Failed to map at least one training observation to the ",
      "spatio-temporal index used for pilot variogram fitting."
    )
  }

  sts_data <- spacetime::STSDF(
    sp = coordinate_data,
    time = year_month_index,
    data = y2,
    index = index_matrix
  )

  variogram_st <- gstat::variogramST(
    level ~ 1,
    data = sts_data,
    tlags = 0,
    width = 10,
    na.omit = TRUE
  )

  variogram_st <- subset(
    variogram_st,
    dist > 1
  )

  y2_spatial <- y2
  sp::coordinates(
    y2_spatial
  ) <- ~ long + lat

  sp::proj4string(
    y2_spatial
  ) <- sp::CRS(
    "+proj=longlat +datum=WGS84"
  )

  variogram_spatial <- gstat::variogram(
    level ~ 1,
    data = y2_spatial,
    width = 10
  )

  variogram_indices <- match(
    as.numeric(
      rownames(
        variogram_st
      )
    ) -
      1,
    rownames(
      variogram_spatial
    )
  )

  keep_variogram_rows <- rep(
    FALSE,
    nrow(
      variogram_spatial
    )
  )

  keep_variogram_rows[
    variogram_indices[
      !is.na(
        variogram_indices
      )
    ]
  ] <- TRUE

  variogram_spatial <- variogram_spatial[
    keep_variogram_rows,
  ]

  if (
    sum(
      variogram_st$dist == 0
    ) >
      0
  ) {
    variogram_spatial[, 1:3] <- variogram_st[
      variogram_indices[
        -1
      ],
      1:3
    ]
  } else {
    variogram_spatial[, 1:3] <- variogram_st[
      variogram_indices,
      1:3
    ]
  }

  variogram_fit <- gstat::fit.variogram(
    variogram_spatial,
    gstat::vgm(
      model = "Exp",
      nugget = NA
    ),
    fit.kappa = TRUE,
    fit.method = 2
  )

  if (
    nrow(variogram_fit) < 2L ||
      any(!is.finite(variogram_fit$psill)) ||
      !is.finite(variogram_fit$range[2])
  ) {
    stop(
      "Pilot exponential variogram fitting failed."
    )
  }

  list(
    psill = variogram_fit$psill,
    range = variogram_fit$range[2],
    training_spatial = y2_spatial,
    variogram_fit = variogram_fit
  )
}

# ------------------------------------------------------------
# 5. Monte Carlo calibration of covariance parameters
# ------------------------------------------------------------

find_calibrated_x <- function(
    x,
    y,
    target_y) {

  calibration_data <- data.frame(
    x = x,
    y = y
  )

  calibration_fit <- mgcv::gam(
    y ~ s(x),
    data = calibration_data
  )

  prediction_grid <- seq(
    min(calibration_data$x),
    max(calibration_data$x),
    length.out = 1000L
  )

  fitted_values <- predict(
    calibration_fit,
    newdata = data.frame(
      x = prediction_grid
    )
  )

  prediction_grid[
    which.min(
      abs(
        fitted_values -
          target_y
      )
    )
  ]
}


round_and_clip_ordinal <- function(
    x) {

  rounded <- round(
    x,
    digits = 0
  )

  pmin(
    pmax(
      rounded,
      -2
    ),
    1
  )
}


calibrate_covariance_parameters <- function(
    pilot_nugget,
    pilot_psill,
    pilot_range,
    n1 = 500L,
    n2 = 200L,
    n3 = 40L) {

  # f1: latent total variance -> variance of rounded categories.
  latent_total_variance <- numeric(
    n1
  )

  rounded_total_variance <- numeric(
    n1
  )

  for (i in seq_len(n1)) {

    target_variance <- 2 * i /
      n1

    latent_sample <- rnorm(
      10000L,
      mean = 0,
      sd = sqrt(
        target_variance
      )
    )

    rounded_sample <- round_and_clip_ordinal(
      latent_sample
    )

    latent_total_variance[i] <- var(
      latent_sample
    )

    rounded_total_variance[i] <- var(
      rounded_sample
    )
  }

  total_variance <- find_calibrated_x(
    x = latent_total_variance,
    y = rounded_total_variance,
    target_y = sum(
      pilot_psill
    )
  )

  # f2: measurement-error variance -> rounded nugget variance.
  epsilon_variance <- numeric(
    n2
  )

  rounded_epsilon_variance <- numeric(
    n2
  )

  for (i in seq_len(n2)) {

    latent_process_variance <- total_variance -
      i *
      total_variance /
      n2

    latent_process <- rnorm(
      10000L,
      mean = 0,
      sd = sqrt(
        latent_process_variance
      )
    )

    epsilon_pair <- mvtnorm::rmvnorm(
      10000L,
      sigma = diag(
        total_variance -
          latent_process_variance,
        2L
      )
    )

    z1 <- round_and_clip_ordinal(
      latent_process +
        epsilon_pair[, 1]
    )

    z2 <- round_and_clip_ordinal(
      latent_process +
        epsilon_pair[, 2]
    )

    epsilon_variance[i] <- var(
      epsilon_pair[, 1] -
        epsilon_pair[, 2]
    ) /
      2

    rounded_epsilon_variance[i] <- var(
      z1 -
        z2
    ) /
      2
  }

  sigma_epsilon2 <- find_calibrated_x(
    x = epsilon_variance,
    y = rounded_epsilon_variance,
    target_y = pilot_nugget
  )

  sigma_Y2 <- total_variance -
    sigma_epsilon2

  if (
    !is.finite(sigma_Y2) ||
      sigma_Y2 <= 0 ||
      !is.finite(sigma_epsilon2) ||
      sigma_epsilon2 < 0
  ) {
    stop(
      "The calibrated latent-process or measurement-error ",
      "variance is invalid."
    )
  }

  # f3: latent range -> range estimated after ordinal rounding.
  candidate_latent_ranges <- 100 +
    10 *
    0:n3

  candidate_rounded_ranges <- numeric(
    length(
      candidate_latent_ranges
    )
  )

  distance_grid <- 5 *
    0:240

  rounded_covariance <- matrix(
    0,
    nrow = length(
      distance_grid
    ),
    ncol = 151L
  )

  for (range_index in seq_len(151L)) {

    rounded_range <- 95 +
      5 *
      range_index

    for (distance_index in seq_along(
      distance_grid
    )) {

      distance_value <- distance_grid[
        distance_index
      ]

      covariance_matrix <- matrix(
        c(
          sigma_Y2 +
            sigma_epsilon2,
          sigma_Y2 *
            exp(
              -distance_value /
                rounded_range
            ),
          sigma_Y2 *
            exp(
              -distance_value /
                rounded_range
            ),
          sigma_Y2 +
            sigma_epsilon2
        ),
        nrow = 2L,
        byrow = TRUE
      )

      latent_pair <- mvtnorm::rmvnorm(
        1000L,
        sigma = covariance_matrix
      )

      rounded_pair <- round_and_clip_ordinal(
        latent_pair
      )

      rounded_covariance[
        distance_index,
        range_index
      ] <- cov(
        rounded_pair[, 1],
        rounded_pair[, 2]
      )
    }
  }

  search_start <- 1L

  for (i in seq_along(
    candidate_latent_ranges
  )) {

    latent_range <- candidate_latent_ranges[i]

    objective_values <- rep(
      Inf,
      151L
    )

    for (range_index in search_start:151L) {

      target_covariance <- sigma_Y2 *
        exp(
          -distance_grid /
            latent_range
        )

      objective_values[
        range_index
      ] <- sum(
        (
          rounded_covariance[
            ,
            range_index
          ] -
            target_covariance
        )^2
      )
    }

    best_index <- which.min(
      objective_values
    )

    candidate_rounded_ranges[i] <- 95 +
      5 *
      best_index

    search_start <- min(
      best_index,
      150L
    )
  }

  calibrated_range <- find_calibrated_x(
    x = candidate_latent_ranges,
    y = candidate_rounded_ranges,
    target_y = pilot_range
  )

  list(
    sigma_epsilon2 = sigma_epsilon2,
    sigma_Y2 = sigma_Y2,
    total_variance = total_variance,
    range = calibrated_range
  )
}

# ------------------------------------------------------------
# 6. Interval-censored best linear predictor
# ------------------------------------------------------------

Ez_discrete <- function(
    sigma,
    cuts = ordinal_cuts,
    values = ordinal_values) {

  lower_standardized <- head(
    cuts,
    -1
  ) /
    sigma

  upper_standardized <- tail(
    cuts,
    -1
  ) /
    sigma

  probabilities <- pnorm(
    upper_standardized
  ) -
    pnorm(
      lower_standardized
    )

  sum(
    values *
      probabilities
  )
}


EZstar_h <- function(
    sigma,
    cuts = ordinal_cuts,
    values = ordinal_values) {

  lower_standardized <- head(
    cuts,
    -1
  ) /
    sigma

  upper_standardized <- tail(
    cuts,
    -1
  ) /
    sigma

  truncated_first_moment <- sigma *
    (
      dnorm(
        lower_standardized
      ) -
        dnorm(
          upper_standardized
        )
    )

  sum(
    values *
      truncated_first_moment
  )
}


cov_Z_pair <- function(
    rho,
    sigma,
    cuts = ordinal_cuts,
    values = ordinal_values,
    expected_Z = NULL) {

  if (is.null(expected_Z)) {
    expected_Z <- Ez_discrete(
      sigma,
      cuts,
      values
    )
  }

  lower_standardized <- head(
    cuts,
    -1
  ) /
    sigma

  upper_standardized <- tail(
    cuts,
    -1
  ) /
    sigma

  correlation_matrix <- matrix(
    c(
      1,
      rho,
      rho,
      1
    ),
    nrow = 2L
  )

  expected_product <- 0

  for (k in seq_along(
    values
  )) {
    for (l in seq_along(
      values
    )) {

      rectangle_probability <- as.numeric(
        mvtnorm::pmvnorm(
          lower = c(
            lower_standardized[k],
            lower_standardized[l]
          ),
          upper = c(
            upper_standardized[k],
            upper_standardized[l]
          ),
          mean = c(
            0,
            0
          ),
          sigma = correlation_matrix
        )
      )

      expected_product <- expected_product +
        values[k] *
        values[l] *
        rectangle_probability
    }
  }

  expected_product -
    expected_Z^2
}


SigmaZ_matrix <- function(
    observation_coordinates,
    sigma_Y2,
    sigma_epsilon2,
    range_parameter,
    cuts = ordinal_cuts,
    values = ordinal_values) {

  number_of_observations <- nrow(
    observation_coordinates
  )

  total_variance <- sigma_Y2 +
    sigma_epsilon2

  total_sd <- sqrt(
    total_variance
  )

  expected_Z <- Ez_discrete(
    total_sd,
    cuts,
    values
  )

  distance_matrix <- as.matrix(
    dist(
      observation_coordinates
    )
  )

  correlation_matrix <- (
    sigma_Y2 *
      exp(
        -distance_matrix /
          range_parameter
      )
  ) /
    total_variance

  diag(
    correlation_matrix
  ) <- 1

  SigmaZ <- matrix(
    NA_real_,
    nrow = number_of_observations,
    ncol = number_of_observations
  )

  for (i in seq_len(
    number_of_observations
  )) {
    for (j in i:number_of_observations) {

      covariance_value <- cov_Z_pair(
        rho = correlation_matrix[i, j],
        sigma = total_sd,
        cuts = cuts,
        values = values,
        expected_Z = expected_Z
      )

      SigmaZ[i, j] <- covariance_value
      SigmaZ[j, i] <- covariance_value
    }
  }

  SigmaZ
}


cZY_vector <- function(
    observation_coordinates,
    prediction_coordinate,
    sigma_Y2,
    sigma_epsilon2,
    range_parameter,
    cuts = ordinal_cuts,
    values = ordinal_values) {

  total_variance <- sigma_Y2 +
    sigma_epsilon2

  total_sd <- sqrt(
    total_variance
  )

  EZh <- EZstar_h(
    total_sd,
    cuts,
    values
  )

  distances <- sqrt(
    rowSums(
      (
        observation_coordinates -
          matrix(
            prediction_coordinate,
            nrow = nrow(
              observation_coordinates
            ),
            ncol = ncol(
              observation_coordinates
            ),
            byrow = TRUE
          )
      )^2
    )
  )

  latent_covariance <- sigma_Y2 *
    exp(
      -distances /
        range_parameter
    )

  (
    latent_covariance /
      total_variance
  ) *
    EZh
}


predict_held_out_sites <- function(
    training_spatial,
    test_data,
    validation_years,
    sigma_Y2,
    sigma_epsilon2,
    calibrated_range,
    solve_tolerance = 1e-2) {

  prediction_results <- vector(
    "list",
    length(
      validation_years
    )
  )

  range_parameter <- calibrated_range /
    100

  for (year_index in seq_along(
    validation_years
  )) {

    current_year <- validation_years[
      year_index
    ]

    training_year <- training_spatial[
      training_spatial$year ==
        current_year,
    ]

    held_out_year <- test_data %>%
      filter(
        year == current_year
      )

    if (
      nrow(training_year) <
        1L ||
        nrow(held_out_year) !=
        1L
    ) {
      stop(
        "Unexpected training/test size for year ",
        current_year,
        "."
      )
    }

    observation_coordinates <- training_year@coords

    prediction_coordinate <- as.numeric(
      held_out_year[
        1,
        c(
          "long",
          "lat"
        )
      ]
    )

    observed_categories <- as.numeric(
      training_year$level
    )

    SigmaZ <- SigmaZ_matrix(
      observation_coordinates = observation_coordinates,
      sigma_Y2 = sigma_Y2,
      sigma_epsilon2 = sigma_epsilon2,
      range_parameter = range_parameter
    )

    SigmaZ <- (
      SigmaZ +
        t(
          SigmaZ
        )
    ) /
      2

    cross_covariance <- cZY_vector(
      observation_coordinates = observation_coordinates,
      prediction_coordinate = prediction_coordinate,
      sigma_Y2 = sigma_Y2,
      sigma_epsilon2 = sigma_epsilon2,
      range_parameter = range_parameter
    )

    expected_Z <- Ez_discrete(
      sqrt(
        sigma_Y2 +
          sigma_epsilon2
      )
    )

    right_hand_sides <- cbind(
      observed_categories -
        expected_Z,
      cross_covariance
    )

    solutions <- tryCatch(
      solve(
        SigmaZ,
        right_hand_sides,
        tol = solve_tolerance
      ),
      error = function(error_condition) {
        stop(
          "The held-out kriging system failed for year ",
          current_year,
          ". Reciprocal condition number = ",
          signif(
            rcond(
              SigmaZ
            ),
            6
          ),
          ". Original error: ",
          conditionMessage(
            error_condition
          )
        )
      }
    )

    prediction_mean <- drop(
      crossprod(
        cross_covariance,
        solutions[, 1]
      )
    )

    prediction_mspe <- sigma_Y2 -
      drop(
        crossprod(
          cross_covariance,
          solutions[, 2]
        )
      )

    prediction_results[[year_index]] <- data.frame(
      year = current_year,
      long = held_out_year$long,
      lat = held_out_year$lat,
      Z_obs = held_out_year$level,
      yhat = prediction_mean,
      mspe = max(
        prediction_mspe,
        0
      )
    )
  }

  bind_rows(
    prediction_results
  )
}

# ------------------------------------------------------------
# 7. Validation metrics
# ------------------------------------------------------------

compute_validation_metrics <- function(
    test_predictions,
    sigma_epsilon2,
    nominal_levels,
    cuts = ordinal_cuts,
    values = ordinal_values) {

  predictive_sd <- sqrt(
    test_predictions$mspe +
      sigma_epsilon2
  )

  category_probabilities <- vapply(
    seq_along(
      values
    ),
    function(category_index) {

      pnorm(
        cuts[category_index + 1L],
        mean = test_predictions$yhat,
        sd = predictive_sd
      ) -
        pnorm(
          cuts[category_index],
          mean = test_predictions$yhat,
          sd = predictive_sd
        )
    },
    numeric(
      nrow(
        test_predictions
      )
    )
  )

  colnames(
    category_probabilities
  ) <- as.character(
    values
  )

  cumulative_probabilities <- t(
    apply(
      category_probabilities,
      1,
      cumsum
    )
  )

  observed_cumulative_indicators <- outer(
    test_predictions$Z_obs,
    values,
    "<="
  ) *
    1

  rps_by_observation <- rowSums(
    (
      cumulative_probabilities[, seq_len(length(values) - 1L), drop = FALSE] -
        observed_cumulative_indicators[, seq_len(length(values) - 1L), drop = FALSE]
    )^2
  )

  predicted_category <- values[max.col( category_probabilities, ties.method = "first")]

  category_index <- match(
    test_predictions$Z_obs,
    values
  )

  observed_lower <- cuts[category_index]

  observed_upper <- cuts[category_index +1L]

  compatibility_by_level <- lapply(
    nominal_levels,
    function(nominal_level) {

      z_value <- qnorm(
        (
          1 +
            nominal_level
        ) /
          2
      )

      predictive_lower <- test_predictions$yhat -
        z_value *
        predictive_sd

      predictive_upper <- test_predictions$yhat +
        z_value *
        predictive_sd

      compatible <- (
        predictive_upper >=
          observed_lower
      ) &
        (
          predictive_lower <
            observed_upper
        )

      data.frame(
        nominal_level = nominal_level,
        compatibility = mean(
          compatible,
          na.rm = TRUE
        )
      )
    }
  )

  list(
    compatibility = bind_rows(
      compatibility_by_level
    ),
    rps = mean(
      rps_by_observation,
      na.rm = TRUE
    ),
    accuracy = mean(
      predicted_category ==
        test_predictions$Z_obs,
      na.rm = TRUE
    ),
    rps_by_observation = rps_by_observation,
    predicted_category = predicted_category,
    category_probabilities = category_probabilities
  )
}

# ------------------------------------------------------------
# 8. Run the repeated validation
# ------------------------------------------------------------

set.seed(
  master_seed
)

start_time <- Sys.time()

split_results <- vector(
  "list",
  number_of_repetitions
)

diagnostic_results <- vector(
  "list",
  number_of_repetitions
)

observation_results <- vector(
  "list",
  number_of_repetitions
)

for (repetition in seq_len(
  number_of_repetitions
)) {

  message(
    "\nValidation repetition ",
    repetition,
    "/",
    number_of_repetitions
  )

  split_data <- split_validation_data(
    data = temperature_site_year,
    minimum_sites = minimum_sites_per_year
  )

  pilot_fit <- estimate_pilot_variogram(
    training_data = split_data$training_data
  )

  calibrated_parameters <- calibrate_covariance_parameters(
    pilot_nugget = pilot_fit$psill[1],
    pilot_psill = pilot_fit$psill,
    pilot_range = pilot_fit$range,
    n1 = calibration_n1,
    n2 = calibration_n2,
    n3 = calibration_n3
  )

  held_out_predictions <- predict_held_out_sites(
    training_spatial = pilot_fit$training_spatial,
    test_data = split_data$test_data,
    validation_years = split_data$validation_years,
    sigma_Y2 = calibrated_parameters$sigma_Y2,
    sigma_epsilon2 = calibrated_parameters$sigma_epsilon2,
    calibrated_range = calibrated_parameters$range,
    solve_tolerance = solve_tolerance
  )

  validation_metrics <- compute_validation_metrics(
    test_predictions = held_out_predictions,
    sigma_epsilon2 = calibrated_parameters$sigma_epsilon2,
    nominal_levels = nominal_levels
  )

  split_results[[repetition]] <- validation_metrics$compatibility %>%
    mutate(
      repetition = repetition,
      number_of_held_out_sites = nrow(
        held_out_predictions
      )
    ) %>%
    select(
      repetition,
      nominal_level,
      compatibility,
      number_of_held_out_sites
    )

  diagnostic_results[[repetition]] <- data.frame(
    repetition = repetition,
    number_of_held_out_sites = nrow(
      held_out_predictions
    ),
    rps = validation_metrics$rps,
    accuracy = validation_metrics$accuracy,
    pilot_nugget = pilot_fit$psill[1],
    pilot_partial_sill = pilot_fit$psill[2],
    pilot_range = pilot_fit$range,
    calibrated_sigma_epsilon2 =
      calibrated_parameters$sigma_epsilon2,
    calibrated_sigma_Y2 =
      calibrated_parameters$sigma_Y2,
    calibrated_total_variance =
      calibrated_parameters$total_variance,
    calibrated_range =
      calibrated_parameters$range
  )

  probability_data <- as.data.frame(
    validation_metrics$category_probabilities
  )

  names(
    probability_data
  ) <- paste0(
    "probability_",
    make.names(
      names(
        probability_data
      )
    )
  )

  observation_results[[repetition]] <- bind_cols(
    data.frame(
      repetition = repetition
    ),
    held_out_predictions,
    data.frame(
      predictive_sd = sqrt(
        held_out_predictions$mspe +
          calibrated_parameters$sigma_epsilon2
      ),
      rps = validation_metrics$rps_by_observation,
      predicted_category =
        validation_metrics$predicted_category
    ),
    probability_data
  )

  # Save a checkpoint after every repetition.
  checkpoint <- list(
    completed_repetitions = repetition,
    split_results = bind_rows(
      split_results[seq_len(repetition)]
    ),
    diagnostic_results = bind_rows(
      diagnostic_results[seq_len(repetition)]
    )
  )

  saveRDS(
    checkpoint,
    file.path(
      intermediate_output_dir,
      "interval_compatibility_checkpoint.rds"
    )
  )
}

compatibility_by_split <- bind_rows(
  split_results
)

diagnostics_by_split <- bind_rows(
  diagnostic_results
)

held_out_observations <- bind_rows(
  observation_results
)

# ------------------------------------------------------------
# 9. Summaries and outputs
# ------------------------------------------------------------

compatibility_summary <- compatibility_by_split %>%
  group_by(
    nominal_level
  ) %>%
  summarise(
    repetitions = n(),
    mean_compatibility = mean(
      compatibility,
      na.rm = TRUE
    ),
    sd_compatibility = sd(
      compatibility,
      na.rm = TRUE
    ),
    minimum_compatibility = min(
      compatibility,
      na.rm = TRUE
    ),
    maximum_compatibility = max(
      compatibility,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(
      nominal_level
    )
  )

compatibility_summary_formatted <- compatibility_summary %>%
  transmute(
    `Nominal level` = paste0(
      round(
        100 *
          nominal_level
      ),
      "%"
    ),
    `Mean compatibility` = sprintf(
      "%.3f",
      mean_compatibility
    ),
    `SD` = sprintf(
      "%.3f",
      sd_compatibility
    )
  )

diagnostic_summary <- data.frame(
  Metric = c(
    "Ranked probability score",
    "Classification accuracy"
  ),
  Mean = c(
    mean(
      diagnostics_by_split$rps,
      na.rm = TRUE
    ),
    mean(
      diagnostics_by_split$accuracy,
      na.rm = TRUE
    )
  ),
  SD = c(
    sd(
      diagnostics_by_split$rps,
      na.rm = TRUE
    ),
    sd(
      diagnostics_by_split$accuracy,
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  compatibility_summary,
  file.path(
    validation_output_dir,
    "interval_compatibility_summary.csv"
  )
)

readr::write_csv(
  compatibility_summary_formatted,
  file.path(
    validation_output_dir,
    "interval_compatibility_summary_formatted.csv"
  )
)

readr::write_csv(
  compatibility_by_split,
  file.path(
    validation_output_dir,
    "interval_compatibility_by_split.csv"
  )
)

readr::write_csv(
  diagnostics_by_split,
  file.path(
    validation_output_dir,
    "interval_compatibility_diagnostics.csv"
  )
)

readr::write_csv(
  diagnostic_summary,
  file.path(
    validation_output_dir,
    "interval_compatibility_diagnostic_summary.csv"
  )
)

readr::write_csv(
  held_out_observations,
  gzfile(
    file.path(
      intermediate_output_dir,
      "interval_compatibility_held_out_predictions.csv.gz"
    )
  )
)

saveRDS(
  list(
    settings = list(
      analysis_years = analysis_years,
      number_of_repetitions =
        number_of_repetitions,
      minimum_sites_per_year =
        minimum_sites_per_year,
      nominal_levels = nominal_levels,
      master_seed = master_seed,
      calibration_n1 = calibration_n1,
      calibration_n2 = calibration_n2,
      calibration_n3 = calibration_n3,
      solve_tolerance = solve_tolerance
    ),
    compatibility_summary =
      compatibility_summary,
    compatibility_by_split =
      compatibility_by_split,
    diagnostics_by_split =
      diagnostics_by_split,
    held_out_observations =
      held_out_observations
  ),
  file.path(
    intermediate_output_dir,
    "interval_compatibility_results.rds"
  )
)

elapsed_minutes <- as.numeric(
  difftime(
    Sys.time(),
    start_time,
    units = "mins"
  )
)

capture.output(
  {
    cat(
      "Held-out interval-compatibility validation\n"
    )

    cat(
      "==========================================\n\n"
    )

    cat(
      "Master seed:",
      master_seed,
      "\n"
    )

    cat(
      "Repetitions:",
      number_of_repetitions,
      "\n"
    )

    cat(
      "Minimum distinct sites per event year:",
      minimum_sites_per_year,
      "\n"
    )

    cat(
      "Elapsed time:",
      round(
        elapsed_minutes,
        2
      ),
      "minutes\n\n"
    )

    print(
      compatibility_summary_formatted
    )

    cat(
      "\nAdditional ordinal diagnostics\n"
    )

    print(
      diagnostic_summary
    )

    cat(
      "\nDuplicate site-year diagnostics\n"
    )

    print(
      duplicate_diagnostics
    )

    cat(
      "\nSession information\n"
    )

    print(
      sessionInfo()
    )
  },
  file = file.path(
    intermediate_output_dir,
    "interval_compatibility_report.txt"
  )
)

message(
  "\nInterval-compatibility validation completed."
)

message(
  "Summary saved to: ",
  file.path(
    validation_output_dir,
    "interval_compatibility_summary_formatted.csv"
  )
)

message(
  "Elapsed time: ",
  round(
    elapsed_minutes,
    2
  ),
  " minutes."
)
