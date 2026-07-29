here::i_am("Code/Analysis/interval_compatibility.R")

# ============================================================
# Held-out interval-compatibility validation for Section 5.2
#
# For each repetition:
#   1. Identify event years containing at least five distinct
#      REACHES sites.
#   2. Randomly hold out one complete site-year from each
#      eligible year.
#   3. Re-estimate the pooled pilot variogram and recalibrate
#      the latent covariance parameters using only training
#      site-years.
#   4. Predict the latent process at every held-out site using
#      the same interval-censored best linear predictor used in
#      the main kriging workflow.
#   5. Approximate
#
#        W(s_test) | training
#
#      by
#
#        N(yhat, MSPE + sigma_epsilon^2).
#
#   6. Evaluate whether central 95%, 90%, 85%, and 80%
#      predictive intervals intersect the held-out censoring
#      interval.
#
# Important consistency choices:
#   - The covariance calibration uses the corrected analytic
#     f1, f2, and Equation-(5) f3 implementation.
#   - Equation (5) targets the pilot partial sill, not the
#     calibrated process variance.
#   - Spatial distances are great-circle distances in km.
#   - All source records at the selected held-out site-year are
#     removed to prevent information leakage.
#   - Pilot calibration uses one modal category per site-year,
#     matching the covariance-calibration dataset.
#   - Prediction retains all remaining source records, matching
#     the main interval-censored kriging implementation.
#   - Only the true covariance-matrix diagonal is set to Var(Z);
#     repeated records at the same location retain the proper
#     off-diagonal covariance.
#
# Main outputs:
#   Output/Validation/
#     interval_compatibility_summary.csv
#     interval_compatibility_summary_formatted.csv
#     interval_compatibility_by_split.csv
#     interval_compatibility_diagnostics.csv
#     interval_compatibility_diagnostic_summary.csv
#
# Detailed outputs:
#   Output/Intermediate/interval_compatibility/
#     interval_compatibility_held_out_predictions.csv.gz
#     interval_compatibility_split_sites.csv.gz
#     interval_compatibility_results.rds
#     interval_compatibility_report.txt
#     interval_compatibility_checkpoint.rds
# ============================================================

library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(sp)
library(spacetime)
library(zoo)
library(gstat)
library(mvtnorm)

# ------------------------------------------------------------
# 1. Paths and analysis settings
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

checkpoint_file <- file.path(
  intermediate_output_dir,
  "interval_compatibility_checkpoint.rds"
)

analysis_years <- 1368:1911

number_of_repetitions <- 30L
minimum_distinct_sites_per_year <- 5L

nominal_levels <- c(
  0.95,
  0.90,
  0.85,
  0.80
)

master_seed <- 10L

# Calibration-grid settings. The default analytic engine is
# deterministic; no Monte Carlo calibration is required.
number_of_f1_grid_points <- 500L
number_of_f2_grid_points <- 200L

distance_grid_km <- seq(
  0,
  1200,
  by = 5
)

alpha_grid_initial_km <- seq(
  50,
  1000,
  by = 10
)

latent_range_bounds_km <- c(
  20,
  2500
)

rho_lookup_size <- 1501L

# Numerical linear-algebra settings.
solver_jitter_values <- c(
  0,
  1e-12,
  1e-10,
  1e-8,
  1e-6
)

resume_from_checkpoint <- TRUE
remove_checkpoint_on_success <- TRUE

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

if (!file.exists(
  reaches_file
)) {
  stop(
    "The REACHES input file was not found: ",
    reaches_file
  )
}

# ------------------------------------------------------------
# 2. Read and clean REACHES data
# ------------------------------------------------------------

temperature_raw <- readxl::read_excel(
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

names(
  temperature_raw
) <- c(
  "level",
  "year",
  "long",
  "lat"
)

temperature_raw <- temperature_raw %>%
  transmute(
    source_record_id = row_number(),
    source_order = row_number(),
    level = as.numeric(
      level
    ),
    year = as.integer(
      year
    ),
    long = as.numeric(
      long
    ),
    lat = as.numeric(
      lat
    )
  ) %>%
  filter(
    year %in%
      analysis_years,
    level %in%
      ordinal_values,
    is.finite(
      long
    ),
    is.finite(
      lat
    )
  ) %>%
  mutate(
    site_key = paste(
      sprintf(
        "%.8f",
        long
      ),
      sprintf(
        "%.8f",
        lat
      ),
      sep = "_"
    ),
    site_year_key = paste(
      year,
      site_key,
      sep = "_"
    )
  ) %>%
  arrange(
    year,
    long,
    lat,
    source_order
  )

if (nrow(
  temperature_raw
) == 0L) {
  stop(
    "No valid REACHES temperature observations were found."
  )
}

mode_first <- function(
    values,
    source_order) {

  counts <- table(
    values
  )

  maximum_count <- max(
    counts
  )

  modal_values <- as.numeric(
    names(
      counts[
        counts ==
          maximum_count
      ]
    )
  )

  if (length(
    modal_values
  ) == 1L) {
    return(
      modal_values
    )
  }

  tied_rows <- which(
    values %in%
      modal_values
  )

  values[tied_rows[which.min(source_order[tied_rows])]]
}

temperature_site_year <- temperature_raw %>%
  group_by(
    year,
    long,
    lat,
    site_key,
    site_year_key
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
    first_source_order = min(
      source_order
    ),
    .groups = "drop"
  ) %>%
  arrange(
    year,
    long,
    lat
  )

duplicate_diagnostics <- temperature_site_year %>%
  summarise(
    number_of_source_records = nrow(
      temperature_raw
    ),
    number_of_site_years = n(),
    duplicated_site_years = sum(
      number_of_source_records >
        1L
    ),
    conflicting_site_years = sum(
      number_of_distinct_levels >
        1L
    ),
    maximum_records_at_one_site_year = max(
      number_of_source_records
    ),
    maximum_distinct_levels_at_one_site_year = max(
      number_of_distinct_levels
    )
  )

eligible_year_table <- temperature_site_year %>%
  count(
    year,
    name = "number_of_distinct_sites"
  ) %>%
  filter(
    number_of_distinct_sites >=
      minimum_distinct_sites_per_year
  ) %>%
  arrange(
    year
  )

eligible_years <- eligible_year_table$year

if (length(
  eligible_years
) == 0L) {
  stop(
    "No event years contain at least ",
    minimum_distinct_sites_per_year,
    " distinct REACHES sites."
  )
}

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

readr::write_csv(
  eligible_year_table,
  file.path(
    intermediate_output_dir,
    "eligible_event_years.csv"
  )
)

message(
  "Eligible validation years: ",
  length(
    eligible_years
  )
)

# ------------------------------------------------------------
# 3. Hold-out splitting
# ------------------------------------------------------------

split_validation_data <- function(
    raw_data,
    site_year_data,
    eligible_years,
    repetition_seed) {

  set.seed(
    repetition_seed
  )

  held_out_site_years <- site_year_data %>%
    filter(
      year %in%
        eligible_years
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
    ) %>%
    mutate(
      test_record_id = row_number()
    )

  if (!identical(
    held_out_site_years$year,
    as.integer(
      eligible_years
    )
  )) {
    stop(
      "The held-out split does not contain exactly one site ",
      "for every eligible event year."
    )
  }

  held_out_keys <- held_out_site_years %>%
    dplyr::select(
      site_year_key
    )

  training_site_year <- anti_join(
    site_year_data,
    held_out_keys,
    by = "site_year_key"
  )

  training_raw <- anti_join(
    raw_data,
    held_out_keys,
    by = "site_year_key"
  )

  held_out_raw <- semi_join(
    raw_data,
    held_out_keys,
    by = "site_year_key"
  )

  leakage_count <- inner_join(
    training_raw %>%
      distinct(
        site_year_key
      ),
    held_out_keys,
    by = "site_year_key"
  ) %>%
    nrow()

  if (leakage_count != 0L) {
    stop(
      "At least one held-out site-year remained in the ",
      "training source records."
    )
  }

  list(
    training_raw = training_raw,
    training_site_year = training_site_year,
    held_out_site_years = held_out_site_years,
    held_out_raw = held_out_raw
  )
}

# ------------------------------------------------------------
# 4. Pooled pilot variogram
# ------------------------------------------------------------

estimate_pilot_variogram <- function(
    training_site_year) {

  spatial_data <- training_site_year %>%
    dplyr::select(
      level,
      year,
      long,
      lat
    ) %>%
    mutate(
      long_index = round(
        long,
        4
      ),
      lat_index = round(
        lat,
        4
      )
    )

  coordinate_data <- spatial_data %>%
    distinct(
      long_index,
      lat_index
    ) %>%
    transmute(
      long = long_index,
      lat = lat_index
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
      spatial_data$year
    )
  )

  time_index <- zoo::as.yearmon(
    years
  )

  location_keys <- paste(
    spatial_data$long_index,
    spatial_data$lat_index,
    sep = "_"
  )

  coordinate_keys <- paste(
    coordinate_data@coords[, 1],
    coordinate_data@coords[, 2],
    sep = "_"
  )

  index_matrix <- cbind(
    match(
      location_keys,
      coordinate_keys
    ),
    match(
      spatial_data$year,
      years
    )
  )

  if (anyNA(
    index_matrix
  )) {
    stop(
      "At least one training site-year could not be mapped ",
      "to the pooled variogram index."
    )
  }

  sts_data <- spacetime::STSDF(
    sp = coordinate_data,
    time = time_index,
    data = spatial_data %>%
      dplyr::select(
        level,
        year,
        long,
        lat
      ),
    index = index_matrix
  )

  pooled_variogram_st <- gstat::variogramST(
    level ~ 1,
    data = sts_data,
    tlags = 0,
    width = 10,
    na.omit = TRUE
  )

  pooled_variogram <- pooled_variogram_st %>%
    as.data.frame() %>%
    filter(
      dist > 1,
      np > 0,
      is.finite(
        dist
      ),
      is.finite(
        gamma
      )
    ) %>%
    transmute(
      np = np,
      dist = dist,
      gamma = gamma,
      dir.hor = 0,
      dir.ver = 0,
      id = "var1"
    )

  if (nrow(
    pooled_variogram
  ) < 3L) {
    stop(
      "Too few empirical variogram bins remain after cleaning."
    )
  }

  class(
    pooled_variogram
  ) <- c(
    "gstatVariogram",
    "data.frame"
  )

  variogram_fit <- gstat::fit.variogram(
    pooled_variogram,
    gstat::vgm(
      model = "Exp",
      nugget = NA
    ),
    fit.kappa = TRUE,
    fit.method = 2
  )

  nugget_row <- which(
    variogram_fit$model ==
      "Nug"
  )

  exponential_row <- which(
    variogram_fit$model ==
      "Exp"
  )

  if (
    length(
      nugget_row
    ) != 1L ||
      length(
        exponential_row
      ) != 1L
  ) {
    stop(
      "The fitted pilot variogram did not contain exactly one ",
      "nugget and one exponential component."
    )
  }

  pilot_nugget <- variogram_fit$psill[
    nugget_row
  ]

  pilot_partial_sill <- variogram_fit$psill[
    exponential_row
  ]

  pilot_range_km <- variogram_fit$range[
    exponential_row
  ]

  if (
    any(
      !is.finite(
        c(
          pilot_nugget,
          pilot_partial_sill,
          pilot_range_km
        )
      )
    ) ||
      pilot_nugget < 0 ||
      pilot_partial_sill <= 0 ||
      pilot_range_km <= 0
  ) {
    stop(
      "The pilot variogram produced invalid parameters."
    )
  }

  list(
    pilot_nugget = pilot_nugget,
    pilot_partial_sill =
      pilot_partial_sill,
    pilot_total_variance =
      pilot_nugget +
      pilot_partial_sill,
    pilot_range_km =
      pilot_range_km,
    empirical_variogram =
      pooled_variogram,
    variogram_fit =
      variogram_fit
  )
}

# ------------------------------------------------------------
# 5. Corrected analytic covariance calibration
# ------------------------------------------------------------

round_to_reaches_index <- function(
    values) {

  rounded <- round(
    values,
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

ordinal_moments <- function(
    total_variance) {

  if (
    !is.finite(
      total_variance
    ) ||
      total_variance <= 0
  ) {
    stop(
      "total_variance must be positive."
    )
  }

  total_sd <- sqrt(
    total_variance
  )

  standardized_cuts <- ordinal_cuts /
    total_sd

  probabilities <- diff(
    pnorm(
      standardized_cuts
    )
  )

  mean_z <- sum(
    ordinal_values *
      probabilities
  )

  second_moment_z <- sum(
    ordinal_values^2 *
      probabilities
  )

  list(
    mean = mean_z,
    variance =
      second_moment_z -
      mean_z^2,
    probabilities =
      probabilities
  )
}

ordinal_covariance <- function(
    rho,
    total_variance,
    marginal_mean = NULL,
    marginal_variance = NULL) {

  if (
    !is.finite(
      rho
    ) ||
      rho < -1 ||
      rho > 1
  ) {
    stop(
      "rho must lie in [-1, 1]."
    )
  }

  if (
    is.null(
      marginal_mean
    ) ||
      is.null(
        marginal_variance
      )
  ) {
    marginal <- ordinal_moments(
      total_variance
    )

    marginal_mean <- marginal$mean
    marginal_variance <- marginal$variance
  }

  if (abs(
    rho
  ) < 1e-12) {
    return(
      0
    )
  }

  if (rho > 1 - 1e-10) {
    return(
      marginal_variance
    )
  }

  total_sd <- sqrt(
    total_variance
  )

  standardized_thresholds <- c(
    -1.5,
    -0.5,
    0.5
  ) /
    total_sd

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
    row_threshold in seq_along(
      standardized_thresholds
    )
  ) {
    for (
      column_threshold in seq_along(
        standardized_thresholds
      )
    ) {

      threshold_1 <- standardized_thresholds[
        row_threshold
      ]

      threshold_2 <- standardized_thresholds[
        column_threshold
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
          corr =
            correlation_matrix,
          algorithm =
            mvtnorm::TVPACK()
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
          row_threshold
        ] *
        exceedance_probabilities[
          column_threshold
        ]
    }
  }

  covariance_sum
}

monotone_inverse <- function(
    x,
    y,
    target,
    curve_name) {

  ordering <- order(
    x
  )

  x <- x[
    ordering
  ]

  y <- y[
    ordering
  ]

  isotonic_fit <- stats::isoreg(
    x,
    y
  )

  monotone_y <- isotonic_fit$yf

  compressed <- data.frame(
    x = x,
    y = monotone_y
  ) %>%
    group_by(
      y
    ) %>%
    summarise(
      x = mean(
        x
      ),
      .groups = "drop"
    ) %>%
    arrange(
      y
    )

  if (
    target < min(
      compressed$y
    ) ||
      target > max(
        compressed$y
      )
  ) {
    stop(
      "The target for ",
      curve_name,
      " lies outside the calibration curve. Target = ",
      target,
      "; curve range = [",
      min(
        compressed$y
      ),
      ", ",
      max(
        compressed$y
      ),
      "]."
    )
  }

  stats::approx(
    x = compressed$y,
    y = compressed$x,
    xout = target,
    ties = "ordered"
  )$y
}

build_f1_curve <- function(
    target_observed_variance,
    number_of_grid_points) {

  maximum_variance <- max(
    2,
    2 *
      target_observed_variance
  )

  repeat {

    total_variance_grid <- seq(
      1e-4,
      maximum_variance,
      length.out =
        number_of_grid_points
    )

    observed_variance <- vapply(
      total_variance_grid,
      function(total_variance) {
        ordinal_moments(
          total_variance
        )$variance
      },
      numeric(1)
    )

    if (
      max(
        observed_variance
      ) >=
        target_observed_variance
    ) {
      break
    }

    maximum_variance <- 1.5 *
      maximum_variance

    if (maximum_variance > 20) {
      stop(
        "Unable to bracket the f1 inverse."
      )
    }
  }

  data.frame(
    total_variance =
      total_variance_grid,
    observed_variance =
      observed_variance
  )
}

build_f2_curve <- function(
    calibrated_total_variance,
    number_of_grid_points) {

  nugget_variance_grid <- seq(
    0,
    calibrated_total_variance,
    length.out =
      number_of_grid_points
  )

  marginal <- ordinal_moments(
    calibrated_total_variance
  )

  observed_nugget <- vapply(
    nugget_variance_grid,
    function(nugget_variance) {

      if (nugget_variance <= 0) {
        return(
          0
        )
      }

      latent_process_variance <-
        calibrated_total_variance -
        nugget_variance

      rho <- latent_process_variance /
        calibrated_total_variance

      marginal$variance -
        ordinal_covariance(
          rho = rho,
          total_variance =
            calibrated_total_variance,
          marginal_mean =
            marginal$mean,
          marginal_variance =
            marginal$variance
        )
    },
    numeric(1)
  )

  data.frame(
    nugget_variance =
      nugget_variance_grid,
    observed_nugget =
      observed_nugget
  )
}

build_covariance_lookup <- function(
    calibrated_total_variance,
    calibrated_process_variance,
    number_of_rho_values) {

  marginal <- ordinal_moments(
    calibrated_total_variance
  )

  maximum_rho <- calibrated_process_variance /
    calibrated_total_variance

  rho_grid <- seq(
    0,
    maximum_rho,
    length.out =
      number_of_rho_values
  )

  covariance_grid <- vapply(
    rho_grid,
    function(rho) {
      ordinal_covariance(
        rho = rho,
        total_variance =
          calibrated_total_variance,
        marginal_mean =
          marginal$mean,
        marginal_variance =
          marginal$variance
      )
    },
    numeric(1)
  )

  isotonic_covariance <- stats::isoreg(
    rho_grid,
    covariance_grid
  )$yf

  list(
    maximum_rho = maximum_rho,
    covariance_function =
      stats::splinefun(
        x = rho_grid,
        y = isotonic_covariance,
        method = "monoH.FC"
      ),
    lookup_table = data.frame(
      rho = rho_grid,
      ordinal_covariance =
        isotonic_covariance
    )
  )
}

fit_latent_range <- function(
    observed_range,
    target_sill,
    calibrated_process_variance,
    covariance_lookup,
    distances_km,
    range_bounds_km) {

  target_covariance <- target_sill *
    exp(
      -distances_km /
        observed_range
    )

  objective <- function(
      latent_range_km) {

    rho_by_distance <-
      covariance_lookup$maximum_rho *
      exp(
        -distances_km /
          latent_range_km
      )

    censored_covariance <-
      covariance_lookup$
        covariance_function(
          rho_by_distance
        )

    sum(
      (
        censored_covariance -
          target_covariance
      )^2
    )
  }

  stats::optimize(
    f = objective,
    interval =
      range_bounds_km
  )$minimum
}

build_f3_curve <- function(
    alpha_grid_km,
    target_sill,
    calibrated_process_variance,
    covariance_lookup,
    distances_km,
    range_bounds_km,
    target_pilot_range_km) {

  current_grid <- alpha_grid_km

  repeat {

    calibrated_latent_range <- vapply(
      current_grid,
      fit_latent_range,
      numeric(1),
      target_sill =
        target_sill,
      calibrated_process_variance =
        calibrated_process_variance,
      covariance_lookup =
        covariance_lookup,
      distances_km =
        distances_km,
      range_bounds_km =
        range_bounds_km
    )

    if (
      target_pilot_range_km >=
        min(
          calibrated_latent_range
        ) &&
        target_pilot_range_km <=
          max(
            calibrated_latent_range
          )
    ) {
      break
    }

    new_maximum <- 1.5 *
      max(
        current_grid
      )

    if (new_maximum > 5000) {
      stop(
        "Unable to bracket the f3 inverse."
      )
    }

    additional_grid <- seq(
      max(
        current_grid
      ) +
        10,
      new_maximum,
      by = 10
    )

    current_grid <- c(
      current_grid,
      additional_grid
    )
  }

  data.frame(
    alpha_km =
      current_grid,
    f3_alpha_km =
      calibrated_latent_range
  )
}

calibrate_covariance_parameters <- function(
    pilot_nugget,
    pilot_partial_sill,
    pilot_range_km) {

  pilot_total_variance <-
    pilot_nugget +
    pilot_partial_sill

  f1_data <- build_f1_curve(
    target_observed_variance =
      pilot_total_variance,
    number_of_grid_points =
      number_of_f1_grid_points
  )

  calibrated_total_variance <-
    monotone_inverse(
      x =
        f1_data$total_variance,
      y =
        f1_data$observed_variance,
      target =
        pilot_total_variance,
      curve_name = "f1"
    )

  f2_data <- build_f2_curve(
    calibrated_total_variance =
      calibrated_total_variance,
    number_of_grid_points =
      number_of_f2_grid_points
  )

  calibrated_nugget_variance <-
    monotone_inverse(
      x =
        f2_data$nugget_variance,
      y =
        f2_data$observed_nugget,
      target =
        pilot_nugget,
      curve_name = "f2"
    )

  calibrated_process_variance <-
    calibrated_total_variance -
    calibrated_nugget_variance

  if (
    !is.finite(
      calibrated_process_variance
    ) ||
      calibrated_process_variance <= 0
  ) {
    stop(
      "The calibrated process variance is not positive."
    )
  }

  covariance_lookup <-
    build_covariance_lookup(
      calibrated_total_variance =
        calibrated_total_variance,
      calibrated_process_variance =
        calibrated_process_variance,
      number_of_rho_values =
        rho_lookup_size
    )

  # Equation (5) targets the pilot partial sill.
  f3_data <- build_f3_curve(
    alpha_grid_km =
      alpha_grid_initial_km,
    target_sill =
      pilot_partial_sill,
    calibrated_process_variance =
      calibrated_process_variance,
    covariance_lookup =
      covariance_lookup,
    distances_km =
      distance_grid_km,
    range_bounds_km =
      latent_range_bounds_km,
    target_pilot_range_km =
      pilot_range_km
  )

  calibrated_range_km <-
    monotone_inverse(
      x =
        f3_data$alpha_km,
      y =
        f3_data$f3_alpha_km,
      target =
        pilot_range_km,
      curve_name = "f3"
    )

  marginal <- ordinal_moments(
    calibrated_total_variance
  )

  list(
    sigma_epsilon2 =
      calibrated_nugget_variance,
    sigma_Y2 =
      calibrated_process_variance,
    total_variance =
      calibrated_total_variance,
    range_km =
      calibrated_range_km,
    mean_z =
      marginal$mean,
    variance_z =
      marginal$variance,
    covariance_lookup =
      covariance_lookup,
    f1_data =
      f1_data,
    f2_data =
      f2_data,
    f3_data =
      f3_data
  )
}

# ------------------------------------------------------------
# 6. Interval-censored best linear prediction
# ------------------------------------------------------------

EZstar_h <- function(
    sigma,
    cuts = ordinal_cuts,
    values = ordinal_values) {

  lower_standardized <- head(
    cuts,
    -1L
  ) /
    sigma

  upper_standardized <- tail(
    cuts,
    -1L
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

solve_symmetric_system <- function(
    covariance_matrix,
    right_hand_side) {

  covariance_matrix <- (
    covariance_matrix +
      t(
        covariance_matrix
      )
  ) /
    2

  for (jitter in solver_jitter_values) {

    current_matrix <-
      covariance_matrix

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
          jitter = jitter,
          reciprocal_condition_number =
            rcond(
              current_matrix
            )
        )
      )
    }
  }

  stop(
    "Unable to obtain a positive-definite ordinal covariance ",
    "matrix after numerical jitter. Reciprocal condition number = ",
    signif(
      rcond(
        covariance_matrix
      ),
      6
    ),
    "."
  )
}

predict_held_out_sites <- function(
    training_raw,
    held_out_site_years,
    calibrated_parameters) {

  total_variance <-
    calibrated_parameters$
      total_variance

  sigma_total <- sqrt(
    total_variance
  )

  mean_z <-
    calibrated_parameters$
      mean_z

  variance_z <-
    calibrated_parameters$
      variance_z

  ezstar_h <- EZstar_h(
    sigma = sigma_total
  )

  covariance_lookup_function <-
    calibrated_parameters$
      covariance_lookup$
      covariance_function

  prediction_results <- vector(
    "list",
    nrow(
      held_out_site_years
    )
  )

  for (
    test_index in seq_len(
      nrow(
        held_out_site_years
      )
    )
  ) {

    held_out_row <-
      held_out_site_years[
        test_index,
        ,
        drop = FALSE
      ]

    current_year <-
      held_out_row$year

    training_year <- training_raw %>%
      filter(
        year ==
          current_year
      )

    if (nrow(
      training_year
    ) < 1L) {
      stop(
        "No training source records remain for year ",
        current_year,
        "."
      )
    }

    if (any(
      training_year$site_year_key ==
        held_out_row$site_year_key
    )) {
      stop(
        "The held-out site-year remains in the prediction ",
        "training set for year ",
        current_year,
        "."
      )
    }

    observation_coordinates <- as.matrix(
      training_year[
        ,
        c(
          "long",
          "lat"
        )
      ]
    )

    prediction_coordinate <- as.numeric(
      held_out_row[
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

    # Great-circle distances in km.
    observation_distance_matrix <-
      sp::spDists(
        observation_coordinates,
        longlat = TRUE
      )

    rho_matrix <- (
      calibrated_parameters$
        sigma_Y2 *
        exp(
          -observation_distance_matrix /
            calibrated_parameters$
              range_km
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

    # Only true self-covariances receive Var(Z).
    diag(
      SigmaZ
    ) <- variance_z

    observation_to_prediction_distance <-
      as.numeric(
        sp::spDists(
          x =
            observation_coordinates,
          y =
            matrix(
              prediction_coordinate,
              nrow = 1L
            ),
          longlat = TRUE
        )
      )

    latent_cross_covariance <-
      calibrated_parameters$
        sigma_Y2 *
      exp(
        -observation_to_prediction_distance /
          calibrated_parameters$
            range_km
      )

    cZY <- (
      latent_cross_covariance /
        total_variance
    ) *
      ezstar_h

    right_hand_side <- cbind(
      observed_categories -
        mean_z,
      cZY
    )

    solved <- solve_symmetric_system(
      covariance_matrix = SigmaZ,
      right_hand_side =
        right_hand_side
    )

    mean_weights <-
      solved$solution[
        ,
        1L
      ]

    variance_weights <-
      solved$solution[
        ,
        2L
      ]

    prediction_mean <- drop(
      crossprod(
        cZY,
        mean_weights
      )
    )

    prediction_mspe <-
      calibrated_parameters$
        sigma_Y2 -
      drop(
        crossprod(
          cZY,
          variance_weights
        )
      )

    prediction_results[[test_index]] <- data.frame(
      test_record_id =
        held_out_row$
          test_record_id,
      year =
        current_year,
      long =
        held_out_row$long,
      lat =
        held_out_row$lat,
      site_key =
        held_out_row$site_key,
      site_year_key =
        held_out_row$
          site_year_key,
      Z_obs =
        held_out_row$level,
      held_out_source_records =
        held_out_row$
          number_of_source_records,
      held_out_distinct_levels =
        held_out_row$
          number_of_distinct_levels,
      number_of_training_source_records =
        nrow(
          training_year
        ),
      number_of_training_distinct_sites =
        n_distinct(
          training_year$
            site_key
        ),
      yhat =
        prediction_mean,
      mspe =
        pmax(
          prediction_mspe,
          0
        ),
      solver_jitter =
        solved$jitter,
      reciprocal_condition_number =
        solved$
          reciprocal_condition_number
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

  predictive_variance <-
    test_predictions$mspe +
    sigma_epsilon2

  if (
    any(!is.finite(
      predictive_variance
    )) ||
      any(
        predictive_variance <= 0
      )
  ) {
    stop(
      "At least one held-out predictive variance is invalid."
    )
  }

  predictive_sd <- sqrt(
    predictive_variance
  )

  category_probabilities <- vapply(
    seq_along(
      values
    ),
    function(category_index) {

      pnorm(
        cuts[
          category_index +
            1L
        ],
        mean =
          test_predictions$yhat,
        sd =
          predictive_sd
      ) -
        pnorm(
          cuts[
            category_index
          ],
          mean =
            test_predictions$yhat,
          sd =
            predictive_sd
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

  probability_row_sums <- rowSums(
    category_probabilities
  )

  maximum_probability_sum_error <- max(
    abs(
      probability_row_sums -
        1
    )
  )

  if (
    !is.finite(
      maximum_probability_sum_error
    ) ||
      maximum_probability_sum_error >
        1e-8
  ) {
    stop(
      "Held-out ordinal category probabilities do not sum to ",
      "one. Maximum absolute error = ",
      maximum_probability_sum_error,
      "."
    )
  }

  category_probabilities <- pmax(
    category_probabilities,
    0
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
      cumulative_probabilities[
        ,
        seq_len(
          length(
            values
          ) -
            1L
        ),
        drop = FALSE
      ] -
        observed_cumulative_indicators[
          ,
          seq_len(
            length(
              values
            ) -
              1L
          ),
          drop = FALSE
        ]
    )^2
  )

  predicted_category <- values[
    max.col(
      category_probabilities,
      ties.method = "first"
    )
  ]

  category_index <- match(
    test_predictions$Z_obs,
    values
  )

  observed_lower <- cuts[
    category_index
  ]

  observed_upper <- cuts[
    category_index +
      1L
  ]

  observed_probability <- category_probabilities[
    cbind(
      seq_len(
        nrow(
          test_predictions
        )
      ),
      category_index
    )
  ]

  negative_log_score <- -log(
    pmax(
      observed_probability,
      .Machine$double.xmin
    )
  )

  compatibility_columns <- list()
  compatibility_summary_rows <- list()

  for (
    level_index in seq_along(
      nominal_levels
    )
  ) {

    nominal_level <- nominal_levels[
      level_index
    ]

    z_value <- qnorm(
      (
        1 +
          nominal_level
      ) /
        2
    )

    predictive_lower <-
      test_predictions$yhat -
      z_value *
        predictive_sd

    predictive_upper <-
      test_predictions$yhat +
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

    level_label <- paste0(
      "compatible_",
      formatC(
        100 *
          nominal_level,
        format = "f",
        digits = 0
      )
    )

    compatibility_columns[[level_label]] <- compatible

    compatibility_summary_rows[[level_index]] <- data.frame(
      nominal_level =
        nominal_level,
      compatibility = mean(
        compatible
      )
    )
  }

  observation_metrics <- bind_cols(
    data.frame(
      predictive_sd =
        predictive_sd,
      observed_interval_lower =
        observed_lower,
      observed_interval_upper =
        observed_upper,
      rps =
        rps_by_observation,
      predicted_category =
        predicted_category,
      observed_category_probability =
        observed_probability,
      negative_log_score =
        negative_log_score
    ),
    as.data.frame(
      compatibility_columns
    ),
    as.data.frame(
      category_probabilities,
      check.names = FALSE
    )
  )

  probability_columns <- paste0(
    "probability_",
    make.names(
      colnames(
        category_probabilities
      )
    )
  )

  probability_start <- ncol(
    observation_metrics
  ) -
    ncol(
      category_probabilities
    ) +
    1L

  names(
    observation_metrics
  )[
    probability_start:ncol(
      observation_metrics
    )
  ] <- probability_columns

  list(
    compatibility =
      bind_rows(
        compatibility_summary_rows
      ),
    rps = mean(
      rps_by_observation
    ),
    accuracy = mean(
      predicted_category ==
        test_predictions$Z_obs
    ),
    mean_negative_log_score = mean(
      negative_log_score
    ),
    maximum_probability_sum_error =
      maximum_probability_sum_error,
    observation_metrics =
      observation_metrics
  )
}

# ------------------------------------------------------------
# 8. Checkpoint helpers
# ------------------------------------------------------------

validation_settings <- list(
  analysis_years =
    analysis_years,
  number_of_repetitions =
    number_of_repetitions,
  minimum_distinct_sites_per_year =
    minimum_distinct_sites_per_year,
  nominal_levels =
    nominal_levels,
  master_seed =
    master_seed,
  number_of_f1_grid_points =
    number_of_f1_grid_points,
  number_of_f2_grid_points =
    number_of_f2_grid_points,
  distance_grid_km =
    distance_grid_km,
  alpha_grid_initial_km =
    alpha_grid_initial_km,
  latent_range_bounds_km =
    latent_range_bounds_km,
  rho_lookup_size =
    rho_lookup_size
)

empty_result_lists <- function() {
  list(
    split_results = vector(
      "list",
      number_of_repetitions
    ),
    diagnostic_results = vector(
      "list",
      number_of_repetitions
    ),
    observation_results = vector(
      "list",
      number_of_repetitions
    ),
    split_site_results = vector(
      "list",
      number_of_repetitions
    )
  )
}

result_lists <- empty_result_lists()
first_repetition <- 1L

if (
  resume_from_checkpoint &&
    file.exists(
      checkpoint_file
    )
) {

  checkpoint <- readRDS(
    checkpoint_file
  )

  if (!identical(
    checkpoint$settings,
    validation_settings
  )) {
    stop(
      "The saved interval-compatibility checkpoint was created ",
      "with different settings. Delete it before rerunning."
    )
  }

  completed_repetitions <-
    checkpoint$
      completed_repetitions

  result_lists$split_results[
    seq_len(
      completed_repetitions
    )
  ] <- split(
    checkpoint$
      split_results,
    checkpoint$
      split_results$
      repetition
  )

  result_lists$diagnostic_results[
    seq_len(
      completed_repetitions
    )
  ] <- split(
    checkpoint$
      diagnostic_results,
    checkpoint$
      diagnostic_results$
      repetition
  )

  result_lists$observation_results[
    seq_len(
      completed_repetitions
    )
  ] <- split(
    checkpoint$
      observation_results,
    checkpoint$
      observation_results$
      repetition
  )

  result_lists$split_site_results[
    seq_len(
      completed_repetitions
    )
  ] <- split(
    checkpoint$
      split_site_results,
    checkpoint$
      split_site_results$
      repetition
  )

  first_repetition <-
    completed_repetitions +
    1L

  message(
    "Resuming after repetition ",
    completed_repetitions,
    "."
  )
}

# ------------------------------------------------------------
# 9. Run repeated validation
# ------------------------------------------------------------

start_time <- Sys.time()

if (first_repetition <= number_of_repetitions) {

  for (
    repetition in seq.int(
      first_repetition,
      number_of_repetitions
    )
  ) {

    repetition_start_time <- Sys.time()

    repetition_seed <-
      master_seed +
      repetition -
      1L

    message(
      "\nValidation repetition ",
      repetition,
      "/",
      number_of_repetitions,
      " (seed ",
      repetition_seed,
      ")"
    )

    split_data <- split_validation_data(
      raw_data =
        temperature_raw,
      site_year_data =
        temperature_site_year,
      eligible_years =
        eligible_years,
      repetition_seed =
        repetition_seed
    )

    pilot_fit <- estimate_pilot_variogram(
      training_site_year =
        split_data$
          training_site_year
    )

    calibrated_parameters <-
      calibrate_covariance_parameters(
        pilot_nugget =
          pilot_fit$
            pilot_nugget,
        pilot_partial_sill =
          pilot_fit$
            pilot_partial_sill,
        pilot_range_km =
          pilot_fit$
            pilot_range_km
      )

    held_out_predictions <-
      predict_held_out_sites(
        training_raw =
          split_data$
            training_raw,
        held_out_site_years =
          split_data$
            held_out_site_years,
        calibrated_parameters =
          calibrated_parameters
      )

    validation_metrics <-
      compute_validation_metrics(
        test_predictions =
          held_out_predictions,
        sigma_epsilon2 =
          calibrated_parameters$
            sigma_epsilon2,
        nominal_levels =
          nominal_levels
      )

    repetition_elapsed_seconds <-
      as.numeric(
        difftime(
          Sys.time(),
          repetition_start_time,
          units = "secs"
        )
      )

    result_lists$split_results[[repetition]] <- validation_metrics$
      compatibility %>%
      mutate(
        repetition =
          repetition,
        repetition_seed =
          repetition_seed,
        number_of_held_out_sites =
          nrow(
            held_out_predictions
          )
      ) %>%
      dplyr::select(
        repetition,
        repetition_seed,
        nominal_level,
        compatibility,
        number_of_held_out_sites
      )

    result_lists$diagnostic_results[[repetition]] <- data.frame(
      repetition =
        repetition,
      repetition_seed =
        repetition_seed,
      number_of_eligible_years =
        length(
          eligible_years
        ),
      number_of_held_out_sites =
        nrow(
          held_out_predictions
        ),
      number_of_training_site_years =
        nrow(
          split_data$
            training_site_year
        ),
      number_of_training_source_records =
        nrow(
          split_data$
            training_raw
        ),
      rps =
        validation_metrics$rps,
      accuracy =
        validation_metrics$accuracy,
      mean_negative_log_score =
        validation_metrics$
          mean_negative_log_score,
      maximum_probability_sum_error =
        validation_metrics$
          maximum_probability_sum_error,
      pilot_nugget =
        pilot_fit$
          pilot_nugget,
      pilot_partial_sill =
        pilot_fit$
          pilot_partial_sill,
      pilot_total_variance =
        pilot_fit$
          pilot_total_variance,
      pilot_range_km =
        pilot_fit$
          pilot_range_km,
      calibrated_sigma_epsilon2 =
        calibrated_parameters$
          sigma_epsilon2,
      calibrated_sigma_Y2 =
        calibrated_parameters$
          sigma_Y2,
      calibrated_total_variance =
        calibrated_parameters$
          total_variance,
      calibrated_range_km =
        calibrated_parameters$
          range_km,
      minimum_solver_rcond = min(
        held_out_predictions$
          reciprocal_condition_number
      ),
      maximum_solver_jitter = max(
        held_out_predictions$
          solver_jitter
      ),
      repetition_elapsed_seconds =
        repetition_elapsed_seconds
    )

    result_lists$observation_results[[repetition]] <- bind_cols(
      data.frame(
        repetition =
          repetition,
        repetition_seed =
          repetition_seed
      ),
      held_out_predictions,
      validation_metrics$
        observation_metrics
    )

    result_lists$split_site_results[[repetition]] <- split_data$
      held_out_site_years %>%
      mutate(
        repetition =
          repetition,
        repetition_seed =
          repetition_seed
      ) %>%
      dplyr::select(
        repetition,
        repetition_seed,
        test_record_id,
        year,
        long,
        lat,
        site_key,
        site_year_key,
        level,
        number_of_source_records,
        number_of_distinct_levels
      )

    checkpoint <- list(
      settings =
        validation_settings,
      completed_repetitions =
        repetition,
      split_results = bind_rows(
        result_lists$
          split_results[
            seq_len(
              repetition
            )
          ]
      ),
      diagnostic_results = bind_rows(
        result_lists$
          diagnostic_results[
            seq_len(
              repetition
            )
          ]
      ),
      observation_results = bind_rows(
        result_lists$
          observation_results[
            seq_len(
              repetition
            )
          ]
      ),
      split_site_results = bind_rows(
        result_lists$
          split_site_results[
            seq_len(
              repetition
            )
          ]
      )
    )

    saveRDS(
      checkpoint,
      checkpoint_file
    )
  }
}

compatibility_by_split <- bind_rows(
  result_lists$
    split_results
)

diagnostics_by_split <- bind_rows(
  result_lists$
    diagnostic_results
)

held_out_observations <- bind_rows(
  result_lists$
    observation_results
)

split_sites <- bind_rows(
  result_lists$
    split_site_results
)

# ------------------------------------------------------------
# 10. Summaries and outputs
# ------------------------------------------------------------

compatibility_summary <-
  compatibility_by_split %>%
  group_by(
    nominal_level
  ) %>%
  summarise(
    repetitions = n(),
    held_out_observations_per_repetition =
      first(
        number_of_held_out_sites
      ),
    mean_compatibility = mean(
      compatibility
    ),
    sd_compatibility = sd(
      compatibility
    ),
    minimum_compatibility = min(
      compatibility
    ),
    maximum_compatibility = max(
      compatibility
    ),
    pooled_compatibility = {
      compatibility_column <- paste0(
        "compatible_",
        formatC(
          100 *
            first(
              nominal_level
            ),
          format = "f",
          digits = 0
        )
      )

      mean(
        held_out_observations[[compatibility_column]]
      )
    },
    .groups = "drop"
  ) %>%
  arrange(
    desc(
      nominal_level
    )
  )

compatibility_summary_formatted <-
  compatibility_summary %>%
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
    ),
    `Minimum` = sprintf(
      "%.3f",
      minimum_compatibility
    ),
    `Maximum` = sprintf(
      "%.3f",
      maximum_compatibility
    ),
    `Pooled compatibility` = sprintf(
      "%.3f",
      pooled_compatibility
    )
  )

diagnostic_summary <- data.frame(
  Metric = c(
    "Ranked probability score",
    "Classification accuracy",
    "Negative log score"
  ),
  Mean = c(
    mean(
      diagnostics_by_split$rps
    ),
    mean(
      diagnostics_by_split$accuracy
    ),
    mean(
      diagnostics_by_split$
        mean_negative_log_score
    )
  ),
  SD = c(
    sd(
      diagnostics_by_split$rps
    ),
    sd(
      diagnostics_by_split$accuracy
    ),
    sd(
      diagnostics_by_split$
        mean_negative_log_score
    )
  )
)

calibration_summary <- diagnostics_by_split %>%
  summarise(
    across(
      c(
        pilot_nugget,
        pilot_partial_sill,
        pilot_total_variance,
        pilot_range_km,
        calibrated_sigma_epsilon2,
        calibrated_sigma_Y2,
        calibrated_total_variance,
        calibrated_range_km
      ),
      list(
        mean = mean,
        sd = sd,
        minimum = min,
        maximum = max
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
  calibration_summary,
  file.path(
    validation_output_dir,
    "interval_compatibility_calibration_summary.csv"
  )
)

readr::write_csv(
  held_out_observations,
  file.path(
    intermediate_output_dir,
    "interval_compatibility_held_out_predictions.csv.gz"
  )
)

readr::write_csv(
  split_sites,
  file.path(
    intermediate_output_dir,
    "interval_compatibility_split_sites.csv.gz"
  )
)

saveRDS(
  list(
    settings =
      validation_settings,
    duplicate_diagnostics =
      duplicate_diagnostics,
    eligible_year_table =
      eligible_year_table,
    compatibility_summary =
      compatibility_summary,
    compatibility_by_split =
      compatibility_by_split,
    diagnostics_by_split =
      diagnostics_by_split,
    diagnostic_summary =
      diagnostic_summary,
    calibration_summary =
      calibration_summary,
    held_out_observations =
      held_out_observations,
    split_sites =
      split_sites
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
      "Analysis years:",
      min(
        analysis_years
      ),
      "--",
      max(
        analysis_years
      ),
      "\n"
    )

    cat(
      "Eligible event years:",
      length(
        eligible_years
      ),
      "\n"
    )

    cat(
      "Minimum distinct sites per event year:",
      minimum_distinct_sites_per_year,
      "\n"
    )

    cat(
      "Repetitions:",
      number_of_repetitions,
      "\n"
    )

    cat(
      "Master seed:",
      master_seed,
      "\n"
    )

    cat(
      "Calibration engine: analytic Gaussian moments\n"
    )

    cat(
      "Range calibration: Equation (5), pilot partial sill target\n"
    )

    cat(
      "Spatial distance unit: great-circle kilometres\n"
    )

    cat(
      "Elapsed time:",
      round(
        elapsed_minutes,
        2
      ),
      "minutes\n\n"
    )

    cat(
      "Compatibility summary\n"
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
      "\nCalibration-parameter summary\n"
    )

    print(
      calibration_summary
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

if (
  remove_checkpoint_on_success &&
    file.exists(
      checkpoint_file
    )
) {
  unlink(
    checkpoint_file
  )
}

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
  "Held-out predictions saved to: ",
  file.path(
    intermediate_output_dir,
    "interval_compatibility_held_out_predictions.csv.gz"
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
