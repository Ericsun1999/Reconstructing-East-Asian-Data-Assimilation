here::i_am("Code/prepare_calibration.R")

# ============================================================
# Shared covariance-parameter calibration for Figures 4, 5,
# and 7(e).
#
# Manuscript definitions:
#   1. Fit a pooled exponential variogram with nugget while
#      temporarily ignoring interval censoring.
#   2. Calibrate total variance through f1.
#   3. Calibrate nugget variance through f2.
#   4. Calibrate the spatial range through f3 and Equation (5).
#
# Default numerical engine:
#   "analytic" evaluates the same Gaussian moments using the
#   univariate/bivariate normal formulas in Supplement S1.
#   This is much faster and deterministic.
#
# To follow the Monte Carlo wording in Supplement S2 literally,
# change calibration_engine to "monte_carlo".
#
# Main output:
#   Output/Intermediate/calibration_parameters.rds
# ============================================================

library(here)
library(readxl)
library(readr)
library(dplyr)
library(sp)
library(spacetime)
library(zoo)
library(gstat)
library(ggplot2)
library(mvtnorm)

# ------------------------------------------------------------
# 1. Numerical settings
# ------------------------------------------------------------

calibration_engine <- "monte_carlo"
stopifnot(
  calibration_engine %in% c(
    "analytic",
    "monte_carlo"
  )
)

master_seed <- 10L

# Number of values used to draw the calibration curves.
n_f1 <- 500L
n_f2 <- 200L

# Supplement S2 states that B around 10^6 produces a smooth
# Monte Carlo curve. This value is used only when
# calibration_engine == "monte_carlo".
mc_sample_size <- 1000000L

# f3 settings.
distance_grid <- seq(
  0,
  1200,
  by = 5
)

alpha_grid_initial <- seq(
  50,
  1000,
  by = 10
)

latent_range_bounds <- c(
  20,
  2500
)

rho_lookup_size <- 1501L

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

input_file <- here::here(
  "Data",
  "temperature index value.v1.xlsx"
)

output_directory <- here::here(
  "Output",
  "Intermediate"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

calibration_file <- file.path(
  output_directory,
  "calibration_parameters.rds"
)

# ------------------------------------------------------------
# 2. Read and clean REACHES observations
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop(
    "The REACHES input file was not found: ",
    input_file
  )
}

temperature_raw <- readxl::read_excel(
  input_file,
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

names(temperature_raw) <- c(
  "level",
  "year",
  "long",
  "lat"
)

temperature_raw <- temperature_raw %>%
  transmute(
    source_order = row_number(),
    level = as.numeric(level),
    year = as.integer(year),
    long = as.numeric(long),
    lat = as.numeric(lat)
  ) %>%
  filter(
    year >= 1368L,
    year <= 1911L,
    level %in% ordinal_values,
    is.finite(long),
    is.finite(lat)
  )

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

# The spatial model has one ordinal observation per site-year.
# Repeated documentary entries at the same site and year are
# therefore collapsed to their modal category.
temp2 <- temperature_raw %>%
  group_by(
    year,
    long,
    lat
  ) %>%
  summarise(
    level = mode_first(
      level,
      source_order
    ),
    number_of_source_records = n(),
    .groups = "drop"
  ) %>%
  dplyr::select(
    level,
    year,
    long,
    lat,
    number_of_source_records
  ) %>%
  arrange(
    year,
    long,
    lat
  )

event_years <- sort(
  unique(
    temp2$year
  )
)

if (length(event_years) != 524L) {
  warning(
    "The manuscript describes 524 REACHES event years, but ",
    length(event_years),
    " were found after cleaning."
  )
}

# ------------------------------------------------------------
# 3. Pooled pilot variogram
# ------------------------------------------------------------

estimate_pilot_variogram <- function(
    data) {

  spatial_data <- data %>%
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

  if (anyNA(index_matrix)) {
    stop(
      "At least one observation could not be mapped to the ",
      "spatio-temporal variogram index."
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
      is.finite(dist),
      is.finite(gamma)
    ) %>%
    transmute(
      np = np,
      dist = dist,
      gamma = gamma,
      dir.hor = 0,
      dir.ver = 0,
      id = "var1"
    )

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
    length(nugget_row) != 1L ||
      length(exponential_row) != 1L
  ) {
    stop(
      "The fitted variogram did not contain one nugget and ",
      "one exponential component."
    )
  }

  pilot_nugget <- variogram_fit$psill[nugget_row]

  pilot_partial_sill <- variogram_fit$psill[exponential_row]

  pilot_range <- variogram_fit$range[exponential_row]

  if (
    any(
      !is.finite(
        c(
          pilot_nugget,
          pilot_partial_sill,
          pilot_range
        )
      )
    ) ||
      pilot_nugget < 0 ||
      pilot_partial_sill <= 0 ||
      pilot_range <= 0
  ) {
    stop(
      "The pilot variogram produced invalid parameter values."
    )
  }

  spatial_output <- data %>%
    dplyr::select(
      level,
      year,
      long,
      lat
    )

  sp::coordinates(
    spatial_output
  ) <- ~ long + lat

  sp::proj4string(
    spatial_output
  ) <- sp::CRS(
    "+proj=longlat +datum=WGS84"
  )

  list(
    psill = c(
      pilot_nugget,
      pilot_partial_sill
    ),
    range = pilot_range,
    y2 = spatial_output,
    empirical_variogram = pooled_variogram,
    variogram_fit = variogram_fit
  )
}

var_fit2 <- estimate_pilot_variogram(
  temp2
)

pilot_nugget <- var_fit2$psill[1]
pilot_partial_sill <- var_fit2$psill[2]
pilot_total_variance <- sum(
  var_fit2$psill
)
pilot_range <- var_fit2$range

message(
  "Pilot estimates: range = ",
  signif(
    pilot_range,
    6
  ),
  ", partial sill = ",
  signif(
    pilot_partial_sill,
    6
  ),
  ", nugget = ",
  signif(
    pilot_nugget,
    6
  )
)

# ------------------------------------------------------------
# 4. Ordinal Gaussian moments
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
    !is.finite(total_variance) ||
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
    variance = second_moment_z -
      mean_z^2,
    probabilities = probabilities
  )
}

ordinal_covariance <- function(
    rho,
    total_variance,
    marginal_mean = NULL,
    marginal_variance = NULL) {

  if (
    !is.finite(rho) ||
      rho < -1 ||
      rho > 1
  ) {
    stop(
      "rho must lie in [-1, 1]."
    )
  }

  if (is.null(marginal_mean) ||
      is.null(marginal_variance)) {

    marginal <- ordinal_moments(
      total_variance
    )

    marginal_mean <- marginal$mean
    marginal_variance <- marginal$variance
  }

  if (abs(rho) < 1e-12) {
    return(0)
  }

  if (rho > 1 - 1e-10) {
    return(
      marginal_variance
    )
  }

  # For the four REACHES categories,
  #
  #   Z = -2 +
  #       I(W >= -1.5) +
  #       I(W >= -0.5) +
  #       I(W >=  0.5).
  #
  # The constant does not affect covariance, so Cov(Z1, Z2)
  # is a sum of nine threshold-indicator covariances. This
  # requires only bivariate normal CDF evaluations over
  # semi-infinite regions and is faster than evaluating all
  # sixteen category rectangles separately.

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

  for (row_threshold in seq_along(
    standardized_thresholds
  )) {
    for (column_threshold in seq_along(
      standardized_thresholds
    )) {

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
        exceedance_probabilities[row_threshold] *
        exceedance_probabilities[column_threshold]
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
    target < min(compressed$y) ||
      target > max(compressed$y)
  ) {
    stop(
      "The target for ",
      curve_name,
      " lies outside the computed calibration curve. ",
      "Target = ",
      target,
      "; curve range = [",
      min(compressed$y),
      ", ",
      max(compressed$y),
      "]."
    )
  }

  approx(
    x = compressed$y,
    y = compressed$x,
    xout = target,
    ties = "ordered"
  )$y
}

# ------------------------------------------------------------
# 5. f1: total latent variance -> observed ordinal variance
# ------------------------------------------------------------

build_f1_curve <- function(
    target_observed_variance,
    number_of_grid_points,
    engine,
    mc_size) {

  maximum_variance <- max(
    2,
    2 *
      target_observed_variance
  )

  repeat {

    total_variance_grid <- seq(
      1e-4,
      maximum_variance,
      length.out = number_of_grid_points
    )

    if (engine == "analytic") {

      observed_variance <- vapply(
        total_variance_grid,
        function(total_variance) {
          ordinal_moments(
            total_variance
          )$variance
        },
        numeric(1)
      )

    } else {

      standard_normal <- rnorm(
        mc_size
      )

      observed_variance <- vapply(
        total_variance_grid,
        function(total_variance) {

          ordinal_sample <- round_to_reaches_index(
            sqrt(
              total_variance
            ) *
              standard_normal
          )

          var(
            ordinal_sample
          )
        },
        numeric(1)
      )
    }

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
    total_variance = total_variance_grid,
    observed_variance = observed_variance
  )
}

# ------------------------------------------------------------
# 6. f2: latent nugget variance -> observed nugget
# ------------------------------------------------------------

build_f2_curve <- function(
    calibrated_total_variance,
    number_of_grid_points,
    engine,
    mc_size) {

  nugget_variance_grid <- seq(
    0,
    calibrated_total_variance,
    length.out = number_of_grid_points
  )

  marginal <- ordinal_moments(
    calibrated_total_variance
  )

  if (engine == "analytic") {

    observed_nugget <- vapply(
      nugget_variance_grid,
      function(nugget_variance) {

        if (nugget_variance <= 0) {
          return(0)
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

  } else {

    standard_y <- rnorm(
      mc_size
    )

    standard_epsilon_1 <- rnorm(
      mc_size
    )

    standard_epsilon_2 <- rnorm(
      mc_size
    )

    observed_nugget <- vapply(
      nugget_variance_grid,
      function(nugget_variance) {

        latent_process_variance <-
          calibrated_total_variance -
          nugget_variance

        latent_process <- sqrt(
          max(
            latent_process_variance,
            0
          )
        ) *
          standard_y

        epsilon_sd <- sqrt(
          max(
            nugget_variance,
            0
          )
        )

        z1 <- round_to_reaches_index(
          latent_process +
            epsilon_sd *
            standard_epsilon_1
        )

        z2 <- round_to_reaches_index(
          latent_process +
            epsilon_sd *
            standard_epsilon_2
        )

        var(
          z1 -
            z2
        ) /
          2
      },
      numeric(1)
    )
  }

  data.frame(
    nugget_variance = nugget_variance_grid,
    observed_nugget = observed_nugget
  )
}

# ------------------------------------------------------------
# 7. f3: observed exponential range -> latent range
# ------------------------------------------------------------

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
    length.out = number_of_rho_values
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
    covariance_function = stats::splinefun(
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
    calibrated_total_variance,
    calibrated_process_variance,
    covariance_lookup,
    distances,
    range_bounds) {

  target_covariance <- target_sill *
    exp(
      -distances /
        observed_range
    )

  objective <- function(
      latent_range) {

    rho_by_distance <-
      covariance_lookup$maximum_rho *
      exp(
        -distances /
          latent_range
      )

    censored_covariance <-
      covariance_lookup$covariance_function(
        rho_by_distance
      )

    sum(
      (
        censored_covariance -
          target_covariance
      )^2
    )
  }

  optimisation <- optimize(
    f = objective,
    interval = range_bounds
  )

  optimisation$minimum
}


build_f3_curve <- function(
    alpha_grid,
    target_sill,
    calibrated_total_variance,
    calibrated_process_variance,
    covariance_lookup,
    distances,
    range_bounds,
    target_pilot_range) {

  current_grid <- alpha_grid

  repeat {

    calibrated_latent_range <- vapply(
      current_grid,
      fit_latent_range,
      numeric(1),
      target_sill = target_sill,
      calibrated_total_variance =
        calibrated_total_variance,
      calibrated_process_variance =
        calibrated_process_variance,
      covariance_lookup =
        covariance_lookup,
      distances = distances,
      range_bounds = range_bounds
    )

    if (
      target_pilot_range >=
        min(
          calibrated_latent_range
        ) &&
        target_pilot_range <=
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
    alpha = current_grid,
    f3_alpha = calibrated_latent_range
  )
}

# ------------------------------------------------------------
# 8. Run the three calibrations
# ------------------------------------------------------------

set.seed(
  master_seed
)

start_time <- Sys.time()

f1_data <- build_f1_curve(
  target_observed_variance =
    pilot_total_variance,
  number_of_grid_points = n_f1,
  engine = calibration_engine,
  mc_size = mc_sample_size
)

calibrated_total_variance <- monotone_inverse(
  x = f1_data$total_variance,
  y = f1_data$observed_variance,
  target = pilot_total_variance,
  curve_name = "f1"
)

f2_data <- build_f2_curve(
  calibrated_total_variance =
    calibrated_total_variance,
  number_of_grid_points = n_f2,
  engine = calibration_engine,
  mc_size = mc_sample_size
)

calibrated_nugget_variance <- monotone_inverse(
  x = f2_data$nugget_variance,
  y = f2_data$observed_nugget,
  target = pilot_nugget,
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

covariance_lookup <- build_covariance_lookup(
  calibrated_total_variance =
    calibrated_total_variance,
  calibrated_process_variance =
    calibrated_process_variance,
  number_of_rho_values =
    rho_lookup_size
)

# Equation (5) uses the pilot partial sill in the target
# exponential covariance curve.
f3_data <- build_f3_curve(
  alpha_grid = alpha_grid_initial,
  target_sill = pilot_partial_sill,
  calibrated_total_variance =
    calibrated_total_variance,
  calibrated_process_variance =
    calibrated_process_variance,
  covariance_lookup = covariance_lookup,
  distances = distance_grid,
  range_bounds = latent_range_bounds,
  target_pilot_range = pilot_range
)

calibrated_range <- monotone_inverse(
  x = f3_data$alpha,
  y = f3_data$f3_alpha,
  target = pilot_range,
  curve_name = "f3"
)

# Legacy diagnostic:
# The original code used the calibrated latent process variance
# rather than the pilot partial sill in Equation (5). This curve
# is saved only to quantify the consequence of that difference.
f3_legacy_data <- build_f3_curve(
  alpha_grid = alpha_grid_initial,
  target_sill =
    calibrated_process_variance,
  calibrated_total_variance =
    calibrated_total_variance,
  calibrated_process_variance =
    calibrated_process_variance,
  covariance_lookup = covariance_lookup,
  distances = distance_grid,
  range_bounds = latent_range_bounds,
  target_pilot_range = pilot_range
)

legacy_calibrated_range <- monotone_inverse(
  x = f3_legacy_data$alpha,
  y = f3_legacy_data$f3_alpha,
  target = pilot_range,
  curve_name = "legacy f3"
)

elapsed_minutes <- as.numeric(
  difftime(
    Sys.time(),
    start_time,
    units = "mins"
  )
)

# ------------------------------------------------------------
# 9. Figure 4 plot objects
# ------------------------------------------------------------

plot1 <- ggplot(
  f1_data,
  aes(
    x = total_variance,
    y = observed_variance
  )
) +
  geom_line(
    colour = "firebrick",
    linewidth = 0.8
  ) +
  geom_point(
    data = data.frame(
      x = calibrated_total_variance,
      y = pilot_total_variance
    ),
    aes(
      x = x,
      y = y
    ),
    colour = "red",
    size = 3,
    inherit.aes = FALSE
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "blue"
  ) +
  geom_hline(
    yintercept = pilot_total_variance,
    linetype = "longdash",
    colour = "red"
  ) +
  labs(
    x = expression(
      sigma[Y]^2 +
        sigma[epsilon]^2
    ),
    y = expression(
      f[1](
        sigma[Y]^2 +
          sigma[epsilon]^2
      )
    )
  ) +
  theme_minimal(
    base_size = 14
  )

plot2 <- ggplot(
  f2_data,
  aes(
    x = nugget_variance,
    y = observed_nugget
  )
) +
  geom_line(
    colour = "firebrick",
    linewidth = 0.8
  ) +
  geom_point(
    data = data.frame(
      x = calibrated_nugget_variance,
      y = pilot_nugget
    ),
    aes(
      x = x,
      y = y
    ),
    colour = "red",
    size = 3,
    inherit.aes = FALSE
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "blue"
  ) +
  geom_hline(
    yintercept = pilot_nugget,
    linetype = "longdash",
    colour = "red"
  ) +
  labs(
    x = expression(
      sigma[epsilon]^2
    ),
    y = expression(
      f[2](
        sigma[epsilon]^2
      )
    )
  ) +
  theme_minimal(
    base_size = 14
  )

plot3 <- ggplot(
  f3_data,
  aes(
    x = alpha,
    y = f3_alpha
  )
) +
  geom_line(
    colour = "firebrick",
    linewidth = 0.8
  ) +
  geom_point(
    data = data.frame(
      x = calibrated_range,
      y = pilot_range
    ),
    aes(
      x = x,
      y = y
    ),
    colour = "red",
    size = 3,
    inherit.aes = FALSE
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "blue"
  ) +
  geom_hline(
    yintercept = pilot_range,
    linetype = "longdash",
    colour = "red"
  ) +
  labs(
    x = expression(
      alpha
    ),
    y = expression(
      f[3](
        alpha
      )
    )
  ) +
  theme_minimal(
    base_size = 14
  )

# ------------------------------------------------------------
# 10. Save results
# ------------------------------------------------------------

vario_fit2 <- list(
  psill1 = calibrated_nugget_variance,
  psill2 = calibrated_process_variance,
  range = calibrated_range,
  range_legacy = legacy_calibrated_range,
  total_variance = calibrated_total_variance,
  calibration_engine = calibration_engine,
  plot1 = plot1,
  plot2 = plot2,
  plot3 = plot3,
  calibration_data1 = f1_data,
  calibration_data2 = f2_data,
  calibration_data3 = f3_data,
  calibration_data3_legacy =
    f3_legacy_data,
  covariance_lookup =
    covariance_lookup$lookup_table
)

results <- list(
  temp2 = temp2,
  var_fit2 = var_fit2,
  vario_fit2 = vario_fit2,
  diagnostics = list(
    number_of_event_years =
      length(
        event_years
      ),
    pilot = c(
      range = pilot_range,
      process_variance =
        pilot_partial_sill,
      nugget_variance =
        pilot_nugget,
      total_variance =
        pilot_total_variance
    ),
    calibrated_manuscript = c(
      range = calibrated_range,
      process_variance =
        calibrated_process_variance,
      nugget_variance =
        calibrated_nugget_variance,
      total_variance =
        calibrated_total_variance
    ),
    calibrated_legacy_f3 = c(
      range = legacy_calibrated_range,
      process_variance =
        calibrated_process_variance,
      nugget_variance =
        calibrated_nugget_variance,
      total_variance =
        calibrated_total_variance
    ),
    elapsed_minutes =
      elapsed_minutes
  )
)

saveRDS(
  results,
  calibration_file
)

readr::write_csv(
  f1_data,
  file.path(
    output_directory,
    "calibration_f1.csv"
  )
)

readr::write_csv(
  f2_data,
  file.path(
    output_directory,
    "calibration_f2.csv"
  )
)

readr::write_csv(
  f3_data,
  file.path(
    output_directory,
    "calibration_f3_manuscript.csv"
  )
)

readr::write_csv(
  f3_legacy_data,
  file.path(
    output_directory,
    "calibration_f3_legacy.csv"
  )
)

message(
  "Calibration results saved to: ",
  calibration_file
)

message(
  "Calibrated manuscript parameters: range = ",
  signif(
    calibrated_range,
    6
  ),
  ", process variance = ",
  signif(
    calibrated_process_variance,
    6
  ),
  ", nugget variance = ",
  signif(
    calibrated_nugget_variance,
    6
  )
)

message(
  "Legacy-f3 range for comparison: ",
  signif(
    legacy_calibrated_range,
    6
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
