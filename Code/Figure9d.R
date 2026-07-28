here::i_am("Code/Figure9d.R")

# ============================================================
# Annual Kalman assimilation for:
#   Figure 9(d)  -- Beijing
#   Figure S3(d) -- Shanghai
#   Figure S4(d) -- Hong Kong
#
# Major implementation points:
#   1. The filter and smoother run on every year from 1368 to
#      1911. Years without REACHES information use a
#      prediction-only Kalman step.
#   2. Latent Gaussian variables are never rounded or clipped
#      to the observed REACHES category range [-2, 1].
#   3. The kriging SD is used in the quantile-mapping CDF, while
#      its square (MSPE) is used to obtain Var(Yhat_t).
#   4. The latent-process variance is read from the kriging
#      metadata rather than hard-coded.
#   5. REACHES, LME, and prior parameters are checked against
#      the same city coordinates and complete annual grid.
#
# Supported input modes:
#   "generated":
#     Output/Intermediate/REACHES/
#     Output/Intermediate/LME/
#     Output/Intermediate/Prior/
#
#   "precomputed":
#     Data/REACHES/precomputed/
#     Data/LME data/precomputed/
#     Data/par/
#
#   "auto" (default):
#     Use the complete generated set when available; otherwise
#     use the complete precomputed set.
#
# Main outputs:
#   Output/Figure9/Figure9d.png
#   Output/Supplementary/FigureS3d.png
#   Output/Supplementary/FigureS4d.png
#
# Numerical and diagnostic outputs:
#   Output/Intermediate/Assimilation/
#     assimilation_Beijing.csv
#     assimilation_Shanghai.csv
#     assimilation_HongKong.csv
#     assimilation_metrics.csv
#     assimilation_metadata.csv
#     diagnostic_Beijing.png
#     diagnostic_Shanghai.png
#     diagnostic_HongKong.png
#     Figure9d_sessionInfo.txt
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(np)
library(ggplot2)

# ------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------

analysis_years <- 1368:1911

input_mode <- "auto"

allowed_input_modes <- c(
  "auto",
  "generated",
  "precomputed"
)

if (!input_mode %in% allowed_input_modes) {
  stop(
    "input_mode must be one of: ",
    paste(
      allowed_input_modes,
      collapse = ", "
    )
  )
}

measurement_mc_size <- 10000L
measurement_mc_seed <- 1L

# The numerical representation of g(y) is constructed over a
# very wide latent-Gaussian range. Increase grid_size for a finer
# interpolation at additional computational cost.
mapping_grid_size <- 2001L
latent_tail_probability <- 1e-14
lme_cdf_grid_size <- 5001L

coordinate_tolerance <- 1e-8
variance_tolerance <- 1e-10
beta_tolerance <- 1e-10
identity_tolerance <- 1e-8

save_diagnostic_plots <- FALSE

figure9_dir <- here::here(
  "Output",
  "Figure9"
)

supplementary_dir <- here::here(
  "Output",
  "Supplementary"
)

assimilation_output_dir <- here::here(
  "Output",
  "Intermediate",
  "Assimilation"
)

dir.create(
  figure9_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  supplementary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  assimilation_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

city_config <- list(
  HongKong = list(
    code = "H",
    long = 113.75,
    lat = 22.25,
    figure_file = file.path(
      supplementary_dir,
      "FigureS4d.png"
    )
  ),
  Shanghai = list(
    code = "S",
    long = 121.25,
    lat = 31.25,
    figure_file = file.path(
      supplementary_dir,
      "FigureS3d.png"
    )
  ),
  Beijing = list(
    code = "B",
    long = 116.25,
    lat = 39.75,
    figure_file = file.path(
      figure9_dir,
      "Figure9d.png"
    )
  )
)

# ------------------------------------------------------------
# 2. Select one complete input set
# ------------------------------------------------------------

generated_inputs <- list(
  reaches_mean = here::here(
    "Output",
    "Intermediate",
    "REACHES",
    "reaches_kriging_city3_mean.csv"
  ),
  reaches_sd = here::here(
    "Output",
    "Intermediate",
    "REACHES",
    "reaches_kriging_city3_sd.csv"
  ),
  reaches_metadata = here::here(
    "Output",
    "Intermediate",
    "REACHES",
    "reaches_kriging_metadata.csv"
  ),
  lme_city3 = here::here(
    "Output",
    "Intermediate",
    "LME",
    "lme_city3_annual_1368_1911.csv"
  ),
  parameter_dir = here::here(
    "Output",
    "Intermediate",
    "Prior"
  )
)

precomputed_inputs <- list(
  reaches_mean = here::here(
    "Data",
    "REACHES",
    "precomputed",
    "reaches_kriging_city3_mean.csv"
  ),
  reaches_sd = here::here(
    "Data",
    "REACHES",
    "precomputed",
    "reaches_kriging_city3_sd.csv"
  ),
  reaches_metadata = here::here(
    "Data",
    "REACHES",
    "precomputed",
    "reaches_kriging_metadata.csv"
  ),
  lme_city3 = here::here(
    "Data",
    "LME data",
    "precomputed",
    "lme_city3_annual_1368_1911.csv"
  ),
  parameter_dir = here::here(
    "Data",
    "par"
  )
)

required_parameter_files <- c(
  "mtB.csv",
  "muB.csv",
  "rtB.csv",
  "mtS.csv",
  "muS.csv",
  "rtS.csv",
  "mtH.csv",
  "muH.csv",
  "rtH.csv"
)

input_set_is_complete <- function(input_set) {

  required_files <- c(
    input_set$reaches_mean,
    input_set$reaches_sd,
    input_set$reaches_metadata,
    input_set$lme_city3,
    file.path(
      input_set$parameter_dir,
      required_parameter_files
    )
  )

  all(
    file.exists(
      required_files
    )
  )
}

missing_input_files <- function(input_set) {

  required_files <- c(
    input_set$reaches_mean,
    input_set$reaches_sd,
    input_set$reaches_metadata,
    input_set$lme_city3,
    file.path(
      input_set$parameter_dir,
      required_parameter_files
    )
  )

  required_files[
    !file.exists(
      required_files
    )
  ]
}

select_input_set <- function(
    input_mode,
    generated_inputs,
    precomputed_inputs) {

  generated_complete <- input_set_is_complete(
    generated_inputs
  )

  precomputed_complete <- input_set_is_complete(
    precomputed_inputs
  )

  if (input_mode == "generated") {
    if (!generated_complete) {
      stop(
        "The generated Figure 9(d) input set is incomplete. ",
        "Missing:\n  ",
        paste(
          missing_input_files(
            generated_inputs
          ),
          collapse = "\n  "
        )
      )
    }

    return(
      generated_inputs
    )
  }

  if (input_mode == "precomputed") {
    if (!precomputed_complete) {
      stop(
        "The precomputed Figure 9(d) input set is incomplete. ",
        "Missing:\n  ",
        paste(
          missing_input_files(
            precomputed_inputs
          ),
          collapse = "\n  "
        )
      )
    }

    return(
      precomputed_inputs
    )
  }

  if (generated_complete) {
    return(
      generated_inputs
    )
  }

  if (precomputed_complete) {
    return(
      precomputed_inputs
    )
  }

  stop(
    "Neither a complete generated nor a complete precomputed ",
    "Figure 9(d) input set was found."
  )
}

input_files <- select_input_set(
  input_mode = input_mode,
  generated_inputs = generated_inputs,
  precomputed_inputs = precomputed_inputs
)

message(
  "Figure 9(d) input mode: ",
  input_mode
)

message(
  "REACHES input selected from: ",
  dirname(
    input_files$reaches_mean
  )
)

message(
  "Prior parameters selected from: ",
  input_files$parameter_dir
)

# ------------------------------------------------------------
# 3. General input helpers
# ------------------------------------------------------------

coordinate_key <- function(
    long,
    lat) {

  sprintf(
    "%.8f_%.8f",
    as.numeric(
      long
    ),
    as.numeric(
      lat
    )
  )
}

extract_year_columns <- function(data) {

  column_names <- names(
    data
  )

  year_values <- suppressWarnings(
    as.integer(
      sub(
        "^[Xx]",
        "",
        column_names
      )
    )
  )

  valid <- !is.na(
    year_values
  ) &
    year_values >= 1000L &
    year_values <= 3000L

  if (!any(
    valid
  )) {
    stop(
      "No annual columns were found. Expected names such as ",
      "'1368', 'X1368', or 'x1368'."
    )
  }

  output <- column_names[valid]

  names(
    output
  ) <- as.character(
    year_values[valid]
  )

  output
}

normalize_city_name <- function(city) {

  normalized <- tolower(
    gsub(
      "[[:space:]_-]",
      "",
      as.character(
        city
      )
    )
  )

  dplyr::recode(
    normalized,
    hongkong = "HongKong",
    shanghai = "Shanghai",
    beijing = "Beijing",
    .default = NA_character_
  )
}

# ------------------------------------------------------------
# 4. Read REACHES means, SDs, and process variance
# ------------------------------------------------------------

read_kriging_city_file <- function(
    input_file,
    value_type) {

  data <- readr::read_csv(
    input_file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  if (
    !"lat" %in%
      names(
        data
      ) &&
      "lati" %in%
        names(
          data
        )
  ) {
    data <- data %>%
      rename(
        lat = lati
      )
  }

  required_columns <- c(
    "long",
    "lat"
  )

  missing_columns <- setdiff(
    required_columns,
    names(
      data
    )
  )

  if (length(
    missing_columns
  ) > 0L) {
    stop(
      input_file,
      " is missing coordinate columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  year_map <- extract_year_columns(
    data
  )

  data <- data %>%
    mutate(
      long = as.numeric(
        long
      ),
      lat = as.numeric(
        lat
      )
    )

  if (
    any(!is.finite(
      data$long
    )) ||
      any(!is.finite(
        data$lat
      )) ||
      anyDuplicated(
        coordinate_key(
          data$long,
          data$lat
        )
      )
  ) {
    stop(
      "Invalid or duplicated coordinates were found in ",
      input_file,
      "."
    )
  }

  list(
    data = data,
    year_map = year_map,
    value_type = value_type
  )
}

kriging_mean_object <- read_kriging_city_file(
  input_files$reaches_mean,
  value_type = "mean"
)

kriging_sd_object <- read_kriging_city_file(
  input_files$reaches_sd,
  value_type = "standard deviation"
)

common_reaches_years <- intersect(
  names(
    kriging_mean_object$year_map
  ),
  names(
    kriging_sd_object$year_map
  )
)

common_reaches_years <- sort(
  as.integer(
    common_reaches_years
  )
)

common_reaches_years <- intersect(
  common_reaches_years,
  analysis_years
)

if (length(
  common_reaches_years
) < 2L) {
  stop(
    "The REACHES mean and SD files do not share at least two ",
    "analysis years."
  )
}

kriging_metadata <- readr::read_csv(
  input_files$reaches_metadata,
  show_col_types = FALSE
)

if (
  !"process_variance" %in%
    names(
      kriging_metadata
    ) ||
    nrow(
      kriging_metadata
    ) != 1L
) {
  stop(
    "The REACHES metadata file must contain one ",
    "'process_variance' value."
  )
}

sigmaY2 <- as.numeric(kriging_metadata$process_variance[1])

if (
  length(
    sigmaY2
  ) != 1L ||
    !is.finite(
      sigmaY2
    ) ||
    sigmaY2 <= 0
) {
  stop(
    "The process variance in the REACHES metadata is invalid."
  )
}

message(
  "Latent REACHES process variance: ",
  signif(
    sigmaY2,
    8
  )
)

extract_city_kriging <- function(
    city_name,
    long,
    lat,
    mean_object,
    sd_object,
    years) {

  target_key <- coordinate_key(
    long,
    lat
  )

  mean_keys <- coordinate_key(
    mean_object$data$long,
    mean_object$data$lat
  )

  sd_keys <- coordinate_key(
    sd_object$data$long,
    sd_object$data$lat
  )

  mean_row <- match(
    target_key,
    mean_keys
  )

  sd_row <- match(
    target_key,
    sd_keys
  )

  if (
    is.na(
      mean_row
    ) ||
      is.na(
        sd_row
      )
  ) {
    stop(
      "The REACHES location for ",
      city_name,
      " was not found at (",
      long,
      ", ",
      lat,
      ")."
    )
  }

  mean_columns <- unname(
    mean_object$year_map[
      as.character(
        years
      )
    ]
  )

  sd_columns <- unname(
    sd_object$year_map[
      as.character(
        years
      )
    ]
  )

  if (
    anyNA(
      mean_columns
    ) ||
      anyNA(
        sd_columns
      )
  ) {
    stop(
      "The REACHES files could not be aligned to all requested ",
      "event years for ",
      city_name,
      "."
    )
  }

  yhat <- as.numeric(
    mean_object$data[
      mean_row,
      mean_columns,
      drop = TRUE
    ]
  )

  nu <- as.numeric(
    sd_object$data[
      sd_row,
      sd_columns,
      drop = TRUE
    ]
  )

  if (
    any(!is.finite(
      yhat
    )) ||
      any(!is.finite(
        nu
      )) ||
      any(
        nu < 0
      )
  ) {
    stop(
      "Invalid REACHES means or standard deviations were found ",
      "for ",
      city_name,
      "."
    )
  }

  data.frame(
    year = as.integer(
      years
    ),
    reaches_index_mean = yhat,
    reaches_index_sd = nu
  )
}

# ------------------------------------------------------------
# 5. Read the common long-format LME city file
# ------------------------------------------------------------

lme_city_data <- readr::read_csv(
  input_files$lme_city3,
  show_col_types = FALSE
)

required_lme_columns <- c(
  "city",
  "long",
  "lati",
  "member",
  "year",
  "temperature_kelvin"
)

missing_lme_columns <- setdiff(
  required_lme_columns,
  names(
    lme_city_data
  )
)

if (length(
  missing_lme_columns
) > 0L) {
  stop(
    "The LME city file is missing columns: ",
    paste(
      missing_lme_columns,
      collapse = ", "
    )
  )
}

lme_city_data <- lme_city_data %>%
  transmute(
    city = normalize_city_name(
      city
    ),
    long = as.numeric(
      long
    ),
    lat = as.numeric(
      lati
    ),
    member = as.character(
      member
    ),
    year = as.integer(
      year
    ),
    temperature_celsius =
      as.numeric(
        temperature_kelvin
      ) -
      273.15
  ) %>%
  filter(
    year %in%
      analysis_years
  ) %>%
  arrange(
    city,
    member,
    year
  )

if (
  anyNA(
    lme_city_data
  ) ||
    any(!is.finite(
      lme_city_data$temperature_celsius
    ))
) {
  stop(
    "The LME city file contains missing or invalid values."
  )
}

duplicate_lme_rows <- lme_city_data %>%
  count(
    city,
    member,
    year,
    name = "number_of_rows"
  ) %>%
  filter(
    number_of_rows != 1L
  )

if (nrow(
  duplicate_lme_rows
) > 0L) {
  stop(
    "Each LME city-member-year combination must occur exactly ",
    "once."
  )
}

read_lme_city <- function(
    city_name,
    config,
    lme_city_data) {

  city_data <- lme_city_data %>%
    filter(
      city == city_name
    )

  if (nrow(
    city_data
  ) == 0L) {
    stop(
      "No LME data were found for ",
      city_name,
      "."
    )
  }

  city_coordinates <- city_data %>%
    distinct(
      long,
      lat
    )

  if (nrow(
    city_coordinates
  ) != 1L) {
    stop(
      "The LME input contains multiple coordinates for ",
      city_name,
      "."
    )
  }

  coordinate_difference <- max(
    abs(
      city_coordinates$long -
        config$long
    ),
    abs(
      city_coordinates$lat -
        config$lat
    )
  )

  if (
    !is.finite(
      coordinate_difference
    ) ||
      coordinate_difference >
        coordinate_tolerance
  ) {
    stop(
      "The LME coordinate for ",
      city_name,
      " is (",
      city_coordinates$long,
      ", ",
      city_coordinates$lat,
      "), but the required common coordinate is (",
      config$long,
      ", ",
      config$lat,
      ")."
    )
  }

  member_names <- sort(
    unique(
      city_data$member
    )
  )

  if (length(
    member_names
  ) != 13L) {
    stop(
      city_name,
      " should contain 13 LME ensemble members, but contains ",
      length(
        member_names
      ),
      "."
    )
  }

  annual_mean <- city_data %>%
    group_by(
      year
    ) %>%
    summarise(
      lme_ensemble_mean_celsius = mean(
        temperature_celsius
      ),
      .groups = "drop"
    ) %>%
    arrange(
      year
    )

  if (!identical(
    annual_mean$year,
    as.integer(
      analysis_years
    )
  )) {
    stop(
      "The LME annual grid for ",
      city_name,
      " is incomplete or misaligned."
    )
  }

  list(
    sample = city_data$temperature_celsius,
    annual_mean = annual_mean,
    members = member_names,
    coordinates = city_coordinates
  )
}

# ------------------------------------------------------------
# 6. Read penalized annual prior parameters
# ------------------------------------------------------------

read_penalized_parameter <- function(
    filename,
    parameter_name,
    expected_years) {

  if (!file.exists(
    filename
  )) {
    stop(
      "Parameter file was not found: ",
      filename
    )
  }

  data <- readr::read_csv(
    filename,
    show_col_types = FALSE
  )

  required_columns <- c(
    "year",
    "coefficient",
    "value"
  )

  missing_columns <- setdiff(
    required_columns,
    names(
      data
    )
  )

  if (length(
    missing_columns
  ) > 0L) {
    stop(
      filename,
      " is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  penalized <- data %>%
    transmute(
      year = as.integer(
        year
      ),
      coefficient = as.character(
        coefficient
      ),
      value = as.numeric(
        value
      )
    ) %>%
    filter(
      coefficient ==
        "Penalized ML"
    )

  duplicate_rows <- penalized %>%
    count(
      year,
      name = "number_of_rows"
    ) %>%
    filter(
      number_of_rows != 1L
    )

  if (
    nrow(
      duplicate_rows
    ) > 0L ||
      !identical(
        sort(
          penalized$year
        ),
        as.integer(
          expected_years
        )
      )
  ) {
    stop(
      parameter_name,
      " in ",
      filename,
      " is not aligned to the expected annual grid."
    )
  }

  values <- penalized$value[
    match(
      expected_years,
      penalized$year
    )
  ]

  if (
    length(
      values
    ) !=
      length(
        expected_years
      ) ||
      any(!is.finite(
        values
      ))
  ) {
    stop(
      "Invalid ",
      parameter_name,
      " values were found in ",
      filename,
      "."
    )
  }

  values
}

read_city_prior <- function(
    city_code,
    parameter_dir) {

  M <- read_penalized_parameter(
    file.path(
      parameter_dir,
      paste0(
        "mt",
        city_code,
        ".csv"
      )
    ),
    parameter_name = "M",
    expected_years =
      analysis_years[
        -length(
          analysis_years
        )
      ]
  )

  mu <- read_penalized_parameter(
    file.path(
      parameter_dir,
      paste0(
        "mu",
        city_code,
        ".csv"
      )
    ),
    parameter_name = "mu",
    expected_years = analysis_years
  )

  r2 <- read_penalized_parameter(
    file.path(
      parameter_dir,
      paste0(
        "rt",
        city_code,
        ".csv"
      )
    ),
    parameter_name = "r2",
    expected_years = analysis_years
  )

  if (any(
    r2 <= 0
  )) {
    stop(
      "All prior innovation variances must be positive."
    )
  }

  list(
    M = M,
    mu = mu,
    r2 = r2
  )
}

# ------------------------------------------------------------
# 7. Quantile-mapping construction
# ------------------------------------------------------------

build_quantile_mapping <- function(
    lme_sample,
    yhat,
    nu,
    sigmaY2,
    bandwidth_seed,
    mapping_grid_size,
    lme_cdf_grid_size,
    latent_tail_probability) {

  lme_sample <- as.numeric(
    lme_sample
  )

  yhat <- as.numeric(
    yhat
  )

  nu <- as.numeric(
    nu
  )

  if (
    length(
      yhat
    ) !=
      length(
        nu
      ) ||
      any(!is.finite(
        lme_sample
      )) ||
      any(!is.finite(
        yhat
      )) ||
      any(!is.finite(
        nu
      )) ||
      any(
        nu < 0
      )
  ) {
    stop(
      "Invalid inputs were supplied to the quantile-mapping ",
      "construction."
    )
  }

  lme_sd <- stats::sd(
    lme_sample
  )

  if (
    !is.finite(
      lme_sd
    ) ||
      lme_sd <= 0
  ) {
    stop(
      "The LME calibration sample has zero or invalid ",
      "variation."
    )
  }

  set.seed(
    bandwidth_seed
  )

  fx_bandwidth <- np::npudistbw(
    dat = data.frame(
      x = lme_sample
    )
  )

  lme_grid <- seq(
    min(
      lme_sample
    ) -
      8 *
        lme_sd,
    max(
      lme_sample
    ) +
      8 *
        lme_sd,
    length.out = lme_cdf_grid_size
  )

  fx_grid <- as.numeric(
    fitted(
      np::npudist(
        bws = fx_bandwidth,
        edat = data.frame(
          x = lme_grid
        )
      )
    )
  )

  # Numerical monotonicity protection for the estimated CDF.
  fx_grid <- cummax(
    pmin(
      pmax(
        fx_grid,
        0
      ),
      1
    )
  )

  keep_unique_probability <- !duplicated(
    fx_grid
  )

  fx_probability_unique <- fx_grid[
    keep_unique_probability
  ]

  lme_grid_unique <- lme_grid[
    keep_unique_probability
  ]

  if (length(
    fx_probability_unique
  ) < 2L) {
    stop(
      "The estimated LME CDF could not be inverted."
    )
  }

  Fx_inv <- stats::approxfun(
    x = fx_probability_unique,
    y = lme_grid_unique,
    rule = 2,
    ties = "ordered"
  )

  nu_safe <- pmax(
    nu,
    1e-8
  )

  FY_hat <- function(y) {

    y <- as.numeric(
      y
    )

    standardized <- outer(
      y,
      yhat,
      "-"
    )

    standardized <- sweep(
      standardized,
      2,
      nu_safe,
      "/"
    )

    rowMeans(
      pnorm(
        standardized
      )
    )
  }

  latent_sd <- sqrt(
    sigmaY2
  )

  gaussian_limit <- qnorm(
    1 -
      latent_tail_probability /
        2
  ) *
    latent_sd

  proxy_limit_lower <- min(
    yhat -
      8 *
        nu_safe
  )

  proxy_limit_upper <- max(
    yhat +
      8 *
        nu_safe
  )

  latent_lower <- min(
    -gaussian_limit,
    proxy_limit_lower
  )

  latent_upper <- max(
    gaussian_limit,
    proxy_limit_upper
  )

  latent_grid <- seq(
    latent_lower,
    latent_upper,
    length.out = mapping_grid_size
  )

  mapping_probability <- FY_hat(
    latent_grid
  )

  mapped_grid <- Fx_inv(
    pmin(
      pmax(
        mapping_probability,
        1e-12
      ),
      1 -
        1e-12
    )
  )

  g <- stats::approxfun(
    x = latent_grid,
    y = mapped_grid,
    rule = 2,
    ties = "ordered"
  )

  corrected <- g(
    yhat
  )

  list(
    corrected = corrected,
    g = g,
    latent_range = range(
      latent_grid
    ),
    fx_bandwidth = as.numeric(
      fx_bandwidth$bw
    ),
    lme_range = range(
      lme_sample
    )
  )
}

# ------------------------------------------------------------
# 8. Measurement-equation Monte Carlo
# ------------------------------------------------------------

compute_measurement_parameters_mc <- function(
    sigmaY2,
    vhat,
    g,
    n_mc,
    seed) {

  vhat <- as.numeric(
    vhat
  )

  if (
    any(!is.finite(
      vhat
    )) ||
      any(
        vhat <
          -variance_tolerance
      ) ||
      any(
        vhat >
          sigmaY2 +
            variance_tolerance
      )
  ) {
    stop(
      "All Var(Yhat_t) values must lie in [0, sigmaY2]."
    )
  }

  vhat <- pmin(
    pmax(
      vhat,
      0
    ),
    sigmaY2
  )

  set.seed(
    seed
  )

  number_of_years <- length(
    vhat
  )

  alpha <- numeric(
    number_of_years
  )

  beta <- numeric(
    number_of_years
  )

  vdelta <- numeric(
    number_of_years
  )

  maximum_absolute_latent_draw <- 0

  for (
    t in seq_len(
      number_of_years
    )
  ) {

    current_vhat <- vhat[t]

    # This exact representation gives:
    #   Var(Yhat) = current_vhat
    #   Var(Y)    = sigmaY2
    #   Cov(Y, Yhat) = current_vhat
    #
    # No rounding or clipping is applied because Y and Yhat are
    # continuous latent Gaussian variables.
    Yhat <- sqrt(
      current_vhat
    ) *
      rnorm(
        n_mc
      )

    Y <- Yhat +
      sqrt(
        pmax(
          sigmaY2 -
            current_vhat,
          0
        )
      ) *
      rnorm(
        n_mc
      )

    maximum_absolute_latent_draw <- max(
      maximum_absolute_latent_draw,
      abs(
        Y
      ),
      abs(
        Yhat
      )
    )

    X <- g(
      Y
    )

    Xstar <- g(
      Yhat
    )

    variance_X <- stats::var(
      X
    )

    if (
      !is.finite(
        variance_X
      ) ||
        variance_X <= 0
    ) {
      stop(
        "A non-positive Monte Carlo variance was obtained at ",
        "measurement year index ",
        t,
        "."
      )
    }

    beta_t <- stats::cov(
      Xstar,
      X
    ) /
      variance_X

    alpha_t <- mean(
      Xstar
    ) -
      beta_t *
        mean(
          X
        )

    delta <- Xstar -
      alpha_t -
      beta_t *
        X

    alpha[
      t
    ] <- alpha_t

    beta[
      t
    ] <- beta_t

    vdelta[
      t
    ] <- pmax(
      stats::var(
        delta
      ),
      1e-10
    )
  }

  list(
    alpha = alpha,
    beta = beta,
    vdelta = vdelta,
    maximum_absolute_latent_draw =
      maximum_absolute_latent_draw
  )
}

# ------------------------------------------------------------
# 9. Full-annual Kalman filter and RTS smoother
# ------------------------------------------------------------

kalman_filter_smoother_annual <- function(
    mu,
    M,
    r2,
    Xstar,
    alpha,
    beta,
    vdelta) {

  number_of_years <- length(
    mu
  )

  if (
    length(
      M
    ) !=
      number_of_years -
        1L ||
      length(
        r2
      ) !=
        number_of_years ||
      length(
        Xstar
      ) !=
        number_of_years ||
      length(
        alpha
      ) !=
        number_of_years ||
      length(
        beta
      ) !=
        number_of_years ||
      length(
        vdelta
      ) !=
        number_of_years
  ) {
    stop(
      "The state and observation vectors are not aligned to the ",
      "same annual grid."
    )
  }

  has_observation <- is.finite(
    Xstar
  )

  if (
    any(
      has_observation &
        (
          !is.finite(
            alpha
          ) |
            !is.finite(
              beta
            ) |
            !is.finite(
              vdelta
            ) |
            vdelta <= 0
        )
    )
  ) {
    stop(
      "At least one observed year has invalid measurement ",
      "parameters."
    )
  }

  X_pred <- numeric(
    number_of_years
  )

  P_pred <- numeric(
    number_of_years
  )

  X_filt <- numeric(
    number_of_years
  )

  P_filt <- numeric(
    number_of_years
  )

  kalman_gain <- rep(
    NA_real_,
    number_of_years
  )

  effective_proxy <- rep(
    NA_real_,
    number_of_years
  )

  observation_weight <- rep(
    NA_real_,
    number_of_years
  )

  innovation <- rep(
    NA_real_,
    number_of_years
  )

  X_pred[
    1
  ] <- mu[
    1
  ]

  P_pred[
    1
  ] <- r2[
    1
  ]

  for (
    t in seq_len(
      number_of_years
    )
  ) {

    if (t >= 2L) {

      transition <- M[
        t -
          1L
      ]

      X_pred[
        t
      ] <- mu[
        t
      ] +
        transition *
          (
            X_filt[
              t -
                1L
            ] -
              mu[
                t -
                  1L
              ]
          )

      P_pred[
        t
      ] <- transition^2 *
        P_filt[
          t -
            1L
        ] +
        r2[
          t
        ]
    }

    if (!has_observation[
      t
    ]) {

      # Prediction-only step for a year without documentary
      # information.
      X_filt[
        t
      ] <- X_pred[
        t
      ]

      P_filt[
        t
      ] <- P_pred[
        t
      ]

      next
    }

    denominator <- beta[
      t
    ]^2 *
      P_pred[
        t
      ] +
      vdelta[
        t
      ]

    if (
      !is.finite(
        denominator
      ) ||
        denominator <= 0
    ) {
      stop(
        "A non-positive Kalman update denominator was obtained ",
        "at year index ",
        t,
        "."
      )
    }

    kalman_gain[
      t
    ] <- beta[
      t
    ] *
      P_pred[
        t
      ] /
      denominator

    innovation[
      t
    ] <- Xstar[
      t
    ] -
      alpha[
        t
      ] -
      beta[
        t
      ] *
        X_pred[
          t
        ]

    X_filt[
      t
    ] <- X_pred[
      t
    ] +
      kalman_gain[
        t
      ] *
        innovation[
          t
        ]

    P_filt[
      t
    ] <- (
      1 -
        beta[
          t
        ] *
          kalman_gain[
            t
          ]
    ) *
      P_pred[
        t
      ]

    P_filt[
      t
    ] <- pmax(
      P_filt[
        t
      ],
      0
    )

    if (
      abs(
        beta[
          t
        ]
      ) >
        beta_tolerance
    ) {
      effective_proxy[
        t
      ] <- (
        Xstar[
          t
        ] -
          alpha[
            t
          ]
      ) /
        beta[
          t
        ]

      observation_weight[
        t
      ] <- beta[
        t
      ]^2 *
        P_pred[
          t
        ] /
        denominator
    }
  }

  X_smooth <- numeric(
    number_of_years
  )

  P_smooth <- numeric(
    number_of_years
  )

  smoother_gain <- numeric(
    number_of_years -
      1L
  )

  X_smooth[
    number_of_years
  ] <- X_filt[
    number_of_years
  ]

  P_smooth[
    number_of_years
  ] <- P_filt[
    number_of_years
  ]

  if (number_of_years >= 2L) {

    for (
      t in seq(
        from = number_of_years -
          1L,
        to = 1L,
        by = -1L
      )
    ) {

      if (
        !is.finite(
          P_pred[
            t +
              1L
          ]
        ) ||
          P_pred[
            t +
              1L
          ] <= 0
      ) {
        stop(
          "A non-positive predicted variance was obtained before ",
          "RTS smoothing at year index ",
          t +
            1L,
          "."
        )
      }

      smoother_gain[
        t
      ] <- P_filt[
        t
      ] *
        M[
          t
        ] /
        P_pred[
          t +
            1L
        ]

      X_smooth[
        t
      ] <- X_filt[
        t
      ] +
        smoother_gain[
          t
        ] *
          (
            X_smooth[
              t +
                1L
            ] -
              X_pred[
                t +
                  1L
              ]
          )

      P_smooth[
        t
      ] <- P_filt[
        t
      ] +
        smoother_gain[
          t
        ]^2 *
          (
            P_smooth[
              t +
                1L
            ] -
              P_pred[
                t +
                  1L
              ]
          )
    }
  }

  P_smooth <- pmax(
    P_smooth,
    0
  )

  filtered_check <- rep(
    NA_real_,
    number_of_years
  )

  check_years <- which(
    has_observation &
      is.finite(
        effective_proxy
      ) &
      is.finite(
        observation_weight
      )
  )

  filtered_check[
    check_years
  ] <- (
    1 -
      observation_weight[
        check_years
      ]
  ) *
    X_pred[
      check_years
    ] +
    observation_weight[
      check_years
    ] *
    effective_proxy[
      check_years
    ]

  maximum_filter_identity_error <- max(
    abs(
      X_filt[
        check_years
      ] -
        filtered_check[
          check_years
        ]
    ),
    na.rm = TRUE
  )

  if (
    length(
      check_years
    ) > 0L &&
      (
        !is.finite(
          maximum_filter_identity_error
        ) ||
          maximum_filter_identity_error >
            identity_tolerance
      )
  ) {
    stop(
      "The filtered-mean weighted-average identity failed. ",
      "Maximum absolute error: ",
      maximum_filter_identity_error
    )
  }

  list(
    has_observation = has_observation,
    X_pred = X_pred,
    P_pred = P_pred,
    X_filt = X_filt,
    P_filt = P_filt,
    X_smooth = X_smooth,
    P_smooth = P_smooth,
    kalman_gain = kalman_gain,
    smoother_gain = smoother_gain,
    innovation = innovation,
    effective_proxy = effective_proxy,
    observation_weight = observation_weight,
    filtered_check = filtered_check,
    maximum_filter_identity_error =
      maximum_filter_identity_error
  )
}

# ------------------------------------------------------------
# 10. Accuracy and displacement diagnostics
# ------------------------------------------------------------

compute_pair_metrics <- function(
    city,
    reference_name,
    estimate_name,
    reference,
    estimate) {

  valid <- is.finite(
    reference
  ) &
    is.finite(
      estimate
    )

  reference <- reference[
    valid
  ]

  estimate <- estimate[
    valid
  ]

  if (length(
    reference
  ) < 2L) {
    return(
      data.frame(
        city = city,
        reference = reference_name,
        estimate = estimate_name,
        number_of_years = length(
          reference
        ),
        correlation = NA_real_,
        mean_bias_estimate_minus_reference =
          NA_real_,
        RMSE = NA_real_,
        anomaly_correlation = NA_real_,
        anomaly_RMSE = NA_real_
      )
    )
  }

  reference_anomaly <- reference -
    mean(
      reference
    )

  estimate_anomaly <- estimate -
    mean(
      estimate
    )

  data.frame(
    city = city,
    reference = reference_name,
    estimate = estimate_name,
    number_of_years = length(
      reference
    ),
    correlation = stats::cor(
      estimate,
      reference
    ),
    mean_bias_estimate_minus_reference = mean(
      estimate -
        reference
    ),
    RMSE = sqrt(
      mean(
        (
          estimate -
            reference
        )^2
      )
    ),
    anomaly_correlation = stats::cor(
      estimate_anomaly,
      reference_anomaly
    ),
    anomaly_RMSE = sqrt(
      mean(
        (
          estimate_anomaly -
            reference_anomaly
        )^2
      )
    )
  )
}

# ------------------------------------------------------------
# 11. Process one city
# ------------------------------------------------------------

process_city <- function(
    city_name,
    config,
    city_index) {

  message(
    "\nProcessing ",
    city_name,
    "..."
  )

  city_kriging <- extract_city_kriging(
    city_name = city_name,
    long = config$long,
    lat = config$lat,
    mean_object = kriging_mean_object,
    sd_object = kriging_sd_object,
    years = common_reaches_years
  )

  lme <- read_lme_city(
    city_name = city_name,
    config = config,
    lme_city_data = lme_city_data
  )

  prior <- read_city_prior(
    city_code = config$code,
    parameter_dir =
      input_files$parameter_dir
  )

  quantile_mapping <- build_quantile_mapping(
    lme_sample = lme$sample,
    yhat =
      city_kriging$reaches_index_mean,
    nu =
      city_kriging$reaches_index_sd,
    sigmaY2 = sigmaY2,
    bandwidth_seed =
      10L +
      city_index,
    mapping_grid_size =
      mapping_grid_size,
    lme_cdf_grid_size =
      lme_cdf_grid_size,
    latent_tail_probability =
      latent_tail_probability
  )

  city_kriging <- city_kriging %>%
    mutate(
      reaches_mspe =
        reaches_index_sd^2,
      calibrated_reaches_celsius =
        quantile_mapping$corrected
    )

  implied_vhat_raw <- sigmaY2 -
    city_kriging$reaches_mspe

  clipped_vhat <- pmin(
    pmax(
      implied_vhat_raw,
      0
    ),
    sigmaY2
  )

  number_vhat_clipped <- sum(
    abs(
      clipped_vhat -
        implied_vhat_raw
    ) >
      variance_tolerance
  )

  if (
    number_vhat_clipped >
      0L
  ) {
    warning(
      city_name,
      ": ",
      number_vhat_clipped,
      " implied Var(Yhat_t) values were outside [0, sigmaY2] ",
      "and were clipped to that interval."
    )
  }

  measurement_parameters <-
    compute_measurement_parameters_mc(
      sigmaY2 = sigmaY2,
      vhat = clipped_vhat,
      g = quantile_mapping$g,
      n_mc = measurement_mc_size,
      seed =
        measurement_mc_seed +
        city_index
    )

  annual_data <- data.frame(
    year = analysis_years
  ) %>%
    left_join(
      city_kriging,
      by = "year"
    ) %>%
    left_join(
      lme$annual_mean,
      by = "year"
    )

  event_indices <- match(
    city_kriging$year,
    annual_data$year
  )

  annual_data$implied_latent_kriging_variance <-
    NA_real_

  annual_data$alpha <- NA_real_
  annual_data$beta <- NA_real_
  annual_data$vdelta <- NA_real_

  annual_data$implied_latent_kriging_variance[
    event_indices
  ] <- clipped_vhat

  annual_data$alpha[
    event_indices
  ] <- measurement_parameters$alpha

  annual_data$beta[
    event_indices
  ] <- measurement_parameters$beta

  annual_data$vdelta[
    event_indices
  ] <- measurement_parameters$vdelta

  smoother_output <- kalman_filter_smoother_annual(
    mu = prior$mu,
    M = prior$M,
    r2 = prior$r2,
    Xstar =
      annual_data$calibrated_reaches_celsius,
    alpha =
      annual_data$alpha,
    beta =
      annual_data$beta,
    vdelta =
      annual_data$vdelta
  )

  result <- annual_data %>%
    mutate(
      has_reaches_observation =
        smoother_output$has_observation,
      prior_mean = prior$mu,
      prior_innovation_variance =
        prior$r2,
      transition_M = c(
        prior$M,
        NA_real_
      ),
      dynamic_prior_prediction =
        smoother_output$X_pred,
      dynamic_prior_prediction_variance =
        smoother_output$P_pred,
      filtered_temperature_celsius =
        smoother_output$X_filt,
      filtered_variance =
        smoother_output$P_filt,
      smoothed_temperature_celsius =
        smoother_output$X_smooth,
      smoothed_variance =
        smoother_output$P_smooth,
      kalman_gain =
        smoother_output$kalman_gain,
      effective_proxy_temperature_celsius =
        smoother_output$effective_proxy,
      observation_weight =
        smoother_output$observation_weight,
      innovation =
        smoother_output$innovation,
      filtered_identity_check =
        smoother_output$filtered_check
    ) %>%
    relocate(
      year,
      has_reaches_observation,
      reaches_index_mean,
      reaches_index_sd,
      reaches_mspe,
      calibrated_reaches_celsius,
      alpha,
      beta,
      vdelta,
      effective_proxy_temperature_celsius,
      observation_weight,
      dynamic_prior_prediction,
      filtered_temperature_celsius,
      smoothed_temperature_celsius,
      lme_ensemble_mean_celsius
    )

  result_file <- file.path(
    assimilation_output_dir,
    paste0(
      "assimilation_",
      city_name,
      ".csv"
    )
  )

  readr::write_csv(
    result,
    result_file
  )

  main_figure <- ggplot(
    result,
    aes(
      x = year,
      y =
        smoothed_temperature_celsius
    )
  ) +
    geom_ribbon(
      aes(
        ymin =
          smoothed_temperature_celsius -
          sqrt(
            smoothed_variance
          ),
        ymax =
          smoothed_temperature_celsius +
          sqrt(
            smoothed_variance
          )
      ),
      alpha = 0.5,
      fill = "grey30"
    ) +
    geom_line(
      colour = "red"
    ) +
    labs(
      x = "Year",
      y = expression(
        "Temperature (" *
          degree *
          "C)"
      )
    ) +
    theme(
      text = element_text(
        size = 11
      ),
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5
      )
    )

  ggsave(
    filename = config$figure_file,
    plot = main_figure,
    width = 6,
    height = 3,
    units = "in",
    dpi = 300
  )

  diagnostic_long <- result %>%
    dplyr::select(
      year,
      `Calibrated REACHES` =
        calibrated_reaches_celsius,
      `Effective proxy` =
        effective_proxy_temperature_celsius,
      `Dynamic prior prediction` =
        dynamic_prior_prediction,
      `Filtered estimate` =
        filtered_temperature_celsius,
      `Smoothed estimate` =
        smoothed_temperature_celsius,
      `LME ensemble mean` =
        lme_ensemble_mean_celsius
    ) %>%
    pivot_longer(
      cols = -year,
      names_to = "series",
      values_to = "temperature_celsius"
    )

  diagnostic_figure <- ggplot(
    diagnostic_long,
    aes(
      x = year,
      y = temperature_celsius,
      colour = series,
      linetype = series
    )
  ) +
    geom_line(
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    labs(
      x = "Year",
      y = expression(
        "Temperature (" *
          degree *
          "C)"
      ),
      colour = NULL,
      linetype = NULL
    ) +
    theme(
      text = element_text(
        size = 10
      ),
      legend.position = "bottom"
    )

  diagnostic_file <- file.path(
    assimilation_output_dir,
    paste0(
      "diagnostic_",
      city_name,
      ".png"
    )
  )
  
  if (save_diagnostic_plots) {
  ggsave(
    filename = diagnostic_file,
    plot = diagnostic_figure,
    width = 8,
    height = 5,
    units = "in",
    dpi = 300
  )
}
  
  overlap <- result %>%
    filter(
      has_reaches_observation
    )

  pair_metrics <- bind_rows(
    compute_pair_metrics(
      city = city_name,
      reference_name =
        "Calibrated REACHES",
      estimate_name =
        "Smoothed estimate",
      reference =
        overlap$calibrated_reaches_celsius,
      estimate =
        overlap$smoothed_temperature_celsius
    ),
    compute_pair_metrics(
      city = city_name,
      reference_name =
        "LME ensemble mean",
      estimate_name =
        "Smoothed estimate",
      reference =
        overlap$lme_ensemble_mean_celsius,
      estimate =
        overlap$smoothed_temperature_celsius
    ),
    compute_pair_metrics(
      city = city_name,
      reference_name =
        "LME ensemble mean",
      estimate_name =
        "Calibrated REACHES",
      reference =
        overlap$lme_ensemble_mean_celsius,
      estimate =
        overlap$calibrated_reaches_celsius
    )
  )

  displacement_summary <- data.frame(
    city = city_name,
    reference = "Displayed-input interval",
    estimate = "Smoothed estimate",
    number_of_years = nrow(
      overlap
    ),
    correlation = NA_real_,
    mean_bias_estimate_minus_reference =
      NA_real_,
    RMSE = NA_real_,
    anomaly_correlation = NA_real_,
    anomaly_RMSE = NA_real_,
    proportion_below_both_inputs = mean(
      overlap$smoothed_temperature_celsius <
        pmin(
          overlap$calibrated_reaches_celsius,
          overlap$lme_ensemble_mean_celsius
        )
    ),
    proportion_above_both_inputs = mean(
      overlap$smoothed_temperature_celsius >
        pmax(
          overlap$calibrated_reaches_celsius,
          overlap$lme_ensemble_mean_celsius
        )
    )
  )

  pair_metrics$proportion_below_both_inputs <-
    NA_real_

  pair_metrics$proportion_above_both_inputs <-
    NA_real_

  metrics <- bind_rows(
    pair_metrics,
    displacement_summary
  )

  metadata <- data.frame(
    city = city_name,
    longitude = config$long,
    latitude = config$lat,
    process_variance = sigmaY2,
    number_of_annual_years =
      length(
        analysis_years
      ),
    number_of_reaches_event_years =
      nrow(
        city_kriging
      ),
    number_of_prediction_only_years =
      sum(
        !result$has_reaches_observation
      ),
    number_of_vhat_values_clipped =
      number_vhat_clipped,
    measurement_mc_size =
      measurement_mc_size,
    mapping_latent_lower =
      quantile_mapping$latent_range[
        1
      ],
    mapping_latent_upper =
      quantile_mapping$latent_range[
        2
      ],
    maximum_absolute_mc_latent_draw =
      measurement_parameters$
        maximum_absolute_latent_draw,
    maximum_filter_identity_error =
      smoother_output$
        maximum_filter_identity_error,
    minimum_beta = min(
      result$beta,
      na.rm = TRUE
    ),
    median_beta = median(
      result$beta,
      na.rm = TRUE
    ),
    maximum_beta = max(
      result$beta,
      na.rm = TRUE
    ),
    minimum_observation_weight = min(
      result$observation_weight,
      na.rm = TRUE
    ),
    median_observation_weight = median(
      result$observation_weight,
      na.rm = TRUE
    ),
    maximum_observation_weight = max(
      result$observation_weight,
      na.rm = TRUE
    )
  )

  message(
    "Saved figure: ",
    config$figure_file
  )

  message(
    "Saved annual assimilation diagnostics: ",
    result_file
  )

  list(
    result = result,
    metrics = metrics,
    metadata = metadata,
    main_figure = main_figure,
    diagnostic_figure = diagnostic_figure
  )
}

# ------------------------------------------------------------
# 12. Run all cities and save combined diagnostics
# ------------------------------------------------------------

city_names <- names(
  city_config
)

figure9d_results <- vector(
  "list",
  length(
    city_names
  )
)

names(
  figure9d_results
) <- city_names

for (
  city_index in seq_along(
    city_names
  )
) {

  city_name <- city_names[city_index]

  figure9d_results[[city_name]] <- process_city(
    city_name = city_name,
    config = city_config[[city_name]],
    city_index = city_index
  )

  gc()
}

combined_metrics <- bind_rows(
  lapply(
    figure9d_results,
    function(result) {
      result$metrics
    }
  )
)

combined_metadata <- bind_rows(
  lapply(
    figure9d_results,
    function(result) {
      result$metadata
    }
  )
)

readr::write_csv(
  combined_metrics,
  file.path(
    assimilation_output_dir,
    "assimilation_metrics.csv"
  )
)

readr::write_csv(
  combined_metadata,
  file.path(
    assimilation_output_dir,
    "assimilation_metadata.csv"
  )
)

saveRDS(
  figure9d_results,
  file.path(
    assimilation_output_dir,
    "assimilation_results_all_cities.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    assimilation_output_dir,
    "Figure9d_sessionInfo.txt"
  )
)

message(
  "\nCompleted Figure 9(d), Figure S3(d), and Figure S4(d)."
)

message(
  "All assimilation diagnostics were saved to: ",
  assimilation_output_dir
)
