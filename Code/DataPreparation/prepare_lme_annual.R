here::i_am("Code/DataPreparation/prepare_lme_annual.R")

# ============================================================
# Convert the 13 raw monthly LME files to annual data and
# extract annual LME series for Hong Kong, Shanghai, and Beijing.
#
# Raw-file structure after read.csv(..., row.names = 1):
#   column 1: lati
#   column 2: long
#   columns 3:7202: monthly temperatures from
#                   1350.01 through 1949.12
#
# The raw temperatures are retained in Kelvin.
#
# Inputs:
#   Data/LME data/a1.csv.gz, ..., a13.csv.gz
#
# Outputs:
#   Output/Intermediate/LME/lme_annual_1368_1911.rds
#   Output/Intermediate/LME/lme_ensemble_mean_1368_1911.csv
#   Output/Intermediate/LME/lme_annual_diagnostics.csv
#   Output/Intermediate/LME/lme_city3_annual_1368_1911.rds
#   Output/Intermediate/LME/lme_city3_annual_1368_1911.csv
#   Output/Intermediate/LME/
#     lme_city3_ensemble_mean_1368_1911.csv
#   Output/Intermediate/LME/
#     lme_city3_interpolation_diagnostics.csv
# ============================================================

library(here)
library(readr)
library(dplyr)
library(mgcv)

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

raw_first_year <- 1350L
raw_last_year <- 1949L

target_first_year <- 1368L
target_last_year <- 1911L

raw_years <- raw_first_year:raw_last_year
target_years <- target_first_year:target_last_year

number_of_members <- 13L
expected_number_of_locations <- 266L
expected_number_of_months <- length(raw_years) * 12L

# The city coordinates are aligned with the 0.5-degree
# REACHES prediction grid used in the downstream analysis.
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
  lati = c(
    22.25,
    31.25,
    39.75
  )
)

# Use a full-rank thin-plate regression-spline basis and a
# smoothing parameter fixed near zero, so the native LME field
# is represented as an approximately interpolating point-support
# surface. Both choices are saved in the output metadata.
spline_basis_dimension <- expected_number_of_locations
fixed_smoothing_parameter <- 1e-8

output_directory <- here::here(
  "Output",
  "Intermediate",
  "LME"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

annual_rds_file <- file.path(
  output_directory,
  "lme_annual_1368_1911.rds"
)

ensemble_mean_file <- file.path(
  output_directory,
  "lme_ensemble_mean_1368_1911.csv"
)

diagnostic_file <- file.path(
  output_directory,
  "lme_annual_diagnostics.csv"
)

city_rds_file <- file.path(
  output_directory,
  "lme_city3_annual_1368_1911.rds"
)

city_long_file <- file.path(
  output_directory,
  "lme_city3_annual_1368_1911.csv"
)

city_ensemble_mean_file <- file.path(
  output_directory,
  "lme_city3_ensemble_mean_1368_1911.csv"
)

city_diagnostic_file <- file.path(
  output_directory,
  "lme_city3_interpolation_diagnostics.csv"
)

# ------------------------------------------------------------
# 2. Locate the 13 raw member files
# ------------------------------------------------------------

input_directory <- here::here(
  "Data",
  "LME data"
)

member_files <- file.path(
  input_directory,
  paste0(
    "a",
    seq_len(number_of_members),
    ".csv.gz"
  )
)

missing_files <- member_files[
  !file.exists(member_files)
]

if (length(missing_files) > 0L) {
  stop(
    paste0(
      "The following raw LME files were not found:\n",
      paste(
        paste0(
          "  - ",
          missing_files
        ),
        collapse = "\n"
      )
    )
  )
}

message(
  "Using raw LME files from: ",
  input_directory
)


read_lme_member <- function(
    input_file) {

  input_connection <- gzfile(
    input_file,
    open = "rt"
  )

  on.exit(
    close(input_connection),
    add = TRUE
  )

  read.csv(
    input_connection,
    row.names = 1,
    check.names = FALSE
  )
}

# ------------------------------------------------------------
# 3. Convert one member from monthly to annual
# ------------------------------------------------------------

monthly_to_annual <- function(
    raw_data,
    member_name) {

  if (nrow(raw_data) != expected_number_of_locations) {
    stop(
      member_name,
      " contains ",
      nrow(raw_data),
      " locations; expected ",
      expected_number_of_locations,
      "."
    )
  }

  if (ncol(raw_data) < 3L) {
    stop(
      member_name,
      " does not contain coordinate and monthly-data columns."
    )
  }

  coordinates <- raw_data[
    ,
    1:2,
    drop = FALSE
  ]

  names(coordinates) <- c(
    "lati",
    "long"
  )

  coordinates <- coordinates %>%
    mutate(
      lati = as.numeric(lati),
      long = as.numeric(long)
    )

  if (
    any(!is.finite(coordinates$lati)) ||
      any(!is.finite(coordinates$long))
  ) {
    stop(
      member_name,
      " contains non-finite coordinates."
    )
  }

  monthly_values <- data.matrix(
    raw_data[
      ,
      -c(1, 2),
      drop = FALSE
    ]
  )

  if (ncol(monthly_values) != expected_number_of_months) {
    stop(
      member_name,
      " contains ",
      ncol(monthly_values),
      " monthly columns; expected ",
      expected_number_of_months,
      " for 1350.01--1949.12."
    )
  }

  if (anyNA(monthly_values)) {
    stop(
      member_name,
      " contains missing monthly temperatures."
    )
  }

  annual_all <- matrix(
    NA_real_,
    nrow = nrow(monthly_values),
    ncol = length(raw_years),
    dimnames = list(
      NULL,
      as.character(raw_years)
    )
  )

  for (year_index in seq_along(raw_years)) {

    month_start <- (
      year_index -
        1L
    ) *
      12L +
      1L

    month_end <- month_start +
      11L

    annual_all[
      ,
      year_index
    ] <- rowMeans(
      monthly_values[
        ,
        month_start:month_end,
        drop = FALSE
      ]
    )
  }

  target_indices <- match(
    target_years,
    raw_years
  )

  annual_target <- annual_all[
    ,
    target_indices,
    drop = FALSE
  ]

  list(
    coordinates = coordinates,
    annual_kelvin = annual_target
  )
}

# ------------------------------------------------------------
# 4. Read and process all 13 members
# ------------------------------------------------------------

annual_array <- array(
  NA_real_,
  dim = c(
    expected_number_of_locations,
    length(target_years),
    number_of_members
  ),
  dimnames = list(
    location_id = as.character(
      seq_len(
        expected_number_of_locations
      )
    ),
    year = as.character(
      target_years
    ),
    member = paste0(
      "a",
      seq_len(
        number_of_members
      )
    )
  )
)

reference_coordinates <- NULL
diagnostic_rows <- vector(
  "list",
  number_of_members
)

for (member_id in seq_len(number_of_members)) {

  member_name <- paste0(
    "a",
    member_id
  )

  input_file <- member_files[
    member_id
  ]

  message(
    "Processing ",
    member_name,
    " (",
    member_id,
    "/",
    number_of_members,
    ")..."
  )

  raw_data <- read_lme_member(
    input_file
  )

  converted <- monthly_to_annual(
    raw_data = raw_data,
    member_name = member_name
  )

  if (is.null(reference_coordinates)) {

    reference_coordinates <- converted$coordinates

  } else {

    coordinates_match <- isTRUE(
      all.equal(
        reference_coordinates,
        converted$coordinates,
        tolerance = 1e-10,
        check.attributes = FALSE
      )
    )

    if (!coordinates_match) {
      stop(
        "The coordinate ordering in ",
        member_name,
        " does not match a1."
      )
    }
  }

  annual_array[
    ,
    ,
    member_id
  ] <- converted$annual_kelvin

  diagnostic_rows[[member_id]] <- data.frame(
    member = member_name,
    input_file = normalizePath(
      input_file,
      winslash = "/",
      mustWork = TRUE
    ),
    number_of_locations = nrow(
      converted$annual_kelvin
    ),
    number_of_years = ncol(
      converted$annual_kelvin
    ),
    first_year = min(
      target_years
    ),
    last_year = max(
      target_years
    ),
    minimum_kelvin = min(
      converted$annual_kelvin
    ),
    maximum_kelvin = max(
      converted$annual_kelvin
    ),
    mean_kelvin = mean(
      converted$annual_kelvin
    )
  )

  rm(
    raw_data,
    converted
  )

  invisible(
    gc()
  )
}

if (anyNA(annual_array)) {
  stop(
    "Missing values remain in the completed annual LME array."
  )
}

# ------------------------------------------------------------
# 5. Save the annual archive
# ------------------------------------------------------------

lme_annual_archive <- list(
  coordinates = data.frame(
    location_id = seq_len(
      expected_number_of_locations
    ),
    lati = reference_coordinates$lati,
    long = reference_coordinates$long
  ),
  years = target_years,
  members = paste0(
    "a",
    seq_len(
      number_of_members
    )
  ),
  units = "Kelvin",
  source_frequency = "monthly",
  annualization = "arithmetic mean of 12 monthly values",
  annual_kelvin = annual_array
)

saveRDS(
  lme_annual_archive,
  annual_rds_file,
  compress = "xz"
)

# ------------------------------------------------------------
# 6. Save the ensemble-mean annual field
# ------------------------------------------------------------

ensemble_mean_kelvin <- apply(
  annual_array,
  c(
    1,
    2
  ),
  mean
)

ensemble_mean_output <- data.frame(
  location_id = seq_len(
    expected_number_of_locations
  ),
  lati = reference_coordinates$lati,
  long = reference_coordinates$long,
  ensemble_mean_kelvin,
  check.names = FALSE
)

names(
  ensemble_mean_output
)[4:ncol(ensemble_mean_output)] <- paste0(
  "x",
  target_years
)

readr::write_csv(
  ensemble_mean_output,
  ensemble_mean_file
)

# ------------------------------------------------------------
# 7. Save annual-data diagnostics
# ------------------------------------------------------------

diagnostics <- bind_rows(
  diagnostic_rows
)

readr::write_csv(
  diagnostics,
  diagnostic_file
)

# ------------------------------------------------------------
# 8. Build one reusable thin-plate interpolation operator
# ------------------------------------------------------------

coordinate_data <- data.frame(
  long = reference_coordinates$long,
  lati = reference_coordinates$lati
)

if (
  nrow(
    unique(
      coordinate_data
    )
  ) != expected_number_of_locations
) {
  stop(
    "The native LME coordinates are not all unique."
  )
}

smooth_object <- mgcv::smoothCon(
  object = mgcv::s(
    long,
    lati,
    bs = "tp",
    k = spline_basis_dimension
  ),
  data = coordinate_data,
  absorb.cons = TRUE,
  scale.penalty = TRUE
)[[1]]

if (length(smooth_object$S) != 1L) {
  stop(
    "The thin-plate smooth unexpectedly produced more than ",
    "one penalty matrix."
  )
}

training_design <- cbind(
  intercept = 1,
  smooth_object$X
)

city_prediction_design <- cbind(
  intercept = 1,
  mgcv::PredictMat(
    smooth_object,
    city_locations[
      ,
      c(
        "long",
        "lati"
      )
    ]
  )
)

penalty_matrix <- matrix(
  0,
  nrow = ncol(
    training_design
  ),
  ncol = ncol(
    training_design
  )
)

penalty_matrix[
  -1,
  -1
] <- smooth_object$S[[1]]

normal_equation_matrix <- crossprod(
  training_design
) +
  fixed_smoothing_parameter *
    penalty_matrix

reciprocal_condition_number <- rcond(
  normal_equation_matrix
)

if (
  !is.finite(
    reciprocal_condition_number
  ) ||
    reciprocal_condition_number <
      1e-14
) {
  stop(
    "The thin-plate interpolation system is numerically ",
    "singular. Increase fixed_smoothing_parameter slightly."
  )
}

# With fixed coordinates, basis dimension, and smoothing
# parameter, the fitted values at the three cities are a fixed
# linear transformation of every native-grid annual field.
# Construct that transformation once instead of fitting 7,072
# separate GAMs.
city_interpolation_operator <-
  city_prediction_design %*%
  solve(
    normal_equation_matrix,
    t(
      training_design
    )
  )

# ------------------------------------------------------------
# 9. Apply the interpolation operator to all years and members
# ------------------------------------------------------------

annual_field_matrix <- matrix(
  annual_array,
  nrow = expected_number_of_locations,
  ncol = length(
    target_years
  ) *
    number_of_members
)

city_prediction_matrix <-
  city_interpolation_operator %*%
  annual_field_matrix

city_annual_array <- array(
  city_prediction_matrix,
  dim = c(
    nrow(
      city_locations
    ),
    length(
      target_years
    ),
    number_of_members
  ),
  dimnames = list(
    city = city_locations$city,
    year = as.character(
      target_years
    ),
    member = paste0(
      "a",
      seq_len(
        number_of_members
      )
    )
  )
)

if (anyNA(city_annual_array)) {
  stop(
    "Missing values remain in the three-city LME array."
  )
}

# ------------------------------------------------------------
# 10. Validate the fast operator against one direct mgcv fit
# ------------------------------------------------------------

validation_data <- coordinate_data
validation_data$temperature_kelvin <-
  annual_array[
    ,
    1,
    1
  ]

direct_validation_fit <- mgcv::gam(
  temperature_kelvin ~
    s(
      long,
      lati,
      bs = "tp",
      k = spline_basis_dimension,
      sp = fixed_smoothing_parameter
    ),
  data = validation_data,
  method = "REML"
)

direct_validation_prediction <- as.numeric(
  predict(
    direct_validation_fit,
    newdata = city_locations[
      ,
      c(
        "long",
        "lati"
      )
    ]
  )
)

operator_validation_prediction <-
  as.numeric(
    city_annual_array[
      ,
      1,
      1
    ]
  )

maximum_validation_difference <- max(
  abs(
    direct_validation_prediction -
      operator_validation_prediction
  )
)

if (
  !is.finite(
    maximum_validation_difference
  ) ||
    maximum_validation_difference >
      1e-5
) {
  stop(
    "The reusable interpolation operator did not reproduce ",
    "the direct mgcv fit. Maximum difference = ",
    signif(
      maximum_validation_difference,
      6
    ),
    "."
  )
}

# ------------------------------------------------------------
# 11. Save the three-city annual products
# ------------------------------------------------------------

city_archive <- list(
  coordinates = city_locations,
  years = target_years,
  members = paste0(
    "a",
    seq_len(
      number_of_members
    )
  ),
  units = "Kelvin",
  spatial_method =
    "thin-plate regression spline evaluated at city grid centres",
  basis_dimension =
    spline_basis_dimension,
  fixed_smoothing_parameter =
    fixed_smoothing_parameter,
  annual_kelvin =
    city_annual_array
)

saveRDS(
  city_archive,
  city_rds_file,
  compress = "xz"
)

city_long_output <- expand.grid(
  city = city_locations$city,
  year = target_years,
  member = paste0(
    "a",
    seq_len(
      number_of_members
    )
  ),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

city_long_output <- city_long_output %>%
  mutate(
    long = city_locations$long[
      match(
        city,
        city_locations$city
      )
    ],
    lati = city_locations$lati[
      match(
        city,
        city_locations$city
      )
    ],
    temperature_kelvin =
      as.vector(
        city_annual_array
      )
  ) %>%
  dplyr::select(
    city,
    long,
    lati,
    member,
    year,
    temperature_kelvin
  )

readr::write_csv(
  city_long_output,
  city_long_file
)

city_ensemble_mean <- apply(
  city_annual_array,
  c(
    1,
    2
  ),
  mean
)

city_ensemble_mean_output <- expand.grid(
  city = city_locations$city,
  year = target_years,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    long = city_locations$long[
      match(
        city,
        city_locations$city
      )
    ],
    lati = city_locations$lati[
      match(
        city,
        city_locations$city
      )
    ],
    ensemble_mean_kelvin =
      as.vector(
        city_ensemble_mean
      )
  ) %>%
  dplyr::select(
    city,
    long,
    lati,
    year,
    ensemble_mean_kelvin
  )

readr::write_csv(
  city_ensemble_mean_output,
  city_ensemble_mean_file
)

city_interpolation_diagnostics <- data.frame(
  basis_dimension =
    spline_basis_dimension,
  fixed_smoothing_parameter =
    fixed_smoothing_parameter,
  reciprocal_condition_number =
    reciprocal_condition_number,
  validation_member = "a1",
  validation_year =
    target_years[1],
  maximum_validation_difference =
    maximum_validation_difference,
  minimum_city_kelvin = min(
    city_annual_array
  ),
  maximum_city_kelvin = max(
    city_annual_array
  )
)

readr::write_csv(
  city_interpolation_diagnostics,
  city_diagnostic_file
)

# ------------------------------------------------------------
# 12. Completion messages
# ------------------------------------------------------------

message(
  "Saved annual LME archive: ",
  annual_rds_file
)

message(
  "Array dimensions: ",
  paste(
    dim(
      annual_array
    ),
    collapse = " x "
  ),
  " (locations x years x members)"
)

message(
  "Saved ensemble mean: ",
  ensemble_mean_file
)

message(
  "Saved annual diagnostics: ",
  diagnostic_file
)

message(
  "Saved three-city LME archive: ",
  city_rds_file
)

message(
  "Saved three-city long-format data: ",
  city_long_file
)

message(
  "Saved three-city ensemble means: ",
  city_ensemble_mean_file
)

message(
  "Saved interpolation diagnostics: ",
  city_diagnostic_file
)
