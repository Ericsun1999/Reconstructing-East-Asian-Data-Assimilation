here::i_am("Code/DataPreparation/prepare_lme_annual.R")

# ============================================================
# Convert the 13 raw monthly LME files to annual data.
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
#   Data/LME data/population/a1.csv.gz, ..., a13.csv.gz
#
# Also accepted:
#   Data/LME dara/population/
#   uncompressed a1.csv, ..., a13.csv
#
# Outputs:
#   Output/Intermediate/LME/lme_annual_1368_1911.rds
#   Output/Intermediate/LME/lme_ensemble_mean_1368_1911.csv.gz
#   Output/Intermediate/LME/lme_annual_diagnostics.csv
# ============================================================

library(here)
library(readr)
library(dplyr)

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

candidate_input_directories <- c(
  here::here("Data", "LME data", "population"),
  here::here("Data", "LME dara", "population")
)

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
  "lme_ensemble_mean_1368_1911.csv.gz"
)

diagnostic_file <- file.path(
  output_directory,
  "lme_annual_diagnostics.csv"
)

# ------------------------------------------------------------
# 2. File helpers
# ------------------------------------------------------------

find_member_file <- function(
    directory,
    member_id) {

  candidates <- c(
    file.path(
      directory,
      paste0(
        "a",
        member_id,
        ".csv.gz"
      )
    ),
    file.path(
      directory,
      paste0(
        "a",
        member_id,
        ".csv"
      )
    )
  )

  existing <- candidates[
    file.exists(candidates)
  ]

  if (length(existing) == 0L) {
    return(NA_character_)
  }

  # Prefer the compressed file when both exist.
  existing[1]
}


valid_input_directory <- vapply(
  candidate_input_directories,
  function(directory) {

    member_files <- vapply(
      seq_len(number_of_members),
      function(member_id) {
        find_member_file(
          directory,
          member_id
        )
      },
      character(1)
    )

    all(
      !is.na(member_files)
    )
  },
  logical(1)
)

if (!any(valid_input_directory)) {
  stop(
    paste0(
      "Could not find all 13 LME member files in one directory.\n",
      "Each member may be stored as a*.csv.gz or a*.csv.\n",
      "Directories searched:\n",
      paste(
        paste0(
          "  - ",
          candidate_input_directories
        ),
        collapse = "\n"
      )
    )
  )
}

input_directory <- candidate_input_directories[
  which(valid_input_directory)[1]
]

message(
  "Using raw LME files from: ",
  input_directory
)


read_lme_member <- function(
    input_file) {

  input_connection <- if (
    grepl(
      "\\.gz$",
      input_file,
      ignore.case = TRUE
    )
  ) {
    gzfile(
      input_file,
      open = "rt"
    )
  } else {
    file(
      input_file,
      open = "rt"
    )
  }

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

  input_file <- find_member_file(
    input_directory,
    member_id
  )

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

  diagnostic_rows[
    [member_id]
  ] <- data.frame(
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
)[
  4:ncol(
    ensemble_mean_output
  )
] <- paste0(
  "x",
  target_years
)

readr::write_csv(
  ensemble_mean_output,
  ensemble_mean_file
)

# ------------------------------------------------------------
# 7. Save diagnostics
# ------------------------------------------------------------

diagnostics <- bind_rows(
  diagnostic_rows
)

readr::write_csv(
  diagnostics,
  diagnostic_file
)

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
  "Saved diagnostics: ",
  diagnostic_file
)
