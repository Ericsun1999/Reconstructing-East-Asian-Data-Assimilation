here::i_am("Code/Analysis/coverage_population_GAM.R")

# ============================================================
# Documentary coverage and historical-population GAM analysis
# corresponding to manuscript Section 5.1 and Table 2.
#
# Analysis grid:
#   Native LME land grid over East Asia
#
# Analysis years:
#   Every calendar year from 1368 through 1911, including
#   years without REACHES temperature records.
#
# Response:
#   R_gu = 1 when at least one distinct REACHES site in LME
#   grid cell g has a temperature record in year u; 0 otherwise.
#
# Climate summaries:
#   Lbar_g, sigma_g, pcold_g, and pwarm_g
#
# Historical-population covariate:
#   H_gu = log1p(pop_hist)
#
# Models:
#   M1: geography + year
#   M2: geography + year + climate
#   M3: geography + year + population
#   M4: geography + year + climate + population
#
# Expected inputs:
#   Data/temperature index value.v1.xlsx
#
#   Generated LME annual archive:
#     Output/Intermediate/LME/lme_annual_1368_1911.rds
#
#   Precomputed fallback:
#     Data/LME data/precomputed/lme_annual_1368_1911.rds
#
#   Data/population/1776_pd.tif
#   Data/population/1820_pd.tif
#   Data/population/1851_pd.tif
#   Data/population/1880_pd.tif
#   Data/population/1910_pd.tif
#   Data/population/1953_pd.tif
#
# Main outputs:
#   Output/Tables/Table2_coverage_GAM.csv
#   Output/Tables/Table2_coverage_GAM_formatted.csv
#   Output/Tables/Table2_nested_comparisons.csv
#   Output/Tables/Table2_Hgu_coefficient.csv
#
# Intermediate outputs:
#   Output/Intermediate/coverage_population/
#
# Important implementation details:
#   - The already annualized LME archive is reused; the 13 raw
#     monthly files are not reread or reannualized here.
#   - LME temperature means are converted from Kelvin to Celsius.
#   - Only finite REACHES temperature-index records contribute
#     to documentary coverage.
#   - All four GAMs use the same complete-case panel.
#   - ML, rather than REML, is used because models with different
#     parametric terms are compared by AIC, BIC, and nested tests.
# ============================================================

library(here)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(terra)
library(sf)
library(mgcv)

# ------------------------------------------------------------
# 1. Paths and analysis settings
# ------------------------------------------------------------

analysis_years <- 1368:1911
number_of_lme_members <- 13L

reaches_file <- here::here(
  "Data",
  "temperature index value.v1.xlsx"
)

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

generated_lme_archive_file <- here::here(
  "Output",
  "Intermediate",
  "LME",
  "lme_annual_1368_1911.rds"
)

precomputed_lme_archive_file <- here::here(
  "Data",
  "LME data",
  "precomputed",
  "lme_annual_1368_1911.rds"
)

select_lme_archive_file <- function(
    input_mode,
    generated_file,
    precomputed_file) {

  generated_available <- file.exists(
    generated_file
  )

  precomputed_available <- file.exists(
    precomputed_file
  )

  if (input_mode == "generated") {
    if (!generated_available) {
      stop(
        "Generated LME annual archive was not found: ",
        generated_file
      )
    }

    return(
      generated_file
    )
  }

  if (input_mode == "precomputed") {
    if (!precomputed_available) {
      stop(
        "Precomputed LME annual archive was not found: ",
        precomputed_file
      )
    }

    return(
      precomputed_file
    )
  }

  if (generated_available) {
    return(
      generated_file
    )
  }

  if (precomputed_available) {
    return(
      precomputed_file
    )
  }

  stop(
    "Neither the generated nor precomputed LME annual archive ",
    "was found.\nGenerated path:\n  ",
    generated_file,
    "\nPrecomputed path:\n  ",
    precomputed_file
  )
}

lme_archive_file <- select_lme_archive_file(
  input_mode = input_mode,
  generated_file = generated_lme_archive_file,
  precomputed_file = precomputed_lme_archive_file
)

message(
  "Using annual LME archive: ",
  lme_archive_file
)

# ML is used because M1--M4 differ in their parametric terms and
# are compared using AIC, BIC, and nested likelihood-based tests.
gam_fit_method <- "ML"

# Number of locations processed per nearest-grid calculation.
nearest_cell_chunk_size <- 10000L

population_years_available <- c(
  1776L,
  1820L,
  1851L,
  1880L,
  1910L,
  1953L
)

required_population_files <- paste0(
  population_years_available,
  "_pd.tif"
)

valid_population_dir <- vapply(
  candidate_population_dirs,
  function(directory) {
    all(
      file.exists(
        file.path(
          directory,
          required_population_files
        )
      )
    )
  },
  logical(1)
)

if (!any(valid_population_dir)) {
  stop(
    "Could not find all six population rasters together in:\n",
    paste(
      paste0("  - ", candidate_population_dirs),
      collapse = "\n"
    ),
    "\nRequired files:\n",
    paste(
      paste0("  - ", required_population_files),
      collapse = "\n"
    )
  )
}

population_dir <- candidate_population_dirs[which(valid_population_dir)[1]]

if (!file.exists(reaches_file)) {
  stop(
    "The REACHES input file was not found: ",
    reaches_file
  )
}

# ------------------------------------------------------------
# 2. General helper functions
# ------------------------------------------------------------

mean_or_na <- function(x) {

  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(
    x,
    na.rm = TRUE
  )
}


nearest_cell_ids <- function(
    long,
    lat,
    grid_df,
    chunk_size = nearest_cell_chunk_size) {

  long <- as.numeric(
    long
  )

  lat <- as.numeric(
    lat
  )

  if (length(long) != length(lat)) {
    stop(
      "The longitude and latitude vectors have different lengths."
    )
  }

  if (
    any(!is.finite(
      long
    )) ||
      any(!is.finite(
        lat
      ))
  ) {
    stop(
      "Non-finite coordinates were supplied to nearest_cell_ids()."
    )
  }

  required_grid_columns <- c(
    "cell_id",
    "long",
    "lat"
  )

  if (!all(
    required_grid_columns %in%
      names(
        grid_df
      )
  )) {
    stop(
      "grid_df must contain: ",
      paste(
        required_grid_columns,
        collapse = ", "
      ),
      "."
    )
  }

  if (
    nrow(
      grid_df
    ) == 0L ||
      any(!is.finite(
        grid_df$long
      )) ||
      any(!is.finite(
        grid_df$lat
      ))
  ) {
    stop(
      "The LME grid contains invalid coordinates."
    )
  }

  number_of_points <- length(
    long
  )

  nearest_ids <- integer(
    number_of_points
  )

  chunk_starts <- seq.int(
    from = 1L,
    to = number_of_points,
    by = chunk_size
  )

  for (chunk_start in chunk_starts) {

    chunk_end <- min(
      chunk_start +
        chunk_size -
        1L,
      number_of_points
    )

    chunk_indices <- chunk_start:chunk_end

    longitude_difference <- outer(
      long[
        chunk_indices
      ],
      grid_df$long,
      "-"
    )

    latitude_difference <- outer(
      lat[
        chunk_indices
      ],
      grid_df$lat,
      "-"
    )

    squared_distance <-
      longitude_difference^2 +
      latitude_difference^2

    nearest_columns <- max.col(
      -squared_distance,
      ties.method = "first"
    )

    nearest_ids[
      chunk_indices
    ] <- grid_df$cell_id[
      nearest_columns
    ]
  }

  nearest_ids
}


first_existing_column <- function(
    data,
    patterns) {

  for (pattern in patterns) {

    matched <- grep(
      pattern,
      names(data),
      value = TRUE
    )

    if (length(matched) > 0L) {
      return(matched[1])
    }
  }

  NA_character_
}

# ------------------------------------------------------------
# 3. Read the prepared annual LME archive
# ------------------------------------------------------------

lme_archive <- readRDS(
  lme_archive_file
)

required_archive_components <- c(
  "coordinates",
  "years",
  "members",
  "units",
  "annual_kelvin"
)

missing_archive_components <- setdiff(
  required_archive_components,
  names(
    lme_archive
  )
)

if (length(
  missing_archive_components
) > 0L) {
  stop(
    "The LME annual archive is missing components: ",
    paste(
      missing_archive_components,
      collapse = ", "
    ),
    "."
  )
}

lme_coordinates <- lme_archive$coordinates

required_coordinate_columns <- c(
  "location_id",
  "lati",
  "long"
)

missing_coordinate_columns <- setdiff(
  required_coordinate_columns,
  names(
    lme_coordinates
  )
)

if (length(
  missing_coordinate_columns
) > 0L) {
  stop(
    "The LME coordinate table is missing: ",
    paste(
      missing_coordinate_columns,
      collapse = ", "
    ),
    "."
  )
}

lme_coordinates <- lme_coordinates %>%
  transmute(
    location_id = as.integer(
      location_id
    ),
    lat = as.numeric(
      lati
    ),
    long = as.numeric(
      long
    )
  ) %>%
  arrange(
    location_id
  )

if (
  anyNA(
    lme_coordinates
  ) ||
    any(!is.finite(
      lme_coordinates$lat
    )) ||
    any(!is.finite(
      lme_coordinates$long
    )) ||
    anyDuplicated(
      lme_coordinates$location_id
    ) ||
    anyDuplicated(
      paste(
        sprintf(
          "%.8f",
          lme_coordinates$long
        ),
        sprintf(
          "%.8f",
          lme_coordinates$lat
        ),
        sep = "_"
      )
    )
) {
  stop(
    "The LME annual archive contains invalid or duplicated ",
    "location metadata."
  )
}

archive_years <- as.integer(
  lme_archive$years
)

archive_members <- as.character(
  lme_archive$members
)

annual_kelvin <- lme_archive$annual_kelvin

if (!all(
  analysis_years %in%
    archive_years
)) {
  stop(
    "The LME annual archive does not contain the complete ",
    min(
      analysis_years
    ),
    "--",
    max(
      analysis_years
    ),
    " analysis period."
  )
}

if (length(
  archive_members
) != number_of_lme_members) {
  stop(
    "The LME annual archive contains ",
    length(
      archive_members
    ),
    " members; expected ",
    number_of_lme_members,
    "."
  )
}

if (
  length(
    dim(
      annual_kelvin
    )
  ) != 3L ||
    dim(
      annual_kelvin
    )[
      1
    ] !=
      nrow(
        lme_coordinates
      ) ||
    dim(
      annual_kelvin
    )[
      2
    ] !=
      length(
        archive_years
      ) ||
    dim(
      annual_kelvin
    )[
      3
    ] !=
      length(
        archive_members
      )
) {
  stop(
    "The dimensions of annual_kelvin do not agree with the ",
    "coordinate, year, and member metadata."
  )
}

if (any(!is.finite(
  annual_kelvin
))) {
  stop(
    "The LME annual archive contains non-finite temperatures."
  )
}

analysis_year_indices <- match(
  analysis_years,
  archive_years
)

annual_celsius <- annual_kelvin[
  ,
  analysis_year_indices,
  ,
  drop = FALSE
] -
  273.15

lme_grid <- lme_coordinates %>%
  mutate(
    cell_id = row_number()
  ) %>%
  dplyr::select(
    cell_id,
    location_id,
    lat,
    long
  )

lme_long <- expand.grid(
  location_index = seq_len(
    nrow(
      lme_grid
    )
  ),
  year_index = seq_along(
    analysis_years
  ),
  member_index = seq_along(
    archive_members
  ),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    cell_id = lme_grid$cell_id[
      location_index
    ],
    location_id = lme_grid$location_id[
      location_index
    ],
    lat = lme_grid$lat[
      location_index
    ],
    long = lme_grid$long[
      location_index
    ],
    member = archive_members[
      member_index
    ],
    year = analysis_years[
      year_index
    ],
    temp_celsius = as.numeric(
      annual_celsius
    )
  )

if (
  nrow(
    lme_long
  ) !=
    nrow(
      lme_grid
    ) *
      length(
        analysis_years
      ) *
      number_of_lme_members
) {
  stop(
    "The long-format LME table has an unexpected number of rows."
  )
}

message(
  "LME grid cells found: ",
  nrow(
    lme_grid
  )
)

message(
  "Annual LME rows constructed: ",
  format(
    nrow(
      lme_long
    ),
    big.mark = ","
  )
)

# ------------------------------------------------------------
# 4. Construct LME climate-summary covariates
# ------------------------------------------------------------

# Long-run mean over all 13 members and all analysis years.
# The prepared archive is converted to degrees Celsius before
# constructing the climate summaries.
Lbar_df <- lme_long %>%
  group_by(
    cell_id
  ) %>%
  summarise(
    Lbar_g = mean(
      temp_celsius,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Member-specific anomalies.
lme_anomaly <- lme_long %>%
  group_by(
    cell_id,
    member
  ) %>%
  mutate(
    member_mean = mean(
      temp_celsius,
      na.rm = TRUE
    ),
    A = temp_celsius -
      member_mean
  ) %>%
  ungroup()

sigma_df <- lme_anomaly %>%
  group_by(
    cell_id
  ) %>%
  summarise(
    sigma2_g = var(
      A,
      na.rm = TRUE
    ),
    sigma_g = sd(
      A,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

cold_threshold <- as.numeric(
  quantile(
    lme_anomaly$A,
    probs = 0.10,
    na.rm = TRUE
  )
)

warm_threshold <- as.numeric(
  quantile(
    lme_anomaly$A,
    probs = 0.90,
    na.rm = TRUE
  )
)

extreme_df <- lme_anomaly %>%
  group_by(
    cell_id
  ) %>%
  summarise(
    pcold_g = mean(
      A <= cold_threshold,
      na.rm = TRUE
    ),
    pwarm_g = mean(
      A >= warm_threshold,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

climate_summary <- lme_grid %>%
  left_join(
    Lbar_df,
    by = "cell_id"
  ) %>%
  left_join(
    sigma_df,
    by = "cell_id"
  ) %>%
  left_join(
    extreme_df,
    by = "cell_id"
  )

readr::write_csv(
  climate_summary,
  file.path(
    output_intermediate_dir,
    "lme_climate_summaries.csv"
  )
)

# ------------------------------------------------------------
# 5. Construct documentary-coverage response R_gu
# ------------------------------------------------------------

reaches <- readxl::read_excel(
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

names(reaches) <- c(
  "level",
  "year",
  "long",
  "lat"
)

reaches_use <- reaches %>%
  transmute(
    year = as.integer(year),
    long = as.numeric(long),
    lat = as.numeric(lat),
    level = as.numeric(level)
  ) %>%
  filter(
    year %in% analysis_years,
    is.finite(level),
    is.finite(long),
    is.finite(lat)
  )

site_keys <- paste(
  sprintf(
    "%.8f",
    reaches_use$long
  ),
  sprintf(
    "%.8f",
    reaches_use$lat
  ),
  sep = "_"
)

reaches_use$site_id <- match(
  site_keys,
  unique(site_keys)
)

reaches_use$cell_id <- nearest_cell_ids(
  long = reaches_use$long,
  lat = reaches_use$lat,
  grid_df = lme_grid,
  chunk_size = nearest_cell_chunk_size
)

# Count distinct documentary sites, not duplicate records from
# the same site in the same year.
reaches_site_year <- reaches_use %>%
  distinct(
    site_id,
    cell_id,
    year
  )

coverage_observed <- reaches_site_year %>%
  count(
    cell_id,
    year,
    name = "Ngu"
  ) %>%
  mutate(
    R = as.integer(
      Ngu > 0L
    )
  )

coverage_panel <- tidyr::expand_grid(
  cell_id = lme_grid$cell_id,
  year = analysis_years
) %>%
  left_join(
    coverage_observed,
    by = c(
      "cell_id",
      "year"
    )
  ) %>%
  mutate(
    Ngu = replace_na(
      Ngu,
      0L
    ),
    R = replace_na(
      R,
      0L
    )
  ) %>%
  left_join(
    climate_summary,
    by = "cell_id"
  )

# ------------------------------------------------------------
# 6. Read, transform, and aggregate population rasters
# ------------------------------------------------------------

message(
  "Reading historical-population rasters..."
)

population_rasters <- lapply(
  population_years_available,
  function(population_year) {

    terra::rast(
      file.path(
        population_dir,
        paste0(
          population_year,
          "_pd.tif"
        )
      )
    )
  }
)

names(population_rasters) <- as.character(
  population_years_available
)

reference_raster <- population_rasters[[1]]

geometry_matches <- vapply(
  population_rasters[-1],
  function(current_raster) {

    isTRUE(
      terra::compareGeom(
        reference_raster,
        current_raster,
        stopOnError = FALSE
      )
    )
  },
  logical(1)
)

if (!all(geometry_matches)) {
  stop(
    "The historical-population rasters do not have identical ",
    "geometry and cannot be stacked directly."
  )
}

population_stack <- reference_raster

if (length(population_rasters) > 1L) {
  for (i in 2:length(population_rasters)) {
    population_stack <- c(
      population_stack,
      population_rasters[[i]]
    )
  }
}

names(population_stack) <- paste0(
  "pop_",
  population_years_available
)

population_xy <- as.data.frame(
  population_stack,
  xy = TRUE,
  na.rm = FALSE
) %>%
  filter(
    if_any(
      starts_with("pop_"),
      ~ !is.na(.x)
    )
  )

population_crs <- terra::crs(
  reference_raster
)

if (
  is.na(population_crs) ||
    !nzchar(population_crs)
) {
  stop(
    "The historical-population rasters do not have a defined CRS."
  )
}

population_sf <- sf::st_as_sf(
  population_xy,
  coords = c(
    "x",
    "y"
  ),
  crs = population_crs
)

population_lonlat_sf <- sf::st_transform(
  population_sf,
  4326
)

population_coordinates <- sf::st_coordinates(
  population_lonlat_sf
)

population_lonlat <- cbind(
  data.frame(
    long = population_coordinates[, 1],
    lat = population_coordinates[, 2]
  ),
  sf::st_drop_geometry(
    population_lonlat_sf
  )
)

population_lonlat$cell_id <- nearest_cell_ids(
  long = population_lonlat$long,
  lat = population_lonlat$lat,
  grid_df = lme_grid,
  chunk_size = nearest_cell_chunk_size
)

population_cell <- population_lonlat %>%
  group_by(
    cell_id
  ) %>%
  summarise(
    across(
      starts_with("pop_"),
      mean_or_na
    ),
    .groups = "drop"
  )

readr::write_csv(
  population_cell,
  file.path(
    output_intermediate_dir,
    "population_on_lme_grid.csv"
  )
)

# For years 1368--1911, use the temporally nearest one of
# 1776, 1820, 1851, 1880, and 1910. The 1953 raster is retained
# in the aggregated intermediate file but is not used here.
population_reference_years <- c(
  1776L,
  1820L,
  1851L,
  1880L,
  1910L
)

coverage_population_panel <- coverage_panel %>%
  left_join(
    population_cell,
    by = "cell_id"
  ) %>%
  rowwise() %>%
  mutate(
    nearest_pop_year = population_reference_years[
      which.min(
        abs(
          year -
            population_reference_years
        )
      )
    ],
    pop_hist = dplyr::case_when(
      nearest_pop_year == 1776L ~ pop_1776,
      nearest_pop_year == 1820L ~ pop_1820,
      nearest_pop_year == 1851L ~ pop_1851,
      nearest_pop_year == 1880L ~ pop_1880,
      nearest_pop_year == 1910L ~ pop_1910,
      TRUE ~ NA_real_
    ),
    Hgu = log1p(
      pop_hist
    )
  ) %>%
  ungroup()

# Common complete-case panel S used by all four models.
model_df <- coverage_population_panel %>%
  dplyr::select(
    cell_id,
    R,
    Ngu,
    long,
    lat,
    year,
    Lbar_g,
    sigma_g,
    pcold_g,
    pwarm_g,
    nearest_pop_year,
    pop_hist,
    Hgu
  ) %>%
  filter(
    complete.cases(
      across(
        c(
          R,
          long,
          lat,
          year,
          Lbar_g,
          sigma_g,
          pcold_g,
          pwarm_g,
          pop_hist,
          Hgu
        )
      )
    )
  ) %>%
  arrange(
    cell_id,
    year
  )

message(
  "Common complete-case model panel: ",
  nrow(model_df),
  " cell-years."
)

message(
  "Coverage outcomes: R=0: ",
  sum(model_df$R == 0L),
  "; R=1: ",
  sum(model_df$R == 1L)
)

if (nrow(model_df) != 63104L) {
  warning(
    "The manuscript reports n = 63,104, but this run produced n = ",
    nrow(model_df),
    ". Check input versions, raster CRS/extent, and LME grid ordering."
  )
}

saveRDS(
  model_df,
  file.path(
    output_intermediate_dir,
    "coverage_population_model_panel.rds"
  )
)

readr::write_csv(
  model_df,
  file.path(
    output_intermediate_dir,
    "coverage_population_model_panel.csv.gz"
  )
)

coverage_population_metadata <- data.frame(
  quantity = c(
    "input_mode",
    "lme_archive_file",
    "lme_archive_units",
    "lme_climate_units",
    "gam_fit_method",
    "analysis_first_year",
    "analysis_last_year",
    "number_of_lme_members",
    "number_of_lme_grid_cells",
    "number_of_distinct_reaches_sites",
    "number_of_reaches_site_years",
    "number_of_positive_coverage_cell_years",
    "number_of_model_cell_years",
    "number_of_model_grid_cells",
    "cold_anomaly_threshold",
    "warm_anomaly_threshold"
  ),
  value = c(
    input_mode,
    lme_archive_file,
    as.character(
      lme_archive$units
    ),
    "degrees Celsius",
    gam_fit_method,
    min(
      analysis_years
    ),
    max(
      analysis_years
    ),
    number_of_lme_members,
    nrow(
      lme_grid
    ),
    dplyr::n_distinct(
      reaches_use$site_id
    ),
    nrow(
      reaches_site_year
    ),
    sum(
      coverage_panel$R ==
        1L
    ),
    nrow(
      model_df
    ),
    dplyr::n_distinct(
      model_df$cell_id
    ),
    cold_threshold,
    warm_threshold
  )
)

readr::write_csv(
  coverage_population_metadata,
  file.path(
    output_intermediate_dir,
    "coverage_population_metadata.csv"
  )
)

# ------------------------------------------------------------
# 7. Fit the four comparable logistic GAMs
# ------------------------------------------------------------

message(
  "Fitting coverage GAM M1..."
)

m1 <- mgcv::gam(
  R ~
    s(
      long,
      lat,
      k = 50
    ) +
    s(
      year,
      k = 20
    ),
  family = binomial(),
  data = model_df,
  method = gam_fit_method
)

message(
  "Fitting coverage GAM M2..."
)

m2 <- mgcv::gam(
  R ~
    s(
      long,
      lat,
      k = 50
    ) +
    s(
      year,
      k = 20
    ) +
    Lbar_g +
    sigma_g +
    pcold_g +
    pwarm_g,
  family = binomial(),
  data = model_df,
  method = gam_fit_method
)

message(
  "Fitting coverage GAM M3..."
)

m3 <- mgcv::gam(
  R ~
    s(
      long,
      lat,
      k = 50
    ) +
    s(
      year,
      k = 20
    ) +
    Hgu,
  family = binomial(),
  data = model_df,
  method = gam_fit_method
)

message(
  "Fitting coverage GAM M4..."
)

m4 <- mgcv::gam(
  R ~
    s(
      long,
      lat,
      k = 50
    ) +
    s(
      year,
      k = 20
    ) +
    Lbar_g +
    sigma_g +
    pcold_g +
    pwarm_g +
    Hgu,
  family = binomial(),
  data = model_df,
  method = gam_fit_method
)

coverage_models <- list(
  M1 = m1,
  M2 = m2,
  M3 = m3,
  M4 = m4
)

saveRDS(
  coverage_models,
  file.path(
    output_intermediate_dir,
    "coverage_population_gam_fits.rds"
  )
)

# ------------------------------------------------------------
# 8. Construct manuscript Table 2
# ------------------------------------------------------------

additional_terms <- c(
  M1 = "None",
  M2 = "Climate",
  M3 = "Population",
  M4 = "Climate + Population"
)

model_table <- data.frame(
  Model = names(
    coverage_models
  ),
  Additional_terms = unname(
    additional_terms[
      names(coverage_models)
    ]
  ),
  n = vapply(
    coverage_models,
    stats::nobs,
    numeric(1)
  ),
  AIC = vapply(
    coverage_models,
    stats::AIC,
    numeric(1)
  ),
  BIC = vapply(
    coverage_models,
    stats::BIC,
    numeric(1)
  ),
  Deviance_explained = vapply(
    coverage_models,
    function(model) {
      summary(model)$dev.expl
    },
    numeric(1)
  ),
  row.names = NULL
)

model_table_formatted <- model_table %>%
  transmute(
    Model,
    `Additional terms` = Additional_terms,
    n = format(
      n,
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    ),
    AIC = sprintf(
      "%.2f",
      AIC
    ),
    BIC = sprintf(
      "%.2f",
      BIC
    ),
    `Dev. expl.` = sprintf(
      "%.2f%%",
      100 *
        Deviance_explained
    )
  )

readr::write_csv(
  model_table,
  file.path(
    output_table_dir,
    "Table2_coverage_GAM.csv"
  )
)

readr::write_csv(
  model_table_formatted,
  file.path(
    output_table_dir,
    "Table2_coverage_GAM_formatted.csv"
  )
)

# ------------------------------------------------------------
# 9. Nested model comparisons
# ------------------------------------------------------------

extract_nested_comparison <- function(
    reduced_model,
    full_model,
    comparison_name,
    added_terms) {

  comparison <- as.data.frame(
    anova(
      reduced_model,
      full_model,
      test = "Chisq"
    )
  )

  result_row <- nrow(
    comparison
  )

  df_column <- first_existing_column(
    comparison,
    c(
      "^Df$",
      "^df$"
    )
  )

  deviance_column <- first_existing_column(
    comparison,
    c(
      "^Deviance$",
      "Deviance"
    )
  )

  p_column <- first_existing_column(
    comparison,
    c(
      "Pr\\(>Chi\\)",
      "p.value",
      "p-value"
    )
  )

  data.frame(
    Comparison = comparison_name,
    Added_terms = added_terms,
    Df = if (
      is.na(df_column)
    ) {
      NA_real_
    } else {
      as.numeric(
        comparison[
          result_row,
          df_column
        ]
      )
    },
    Deviance_reduction = if (
      is.na(deviance_column)
    ) {
      NA_real_
    } else {
      as.numeric(
        comparison[
          result_row,
          deviance_column
        ]
      )
    },
    p_value = if (
      is.na(p_column)
    ) {
      NA_real_
    } else {
      as.numeric(
        comparison[
          result_row,
          p_column
        ]
      )
    }
  )
}

nested_comparisons <- bind_rows(
  extract_nested_comparison(
    m1,
    m2,
    "M2 versus M1",
    "Climate"
  ),
  extract_nested_comparison(
    m1,
    m3,
    "M3 versus M1",
    "Population"
  ),
  extract_nested_comparison(
    m2,
    m4,
    "M4 versus M2",
    "Population"
  ),
  extract_nested_comparison(
    m3,
    m4,
    "M4 versus M3",
    "Climate"
  )
)

readr::write_csv(
  nested_comparisons,
  file.path(
    output_table_dir,
    "Table2_nested_comparisons.csv"
  )
)

# ------------------------------------------------------------
# 10. Population coefficient from the full model
# ------------------------------------------------------------

m4_parametric_table <- as.data.frame(
  summary(m4)$p.table
)

if (!"Hgu" %in% rownames(m4_parametric_table)) {
  stop(
    "The Hgu coefficient was not found in the full-model summary."
  )
}

Hgu_row <- m4_parametric_table[
  "Hgu",
  ,
  drop = FALSE
]

Hgu_coefficient <- data.frame(
  Term = "Hgu",
  Estimate = Hgu_row[
    1,
    1
  ],
  Standard_error = Hgu_row[
    1,
    2
  ],
  Test_statistic = Hgu_row[
    1,
    3
  ],
  p_value = Hgu_row[
    1,
    4
  ]
)

readr::write_csv(
  Hgu_coefficient,
  file.path(
    output_table_dir,
    "Table2_Hgu_coefficient.csv"
  )
)

# ------------------------------------------------------------
# 11. Save a human-readable analysis report
# ------------------------------------------------------------

capture.output(
  {
    cat(
      "Documentary coverage and historical-population GAM analysis\n"
    )

    cat(
      "===========================================================\n\n"
    )

    cat(
      "Analysis years:",
      min(analysis_years),
      "--",
      max(analysis_years),
      "\n"
    )

    cat(
      "LME annual archive:",
      lme_archive_file,
      "\n"
    )

    cat(
      "LME grid cells:",
      nrow(lme_grid),
      "\n"
    )

    cat(
      "GAM fitting method:",
      gam_fit_method,
      "\n"
    )

    cat(
      "Common complete-case cell-years:",
      nrow(model_df),
      "\n"
    )

    cat(
      "Cold anomaly threshold:",
      cold_threshold,
      "\n"
    )

    cat(
      "Warm anomaly threshold:",
      warm_threshold,
      "\n\n"
    )

    print(
      model_table_formatted
    )

    cat(
      "\nNested comparisons\n"
    )

    print(
      nested_comparisons
    )

    cat(
      "\nHgu coefficient in M4\n"
    )

    print(
      Hgu_coefficient
    )

    cat(
      "\nM1 summary\n"
    )

    print(
      summary(m1)
    )

    cat(
      "\nM2 summary\n"
    )

    print(
      summary(m2)
    )

    cat(
      "\nM3 summary\n"
    )

    print(
      summary(m3)
    )

    cat(
      "\nM4 summary\n"
    )

    print(
      summary(m4)
    )

    cat(
      "\nSession information\n"
    )

    print(
      sessionInfo()
    )
  },
  file = file.path(
    output_intermediate_dir,
    "coverage_population_GAM_report.txt"
  )
)

message(
  "Coverage-population GAM analysis completed."
)

message(
  "Manuscript table saved to: ",
  file.path(
    output_table_dir,
    "Table2_coverage_GAM_formatted.csv"
  )
)

message(
  "Intermediate files saved to: ",
  output_intermediate_dir
)
