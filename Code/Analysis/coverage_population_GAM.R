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
# Expected raw-data locations:
#   Data/temperature index value.v1.xlsx
#
#   Data/LME data/population/a1.csv.gz, ..., a13.csv.gz
#
#   For convenience, the script also accepts uncompressed
#   a1.csv, ..., a13.csv and several alternate LME folders.
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
lme_first_year <- 1350L
number_of_lme_members <- 13L

reaches_file <- here::here(
  "Data",
  "temperature index value.v1.xlsx"
)

candidate_lme_dirs <- c(
  # Current repository location.
  here::here("Data", "LME data", "population"),

  # Also accept this literal spelling in case the folder is
  # actually named "LME dara" rather than "LME data".
  here::here("Data", "LME dara", "population"),

  # Backward-compatible alternatives.
  here::here("Data", "LME data"),
  here::here("Data", "LME"),
  here::here("Data", "LME_data")
)

candidate_population_dirs <- c(
  here::here("Data", "population"),
  here::here("Data", "Population")
)

output_table_dir <- here::here(
  "Output",
  "Tables"
)

output_intermediate_dir <- here::here(
  "Output",
  "Intermediate",
  "coverage_population"
)

dir.create(
  output_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_intermediate_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

find_lme_member_file <- function(
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

  existing_files <- candidates[
    file.exists(candidates)
  ]

  if (length(existing_files) == 0L) {
    return(NA_character_)
  }

  # Prefer the compressed file when both versions exist.
  existing_files[1]
}


valid_lme_dir <- vapply(
  candidate_lme_dirs,
  function(directory) {

    member_files <- vapply(
      seq_len(number_of_lme_members),
      function(member_id) {
        find_lme_member_file(
          directory = directory,
          member_id = member_id
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

if (!any(valid_lme_dir)) {
  stop(
    paste0(
      "Could not find all 13 LME files in one directory.\n",
      "Each member may be stored as a*.csv.gz or a*.csv.\n",
      "Directories searched:\n",
      paste(
        paste0(
          "  - ",
          candidate_lme_dirs
        ),
        collapse = "\n"
      )
    )
  )
}

lme_dir <- candidate_lme_dirs[
  which(valid_lme_dir)[1]
]

message(
  "Using LME files from: ",
  lme_dir
)

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

population_dir <- candidate_population_dirs[
  which(valid_population_dir)[1]
]

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
    grid_df) {

  if (length(long) != length(lat)) {
    stop(
      "The longitude and latitude vectors have different lengths."
    )
  }

  vapply(
    seq_along(long),
    function(i) {

      squared_distance <- (
        grid_df$long -
          long[i]
      )^2 +
        (
          grid_df$lat -
            lat[i]
        )^2

      grid_df$cell_id[
        which.min(
          squared_distance
        )
      ]
    },
    integer(1)
  )
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
# 3. Read and annualize the 13 LME ensemble members
# ------------------------------------------------------------

annualize_lme_member <- function(
    input_file,
    member_id,
    years_use,
    first_year = 1350L) {

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

  raw_data <- read.csv(
    input_connection,
    row.names = 1,
    check.names = FALSE
  )

  if (ncol(raw_data) < 14L) {
    stop(
      "LME file ",
      input_file,
      " does not contain two coordinate columns followed by ",
      "monthly temperature columns."
    )
  }

  coordinates <- raw_data[
    ,
    1:2,
    drop = FALSE
  ]

  names(coordinates) <- c(
    "lat",
    "long"
  )

  coordinates$lat <- as.numeric(
    coordinates$lat
  )

  coordinates$long <- as.numeric(
    coordinates$long
  )

  monthly_values <- data.matrix(
    raw_data[
      ,
      -c(1, 2),
      drop = FALSE
    ]
  )

  number_of_month_columns <- ncol(
    monthly_values
  )

  if (
    number_of_month_columns %% 12L !=
      0L
  ) {
    stop(
      "The number of monthly columns in ",
      input_file,
      " is not divisible by 12."
    )
  }

  number_of_years <- number_of_month_columns /
    12L

  available_years <- first_year +
    seq_len(number_of_years) -
    1L

  if (!all(years_use %in% available_years)) {
    stop(
      input_file,
      " does not contain all requested years ",
      min(years_use),
      "--",
      max(years_use),
      "."
    )
  }

  annual_values <- t(
    vapply(
      seq_len(
        nrow(monthly_values)
      ),
      function(row_index) {

        monthly_matrix <- matrix(
          monthly_values[
            row_index,
          ],
          nrow = 12L
        )

        colMeans(
          monthly_matrix
        )
      },
      numeric(number_of_years)
    )
  )

  year_indices <- match(
    years_use,
    available_years
  )

  annual_values <- annual_values[
    ,
    year_indices,
    drop = FALSE
  ]

  data.frame(
    lat = rep(
      coordinates$lat,
      each = length(years_use)
    ),
    long = rep(
      coordinates$long,
      each = length(years_use)
    ),
    member = member_id,
    year = rep(
      years_use,
      times = nrow(coordinates)
    ),
    temp = as.numeric(
      t(
        annual_values
      )
    )
  )
}


message(
  "Reading and annualizing ",
  number_of_lme_members,
  " LME ensemble members..."
)

lme_long_list <- vector(
  "list",
  number_of_lme_members
)

for (member_id in seq_len(number_of_lme_members)) {

  input_file <- find_lme_member_file(
    directory = lme_dir,
    member_id = member_id
  )

  if (is.na(input_file)) {
    stop(
      "The LME file for member ",
      member_id,
      " disappeared after the initial input check."
    )
  }

  message(
    "  LME member ",
    member_id,
    "/",
    number_of_lme_members
  )

  lme_long_list[[member_id]] <- annualize_lme_member(
    input_file = input_file,
    member_id = member_id,
    years_use = analysis_years,
    first_year = lme_first_year
  )
}

lme_long <- bind_rows(
  lme_long_list
)

if (
  any(!is.finite(lme_long$lat)) ||
    any(!is.finite(lme_long$long)) ||
    any(!is.finite(lme_long$year)) ||
    any(!is.finite(lme_long$temp))
) {
  stop(
    "Non-finite values were found in the annualized LME data."
  )
}

lme_grid <- lme_long %>%
  distinct(
    lat,
    long
  ) %>%
  arrange(
    lat,
    long
  ) %>%
  mutate(
    cell_id = row_number()
  )

lme_long <- lme_long %>%
  left_join(
    lme_grid,
    by = c(
      "lat",
      "long"
    )
  )

message(
  "LME grid cells found: ",
  nrow(lme_grid)
)

# ------------------------------------------------------------
# 4. Construct LME climate-summary covariates
# ------------------------------------------------------------

# Long-run mean over all 13 members and all analysis years.
#
# The original analysis retained the LME temperature unit in
# the source files here. Subtracting a common Kelvin-to-Celsius
# constant would only shift this linear covariate and intercept.
Lbar_df <- lme_long %>%
  group_by(
    cell_id
  ) %>%
  summarise(
    Lbar_g = mean(
      temp,
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
      temp,
      na.rm = TRUE
    ),
    A = temp -
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
  grid_df = lme_grid
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
  grid_df = lme_grid
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
  select(
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
  gzfile(
    file.path(
      output_intermediate_dir,
      "coverage_population_model_panel.csv.gz"
    )
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
  method = "REML"
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
  method = "REML"
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
  method = "REML"
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
  method = "REML"
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
      "LME grid cells:",
      nrow(lme_grid),
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
