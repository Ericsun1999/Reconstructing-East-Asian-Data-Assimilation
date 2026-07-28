here::i_am("Code/Figure7e.R")

# ============================================================
# Quantile mapping and spatial smoothing for Figure 7(e)
#
# The original b0.csv contained the 13 LME members at native LME
# grid locations over land or within 50 km of China, Taiwan,
# Hong Kong, and Macao. This script reconstructs that selection.
#
# Generated inputs:
#   Output/Intermediate/REACHES/
#     reaches_kriging_lme_grid_mean.csv
#     reaches_kriging_lme_grid_variance.csv
#   Output/Intermediate/LME/lme_annual_1368_1911.rds
#
# Precomputed alternatives:
#   Data/REACHES/precomputed/
#     reaches_kriging_lme_grid_mean.csv
#     reaches_kriging_lme_grid_variance.csv
#   Data/LME data/precomputed/lme_annual_1368_1911.rds
#
# Outputs:
#   Output/Figure7-8/Figure7e.png
#   Output/Figure7-8/Figure7e_qm_locations.csv
#   Output/Figure7-8/Figure7e_smoothed_grid.csv
#   Output/Figure7-8/Figure7e_metadata.csv
# ============================================================

library(dplyr)
library(readr)
library(sf)
library(rnaturalearth)
library(stringr)
library(np)
library(mgcv)
library(ggplot2)
library(maps)

# ------------------------------------------------------------
# 1. Choose one complete input set
# ------------------------------------------------------------

input_mode <- "auto"
allowed_input_modes <- c("auto", "generated", "precomputed")

if (!input_mode %in% allowed_input_modes) {
  stop("input_mode must be one of: ",
       paste(allowed_input_modes, collapse = ", "))
}

generated_files <- c(
  reaches_mean = here::here(
    "Output", "Intermediate", "REACHES",
    "reaches_kriging_lme_grid_mean.csv"
  ),
  reaches_variance = here::here(
    "Output", "Intermediate", "REACHES",
    "reaches_kriging_lme_grid_variance.csv"
  ),
  lme_annual = here::here(
    "Output", "Intermediate", "LME",
    "lme_annual_1368_1911.rds"
  )
)

precomputed_files <- c(
  reaches_mean = here::here(
    "Data", "REACHES", "precomputed",
    "reaches_kriging_lme_grid_mean.csv"
  ),
  reaches_variance = here::here(
    "Data", "REACHES", "precomputed",
    "reaches_kriging_lme_grid_variance.csv"
  ),
  lme_annual = here::here(
    "Data", "LME data", "precomputed",
    "lme_annual_1368_1911.rds"
  )
)

select_input_set <- function(input_mode,
                             generated_files,
                             precomputed_files) {
  generated_complete <- all(file.exists(generated_files))
  precomputed_complete <- all(file.exists(precomputed_files))

  if (input_mode == "generated") {
    if (!generated_complete) {
      stop(
        "The generated Figure 7(e) input set is incomplete:\n  ",
        paste(
          generated_files[!file.exists(generated_files)],
          collapse = "\n  "
        )
      )
    }
    return(generated_files)
  }

  if (input_mode == "precomputed") {
    if (!precomputed_complete) {
      stop(
        "The precomputed Figure 7(e) input set is incomplete:\n  ",
        paste(
          precomputed_files[!file.exists(precomputed_files)],
          collapse = "\n  "
        )
      )
    }
    return(precomputed_files)
  }

  if (generated_complete) return(generated_files)
  if (precomputed_complete) return(precomputed_files)

  stop(
    "Neither a complete generated nor a complete precomputed ",
    "Figure 7(e) input set was found."
  )
}

input_files <- select_input_set(
  input_mode,
  generated_files,
  precomputed_files
)

message(
  "Figure 7(e) input selected from: ",
  dirname(unname(input_files["reaches_mean"]))
)

# ------------------------------------------------------------
# 2. Load the annual LME archive
# ------------------------------------------------------------

lme_archive <- readRDS(unname(input_files["lme_annual"]))

required_lme_objects <- c(
  "coordinates", "years", "members", "annual_kelvin"
)

missing_lme_objects <- setdiff(
  required_lme_objects,
  names(lme_archive)
)

if (length(missing_lme_objects) > 0L) {
  stop(
    "The annual LME archive is missing: ",
    paste(missing_lme_objects, collapse = ", ")
  )
}

lme_coordinates <- lme_archive$coordinates %>%
  transmute(
    location_id = as.integer(location_id),
    long = as.numeric(long),
    lat = as.numeric(lati)
  )

lme_years <- as.integer(lme_archive$years)
lme_members <- as.character(lme_archive$members)
annual_kelvin <- lme_archive$annual_kelvin

expected_dimensions <- c(
  nrow(lme_coordinates),
  length(lme_years),
  length(lme_members)
)

if (!identical(dim(annual_kelvin), expected_dimensions)) {
  stop("annual_kelvin dimensions do not match its metadata.")
}

if (length(lme_members) != 13L) {
  stop("Expected 13 LME members, but found ", length(lme_members), ".")
}

if (anyNA(lme_coordinates) ||
    anyDuplicated(lme_coordinates$location_id) ||
    anyDuplicated(lme_coordinates[c("long", "lat")]) ||
    anyNA(annual_kelvin)) {
  stop("The annual LME archive contains invalid or duplicated data.")
}

# ------------------------------------------------------------
# 3. Reconstruct the original land-plus-50-km location mask
# ------------------------------------------------------------

sf::sf_use_s2(TRUE)

world <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
)

target_regions <- world %>%
  filter(
    stringr::str_detect(
      admin,
      stringr::regex(
        "China|Taiwan|Hong Kong|Macao|Macau",
        ignore_case = TRUE
      )
    )
  ) %>%
  st_make_valid()

if (nrow(target_regions) == 0L) {
  stop("The requested East Asian regions were not found.")
}

target_union <- st_union(target_regions)

lme_points_sf <- lme_coordinates %>%
  st_as_sf(
    coords = c("long", "lat"),
    crs = 4326,
    remove = FALSE
  )

lme_selection_buffer_km <- 50

keep_lme_location <- lengths(
  st_is_within_distance(
    lme_points_sf,
    target_union,
    dist = units::set_units(lme_selection_buffer_km, "km")
  )
) > 0L

selected_lme_locations <- lme_coordinates[
  keep_lme_location,
  ,
  drop = FALSE
]


selected_lme_indices <- match(
  selected_lme_locations$location_id,
  lme_coordinates$location_id
)

# ------------------------------------------------------------
# 4. Read REACHES predictions at native LME locations
# ------------------------------------------------------------

read_reaches_lme_grid <- function(path) {
  data <- readr::read_csv(
    path,
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  if (!"lat" %in% names(data) && "lati" %in% names(data)) {
    data <- rename(data, lat = lati)
  }

  required <- c("location_id", "long", "lat")
  missing <- setdiff(required, names(data))

  if (length(missing) > 0L) {
    stop(
      "The file ", path, " is missing: ",
      paste(missing, collapse = ", ")
    )
  }

  year_columns <- grep(
    "^[Xx]?[0-9]{4}$",
    names(data),
    value = TRUE
  )

  if (length(year_columns) == 0L) {
    stop("No annual columns were detected in: ", path)
  }

  year_values <- as.integer(
    sub("^[Xx]", "", year_columns)
  )

  data <- data %>%
    mutate(
      location_id = as.integer(location_id),
      long = as.numeric(long),
      lat = as.numeric(lat)
    )

  if (anyNA(data[required]) ||
      anyDuplicated(data$location_id)) {
    stop(
      "The file ", path,
      " contains invalid or duplicated locations."
    )
  }

  list(
    data = data,
    year_map = setNames(
      year_columns,
      as.character(year_values)
    )
  )
}

mean_object <- read_reaches_lme_grid(
  unname(input_files["reaches_mean"])
)

variance_object <- read_reaches_lme_grid(
  unname(input_files["reaches_variance"])
)

common_years <- intersect(
  names(mean_object$year_map),
  names(variance_object$year_map)
)

common_years <- as.character(
  sort(as.integer(common_years))
)

if (length(common_years) == 0L) {
  stop("The REACHES mean and variance files have no common years.")
}

mean_year_columns <- unname(
  mean_object$year_map[common_years]
)

variance_year_columns <- unname(
  variance_object$year_map[common_years]
)

mean_rows <- match(
  selected_lme_locations$location_id,
  mean_object$data$location_id
)

variance_rows <- match(
  selected_lme_locations$location_id,
  variance_object$data$location_id
)

if (anyNA(mean_rows) || anyNA(variance_rows)) {
  stop(
    "At least one selected LME location is absent from the ",
    "REACHES-at-LME-grid files."
  )
}

coordinate_difference <- max(
  abs(
    mean_object$data$long[mean_rows] -
      selected_lme_locations$long
  ),
  abs(
    mean_object$data$lat[mean_rows] -
      selected_lme_locations$lat
  ),
  abs(
    variance_object$data$long[variance_rows] -
      selected_lme_locations$long
  ),
  abs(
    variance_object$data$lat[variance_rows] -
      selected_lme_locations$lat
  )
)

if (!is.finite(coordinate_difference) ||
    coordinate_difference > 1e-8) {
  stop(
    "The LME archive and REACHES-at-LME-grid coordinates ",
    "are not aligned."
  )
}

reaches_mean_matrix <- as.matrix(
  mean_object$data[
    mean_rows,
    mean_year_columns,
    drop = FALSE
  ]
)

reaches_variance_matrix <- as.matrix(
  variance_object$data[
    variance_rows,
    variance_year_columns,
    drop = FALSE
  ]
)

storage.mode(reaches_mean_matrix) <- "double"
storage.mode(reaches_variance_matrix) <- "double"

if (any(!is.finite(reaches_mean_matrix)) ||
    any(!is.finite(reaches_variance_matrix))) {
  stop("The selected REACHES matrices contain non-finite values.")
}

if (any(reaches_variance_matrix < -1e-10)) {
  stop("Negative REACHES prediction variances were found.")
}

reaches_sd_matrix <- sqrt(
  pmax(reaches_variance_matrix, 0)
)

# ------------------------------------------------------------
# 5. Local quantile-mapping functions
# ------------------------------------------------------------

build_Fx_inv_local <- function(x_s) {
  x_s <- as.numeric(x_s)
  x_s <- x_s[is.finite(x_s)]

  if (length(x_s) < 2L) {
    stop("Too few finite LME values for local quantile mapping.")
  }

  x_sd <- stats::sd(x_s)

  if (!is.finite(x_sd) || x_sd <= 0) {
    stop("The local LME sample has zero or invalid variation.")
  }

  bandwidth <- np::npudistbw(
    dat = data.frame(x = x_s)
  )

  Fx_hat <- function(q) {
    fit <- np::npudist(
      bws = bandwidth,
      edat = data.frame(x = as.numeric(q))
    )
    as.numeric(fitted(fit))
  }

  function(u) {
    u <- pmin(pmax(as.numeric(u), 1e-8), 1 - 1e-8)
    multiplier <- 5

    repeat {
      lower <- min(x_s) - multiplier * x_sd
      upper <- max(x_s) + multiplier * x_sd

      if (Fx_hat(lower) <= min(u) &&
          Fx_hat(upper) >= max(u)) {
        break
      }

      multiplier <- multiplier * 2

      if (multiplier > 80) {
        stop("Unable to bracket the requested local LME quantile.")
      }
    }

    vapply(
      u,
      function(current_u) {
        uniroot(
          function(q) Fx_hat(q) - current_u,
          interval = c(lower, upper),
          tol = 1e-6
        )$root
      },
      numeric(1)
    )
  }
}

FY_hat_one_location <- function(y, yhat, prediction_sd) {
  valid <- is.finite(yhat) & is.finite(prediction_sd)
  yhat <- as.numeric(yhat[valid])
  prediction_sd <- as.numeric(prediction_sd[valid])

  if (length(yhat) == 0L) {
    stop(
      "No valid paired REACHES means and standard deviations ",
      "were available at a location."
    )
  }

  prediction_sd <- pmax(prediction_sd, 1e-8)

  mean(
    pnorm((y - yhat) / prediction_sd)
  )
}

# ------------------------------------------------------------
# 6. Compute g_s(0) at the selected locations
# ------------------------------------------------------------

set.seed(10)

mapped_temperature_at_zero <- vapply(
  seq_len(nrow(selected_lme_locations)),
  function(selected_row) {
    lme_location_index <- selected_lme_indices[selected_row]

    local_lme_celsius <- as.numeric(
      annual_kelvin[
        lme_location_index,
        ,
        ,
        drop = FALSE
      ]
    ) - 273.15

    local_inverse <- build_Fx_inv_local(
      local_lme_celsius
    )

    zero_index_probability <- FY_hat_one_location(
      y = 0,
      yhat = reaches_mean_matrix[
        selected_row,
        ,
        drop = TRUE
      ],
      prediction_sd = reaches_sd_matrix[
        selected_row,
        ,
        drop = TRUE
      ]
    )

    local_inverse(zero_index_probability)
  },
  numeric(1)
)

qm_locations <- selected_lme_locations %>%
  mutate(
    mapped_temperature_at_index_zero_celsius =
      mapped_temperature_at_zero
  )

# ------------------------------------------------------------
# 7. Smooth to the 0.25-degree display grid
# ------------------------------------------------------------

# fx = TRUE preserves the fixed-degree-of-freedom thin-plate fit
# used in the original Figure 7(e) script.
smooth_fit <- mgcv::gam(
  mapped_temperature_at_index_zero_celsius ~
    s(long, lat, bs = "tp", fx = TRUE),
  data = qm_locations
)

smooth_grid <- expand.grid(
  lon = seq(98, 124.5, by = 0.25),
  lat = seq(18, 42.5, by = 0.25)
)

smooth_grid$temp0_c <- as.numeric(
  predict(
    smooth_fit,
    newdata = data.frame(
      long = smooth_grid$lon,
      lat = smooth_grid$lat
    )
  )
)

# Preserve the original near-land display mask separately from
# the 50-km mask used to choose the calibration locations.
smooth_display_buffer_km <- 1

smooth_points_sf <- smooth_grid %>%
  st_as_sf(
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )

keep_smooth_location <- lengths(
  st_is_within_distance(
    smooth_points_sf,
    target_union,
    dist = units::set_units(
      smooth_display_buffer_km,
      "km"
    )
  )
) > 0L

smooth_grid <- smooth_grid[
  keep_smooth_location,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 8. Save Figure 7(e) and numerical outputs
# ------------------------------------------------------------

output_directory <- here::here(
  "Output",
  "Figure7-8"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_csv(
  qm_locations,
  file.path(
    output_directory,
    "Figure7e_qm_locations.csv"
  )
)

readr::write_csv(
  smooth_grid,
  file.path(
    output_directory,
    "Figure7e_smoothed_grid.csv"
  )
)

metadata <- data.frame(
  input_mode = input_mode,
  number_of_native_lme_locations = nrow(lme_coordinates),
  lme_selection_buffer_km = lme_selection_buffer_km,
  number_of_selected_lme_locations =
    nrow(selected_lme_locations),
  number_of_lme_members = length(lme_members),
  lme_start_year = min(lme_years),
  lme_end_year = max(lme_years),
  number_of_reaches_event_years = length(common_years),
  smooth_display_buffer_km = smooth_display_buffer_km
)

readr::write_csv(
  metadata,
  file.path(
    output_directory,
    "Figure7e_metadata.csv"
  )
)

p_figure7e <- ggplot(
  smooth_grid,
  aes(x = lon, y = lat)
) +
  geom_point(
    aes(colour = temp0_c),
    size = 4.5,
    shape = 15
  ) +
  coord_map(
    xlim = c(98, 124.5),
    ylim = c(18, 42.5)
  ) +
  scale_colour_gradientn(
    colours = c(
      "blue", "cyan", "green", "yellow", "red"
    ),
    limits = c(-10, 25),
    oob = scales::squish,
    na.value = "transparent",
    guide = "colourbar"
  ) +
  borders(
    database = "world",
    xlim = c(76, 126),
    ylim = c(18, 45),
    fill = NA,
    colour = "grey50"
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = NULL
  ) +
  theme(
    text = element_text(size = 15),
    legend.position = "right"
  )

figure7e_output_file <- file.path(
  output_directory,
  "Figure7e.png"
)

ggsave(
  filename = figure7e_output_file,
  plot = p_figure7e,
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)

message("Figure 7(e) saved to: ", figure7e_output_file)
message(
  "Selected ",
  nrow(selected_lme_locations),
  " LME locations using the ",
  lme_selection_buffer_km,
  "-km land/coastal mask."
)
