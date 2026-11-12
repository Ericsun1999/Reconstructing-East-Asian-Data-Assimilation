here::i_am("Code/DataPreparation/prepare_worldclim_bias_correction.R")

# Apply precomputed WorldClim-v2.1 additive climatological offsets to
# historical CESM-LME temperatures. This script does NOT read NetCDFs or
# GeoTIFFs. The upstream preparation is Data/Get_data/GetWorldClim.R.
#
# Inputs:
#   Data/WorldClim/worldclim_offsets_all_269.csv (269 rows; tracked)
#   Output/Intermediate/LME/lme_annual_1368_1911.rds (native grid, Kelvin)
#   Output/Intermediate/LME/lme_city3_annual_1368_1911.rds (spline, Kelvin)
#
# Outputs (original historical LME remains unchanged):
#   Output/Intermediate/LME/WorldClim_Corrected/
#     lme_city3_annual_1368_1911.csv / .rds
#     lme_land_grid_corrected_annual_1368_1911.rds
#     worldclim_grid_correction_status.csv
#     worldclim_correction_diagnostics.csv
#
# IMPORTANT: Only locations with DIRECT WorldClim references are corrected.
# Offshore/masked grid nodes are not filled with IDW or silently treated
# as corrected. The land-grid RDS is a SUBSET product and cannot be used
# as a drop-in replacement for the existing 266-grid assimilation RDS.

suppressPackageStartupMessages(library(mgcv))
# ============================================================
# Local spline functions (self-contained; no external helper file).
# Identical implementation in GetWorldClim.R and this script.
# Native training: all 266 LME nodes; diagnostic scoring: direct land nodes.
# ============================================================

# Shared thin-plate regression-spline implementation for CESM-LME.
# Used by both the 1970-2000 upstream climatology preparation and the
# 1368-1911 historical WorldClim-correction/audit step.
# The independently prepared historical city RDS is explicitly verified
# against this operator before corrections are applied.

lme_coordinate_key <- function(lon, lat) {
  sprintf("%.8f|%.8f", as.numeric(lon), as.numeric(lat))
}

build_lme_spline_operator <- function(train_coords, target_coords) {
  stopifnot(requireNamespace("mgcv", quietly = TRUE))
  train_coords <- as.data.frame(train_coords)
  target_coords <- as.data.frame(target_coords)
  expected <- c("long", "lati")
  if (!all(expected %in% names(train_coords)) ||
      !all(expected %in% names(target_coords)) ||
      nrow(train_coords) != 266L ||
      any(!is.finite(as.matrix(train_coords[, expected]))) ||
      any(!is.finite(as.matrix(target_coords[, expected]))) ||
      anyDuplicated(lme_coordinate_key(train_coords$long, train_coords$lati))) {
    stop("Invalid or duplicated 266-node LME training coordinates.")
  }

  # Untransformed longitude / latitude in decimal degrees, EPSG:4326.
  # The original and reference-period fits use all 266 native grid nodes
  # with equally weighted observations (w_i = 1); no separate sp selection.
  basis_dimension <- 266L
  fixed_sp <- 1e-8
  smooth <- mgcv::smoothCon(
    mgcv::s(long, lati, bs = "tp", k = basis_dimension),
    data = train_coords, absorb.cons = TRUE, scale.penalty = TRUE
  )[[1]]
  if (length(smooth$S) != 1L) {
    stop("Expected one mgcv thin-plate spline penalty matrix.")
  }
  train_design <- cbind(intercept = 1, smooth$X)
  target_design <- cbind(
    intercept = 1,
    mgcv::PredictMat(smooth, target_coords)
  )
  penalty <- matrix(0, ncol(train_design), ncol(train_design))
  penalty[-1, -1] <- smooth$S[[1]]
  # Explicit weight specification: w_i = 1 for all 266 locations.
  # This reduces to crossprod(X), exactly as in prepare_lme_annual.R.
  weights <- rep(1, nrow(train_design))
  normal <- crossprod(train_design, weights * train_design) +
    fixed_sp * penalty
  conditioning <- rcond(normal)
  if (!is.finite(conditioning) || conditioning < 1e-14) {
    stop("LME thin-plate spline system is numerically singular.")
  }
  operator <- target_design %*% solve(
    normal,
    sweep(t(train_design), 2L, weights, "*")
  )
  # Every member and year has the same operator; Kelvin and Celsius
  # temperature differences have identical numerical magnitudes.
  list(
    operator = operator,
    basis = "thin-plate regression spline (mgcv bs=tp)",
    basis_dimension = basis_dimension,
    coefficient_dimension = ncol(train_design),
    smoothing_parameter = fixed_sp,
    fitting_weights = "equal, 1 at each native grid point",
    coordinate_specification =
      "WGS84 EPSG:4326; raw longitude/latitude decimal degrees; no transform",
    training_coordinate_signature = paste(
      lme_coordinate_key(train_coords$long, train_coords$lati),
      collapse = ";"
    ),
    reciprocal_condition_number = conditioning
  )
}

lme_spline_specification <- function(fit, period) {
  data.frame(
    period = as.character(period),
    basis = fit$basis,
    k = fit$basis_dimension,
    coefficient_dimension = fit$coefficient_dimension,
    sp = fit$smoothing_parameter,
    fitting_weights = fit$fitting_weights,
    coordinate_specification = fit$coordinate_specification,
    training_coordinate_signature = fit$training_coordinate_signature,
    number_of_native_grid_nodes = 266L,
    reciprocal_condition_number = fit$reciprocal_condition_number,
    stringsAsFactors = FALSE
  )
}

# Reconstruct with the FULL 266-node fitting grid, but score only the
# subset with a direct WorldClim land reference. We intentionally do not
# remove offshore nodes from the spline training data: doing so would
# change the spline basis and historical/reference comparability.
lme_spline_native_grid_errors <- function(predicted_native,
                                          original_native,
                                          years, members,
                                          period,
                                          diagnostic_rows) {
  predicted_native <- matrix(predicted_native, nrow = 266L)
  original_native <- matrix(original_native, nrow = 266L)
  expected_fields <- length(years) * length(members)
  if (!identical(dim(predicted_native), dim(original_native)) ||
      ncol(predicted_native) != expected_fields ||
      any(!is.finite(predicted_native)) ||
      any(!is.finite(original_native))) {
    stop("Invalid native-grid reconstruction diagnostic arrays.")
  }
  if (!is.numeric(diagnostic_rows) || length(diagnostic_rows) == 0L ||
      anyNA(diagnostic_rows) ||
      any(diagnostic_rows != as.integer(diagnostic_rows)) ||
      any(diagnostic_rows < 1L | diagnostic_rows > 266L) ||
      anyDuplicated(diagnostic_rows)) {
    stop("Expected nonempty, unique land-node diagnostic row indices within 1:266.")
  }
  # Compare interpolated and original UNCORRECTED native-grid temperatures.
  # The scored locations are specifically the direct-WorldClim land nodes.
  errors <- (predicted_native - original_native)[diagnostic_rows, , drop = FALSE]
  fields <- expand.grid(year = years, member = members,
                        KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)
  data.frame(
    period = as.character(period),
    diagnostic_scope = "WorldClim_direct_land_grid_only",
    number_of_training_nodes = 266L,
    number_of_scored_land_nodes = length(diagnostic_rows),
    member = fields$member,
    year = fields$year,
    rmse_C = sqrt(colMeans(errors^2)),
    max_abs_error_C = apply(abs(errors), 2L, max),
    mean_error_C = colMeans(errors),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# End local spline functions
# ============================================================

lme_dir <- here::here("Output", "Intermediate", "LME")
out_dir <- here::here("Output", "Intermediate", "LME", "WorldClim_Corrected")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

offset_override <- Sys.getenv("WORLDCLIM_OFFSETS_FILE", unset = "")
offset_file <- if (nzchar(offset_override)) {
  path.expand(offset_override)
} else {
  here::here("Data", "WorldClim", "worldclim_offsets_all_269.csv")
}
if (!file.exists(offset_file)) {
  stop("Missing prepared WorldClim offsets: ", offset_file,
       "\nRun Data/Get_data/GetWorldClim.R on the upstream inputs or ",
       "add the prepared CSV under Data/WorldClim/.")
}

offsets <- utils::read.csv(offset_file, stringsAsFactors = FALSE)
required <- c("site", "longitude", "latitude", "type",
              "WorldClim_1970_2000_C", "worldclim_available",
              "LME_1970_2000_C", "d_C", "correction_method")
if (!all(required %in% names(offsets))) {
  stop("Offset table missing columns: ",
       paste(setdiff(required, names(offsets)), collapse = ", "))
}
if (nrow(offsets) != 269L || anyDuplicated(offsets$site) ||
    any(!is.finite(offsets$longitude)) ||
    any(!is.finite(offsets$latitude))) {
  stop("Offsets must contain 269 unique, valid target locations.")
}

offsets$worldclim_available <- as.logical(offsets$worldclim_available)
if (anyNA(offsets$worldclim_available)) {
  stop("Invalid worldclim_available values.")
}
grid_offsets <- offsets[offsets$type == "grid", , drop = FALSE]
city_offsets <- offsets[offsets$type == "station", , drop = FALSE]
if (nrow(grid_offsets) != 266L || nrow(city_offsets) != 3L ||
    !setequal(city_offsets$site, c("HongKong", "Shanghai", "Beijing"))) {
  stop("Expected 266 grid nodes and three validation stations.")
}

# An official offset is available only at a direct WorldClim location.
direct <- offsets$worldclim_available &
  offsets$correction_method == "WorldClim_direct"
if (any(offsets$worldclim_available != direct) ||
    any(!is.finite(offsets$d_C[direct])) ||
    any(!is.na(offsets$d_C[!direct]))) {
  stop("Official offset CSV must have finite direct offsets and NA ",
       "for locations without valid WorldClim data (no IDW).")
}
if (any(abs(offsets$d_C[direct] -
            (offsets$WorldClim_1970_2000_C[direct] -
             offsets$LME_1970_2000_C[direct])) > 1e-7)) {
  stop("WorldClim offset does not equal the difference of climatologies.")
}

city_file <- file.path(lme_dir, "lme_city3_annual_1368_1911.rds")
grid_file <- file.path(lme_dir, "lme_annual_1368_1911.rds")
for (f in c(city_file, grid_file)) {
  if (!file.exists(f)) {
    stop("Missing historical LME archive: ", f,
         "\nRun Code/DataPreparation/prepare_lme_annual.R first.")
  }
}
city_archive <- readRDS(city_file)
grid_archive <- readRDS(grid_file)
if (!all(c("coordinates", "years", "members", "units", "annual_kelvin") %in%
         names(city_archive)) ||
    !all(c("coordinates", "years", "members", "units", "annual_kelvin") %in%
         names(grid_archive))) {
  stop("Historical LME archive metadata are incomplete.")
}
if (!identical(city_archive$units, "Kelvin") ||
    !identical(grid_archive$units, "Kelvin")) {
  stop("Original historical LME arrays must be in Kelvin.")
}
years <- as.integer(city_archive$years)
members <- as.character(city_archive$members)
if (!identical(years, 1368:1911) ||
    !identical(years, as.integer(grid_archive$years)) ||
    length(members) != 13L ||
    !identical(members, as.character(grid_archive$members))) {
  stop("Historical annual LME grids have incompatible year/member metadata.")
}

cities <- as.data.frame(city_archive$coordinates)
grid <- as.data.frame(grid_archive$coordinates)
if (nrow(cities) != 3L || nrow(grid) != 266L ||
    !all(c("city", "long", "lati") %in% names(cities)) ||
    !all(c("long", "lati") %in% names(grid))) {
  stop("Historical city or native-grid coordinate structure is invalid.")
}
if (!identical(as.integer(dim(city_archive$annual_kelvin)),
               c(3L, 544L, 13L)) ||
    !identical(as.integer(dim(grid_archive$annual_kelvin)),
               c(266L, 544L, 13L)) ||
    any(!is.finite(city_archive$annual_kelvin)) ||
    any(!is.finite(grid_archive$annual_kelvin))) {
  stop("Historical LME annual arrays have invalid dimensions or values.")
}

# Verify that the same native-grid spline reproduces the independently
# generated historical city archive (prepare_lme_annual.R). This checks
# more than matching k/sp: it checks numerical equivalence of predictions.
train_coords <- data.frame(long = grid$long, lati = grid$lati)
target_coords <- rbind(
  train_coords,
  data.frame(long = cities$long, lati = cities$lati)
)
spline_fit <- build_lme_spline_operator(train_coords, target_coords)
historical_specification <- lme_spline_specification(
  spline_fit, period = "historical_1368_1911"
)
utils::write.csv(
  historical_specification,
  file.path(out_dir, "lme_spline_specification_historical.csv"),
  row.names = FALSE
)

# Compare the modern period's saved audit specification, if it is present.
# Missing metadata does not mean a pass: it means the upstream modern
# extraction must be rerun once to complete the provenance audit.
modern_spec_file <- here::here(
  "Data", "WorldClim", "lme_spline_specification_modern.csv"
)
modern_audit_status <- "modern_audit_not_yet_run"
if (file.exists(modern_spec_file)) {
  modern_specification <- utils::read.csv(
    modern_spec_file, stringsAsFactors = FALSE
  )
  expected_columns <- c(
    "basis", "k", "coefficient_dimension", "sp", "fitting_weights",
    "coordinate_specification", "training_coordinate_signature",
    "number_of_native_grid_nodes"
  )
  if (nrow(modern_specification) != 1L ||
      !all(expected_columns %in% names(modern_specification))) {
    stop("Invalid saved modern LME spline specification.")
  }
  for (column in expected_columns) {
    if (!isTRUE(all.equal(modern_specification[[column]],
                          historical_specification[[column]],
                          tolerance = 1e-12,
                          check.attributes = FALSE))) {
      stop("Modern and historical LME interpolation specifications ",
           "differ in ", column,
           ". Recompute the modern WorldClim offsets using the common ",
           "spline before running historical correction.")
    }
  }
  modern_audit_status <- "modern_and_historical_settings_match"
} else {
  warning("Modern spline audit metadata not yet available. ",
          "Run Data/Get_data/GetWorldClim.R once with the original ",
          "modern NetCDF files to verify the modern/ historical fit ",
          "and record the modern native-grid reconstruction errors.")
}

# Match grid coordinates (rather than trusting numeric site IDs) before
# selecting land nodes for the native-grid diagnostics and bias correction.
# This does NOT drop sea nodes from the 266-node spline training fit.
key_historical <- lme_coordinate_key(grid$long, grid$lati)
key_offsets <- lme_coordinate_key(grid_offsets$longitude,
                                   grid_offsets$latitude)
if (anyDuplicated(key_historical) || anyDuplicated(key_offsets) ||
    !setequal(key_historical, key_offsets)) {
  stop("Modern and historical native-grid coordinates do not match.")
}
gi <- match(key_historical, key_offsets)
if (anyNA(gi)) stop("Unable to align grid offsets to native LME grid.")
matched <- grid_offsets[gi, , drop = FALSE]
land <- matched$worldclim_available &
  matched$correction_method == "WorldClim_direct"
land_diagnostic_rows <- which(land)
if (length(land_diagnostic_rows) == 0L) {
  stop("No directly referenced WorldClim land nodes for spline diagnostics.")
}

raw_history <- grid_archive$annual_kelvin
history_matrix <- matrix(raw_history, nrow = 266L)
reconstructed_native <-
  spline_fit$operator[seq_len(266L), , drop = FALSE] %*% history_matrix
spline_history <- array(reconstructed_native, dim = dim(raw_history))

historical_native_errors <- lme_spline_native_grid_errors(
  predicted_native = reconstructed_native,
  original_native = history_matrix,
  years = years,
  members = members,
  period = "historical_1368_1911",
  diagnostic_rows = land_diagnostic_rows
)
utils::write.csv(
  historical_native_errors,
  file.path(out_dir, "lme_spline_native_grid_errors_historical.csv"),
  row.names = FALSE
)

reconstructed_cities <-
  spline_fit$operator[267:269, , drop = FALSE] %*% history_matrix
historical_city_difference_C <- max(abs(
  reconstructed_cities - matrix(city_archive$annual_kelvin, nrow = 3L)
))
if (!is.finite(historical_city_difference_C) ||
    historical_city_difference_C > 1e-5) {
  stop("Historical city series cannot be reproduced using the common ",
       "LME spline settings. Maximum difference (C): ",
       signif(historical_city_difference_C, 6),
       ". Refit the historical city series before applying offsets.")
}

utils::write.csv(
  data.frame(
    native_grid_diagnostic_scope = "WorldClim_direct_land_grid_only",
    native_grid_training_nodes = 266L,
    native_grid_scored_land_nodes = length(land_diagnostic_rows),
    historical_native_grid_RMSE_C =
      sqrt(mean(historical_native_errors$rmse_C^2)),
    historical_native_grid_max_abs_error_C =
      max(historical_native_errors$max_abs_error_C),
    historical_city_max_abs_difference_C =
      historical_city_difference_C,
    modern_historical_settings_status = modern_audit_status
  ),
  file.path(out_dir, "lme_spline_audit_summary.csv"),
  row.names = FALSE
)
message("Historical LAND-only native-grid spline RMSE (pooled, ",
        length(land_diagnostic_rows), " scored land nodes): ",
        signif(sqrt(mean(historical_native_errors$rmse_C^2)), 5),
        " C; maximum absolute error: ",
        signif(max(historical_native_errors$max_abs_error_C), 5), " C")
message("Historical city spline check max abs difference: ",
        signif(historical_city_difference_C, 5), " C")
message("Historical / modern specification: ", modern_audit_status)

# Match by city names and coordinates; do not rely on row order.
ci <- match(as.character(cities$city), city_offsets$site)
if (anyNA(ci)) stop("A historical station is missing from the offsets.")
cs <- city_offsets[ci, , drop = FALSE]
if (any(!cs$worldclim_available) ||
    any(abs(cities$long - cs$longitude) > 1e-6) ||
    any(abs(cities$lati - cs$latitude) > 1e-6)) {
  stop("City offsets or station coordinates do not match historical LME.")
}

# Three cities were already interpolated with thin-plate splines by
# prepare_lme_annual.R. Add the common ensemble offset once in Kelvin.
original_city <- city_archive$annual_kelvin
corrected_city <- sweep(original_city, MARGIN = 1L, STATS = cs$d_C, FUN = "+")
summary_city <- lapply(seq_len(nrow(cities)), function(i) {
  old <- as.numeric(original_city[i, , ])
  new <- as.numeric(corrected_city[i, , ])
  err <- max(abs((new - old) - cs$d_C[i]))
  vdiff <- abs(stats::var(new) - stats::var(old))
  if (!is.finite(err) || err > 1e-10 ||
      !is.finite(vdiff) || vdiff > 1e-9) {
    stop("Additive city correction failed: ", cities$city[i])
  }
  data.frame(
    city = cities$city[i],
    offset_C = cs$d_C[i],
    original_mean_C = mean(old) - 273.15,
    corrected_mean_C = mean(new) - 273.15,
    maximum_addition_error = err,
    variance_difference = vdiff
  )
})

city_long <- expand.grid(
  city = as.character(cities$city), year = years, member = members,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
city_idx <- match(city_long$city, cities$city)
city_long$long <- cities$long[city_idx]
city_long$lati <- cities$lati[city_idx]
city_long$temperature_kelvin <- as.vector(corrected_city)
city_long <- city_long[
  c("city", "long", "lati", "member", "year", "temperature_kelvin")
]
if (nrow(city_long) != 3L * 544L * 13L ||
    anyDuplicated(city_long[c("city", "member", "year")])) {
  stop("Corrected city panel is incomplete or duplicated.")
}
city_corrected_archive <- city_archive
city_corrected_archive$annual_kelvin <- corrected_city
city_corrected_archive$bias_correction <- list(
  reference = "WorldClim v2.1 tavg (1970-2000)",
  method = "location-specific additive shift, common to all members and years",
  locations = "three validation stations, direct WorldClim only",
  offset_file = "Data/WorldClim/worldclim_offsets_all_269.csv"
)
utils::write.csv(city_long,
                 file.path(out_dir, "lme_city3_annual_1368_1911.csv"),
                 row.names = FALSE)
saveRDS(city_corrected_archive,
        file.path(out_dir, "lme_city3_annual_1368_1911.rds"),
        compress = "xz")
utils::write.csv(do.call(rbind, summary_city),
                 file.path(out_dir, "worldclim_correction_diagnostics.csv"),
                 row.names = FALSE)

# Native-grid correction: apply the SAME near-interpolating thin-plate
# operator to each historical annual grid field, then add WorldClim offsets
# only where the independent WorldClim land reference is available.
# This yields a separate 144-node (data-dependent) LAND-ONLY archive.
# Coordinates and land-mask were already matched and validated above.
# The historical spline field was computed on the full fitting grid.
land_rows <- land_diagnostic_rows
land_history <- spline_history[land_rows, , , drop = FALSE]
land_corrected <- sweep(land_history, 1L, matched$d_C[land_rows], "+")
for (i in seq_along(land_rows)) {
  err <- max(abs((land_corrected[i, , ] - land_history[i, , ]) -
                   matched$d_C[land_rows[i]]))
  if (!is.finite(err) || err > 1e-10) {
    stop("Native-grid WorldClim correction failed for: ",
         matched$site[land_rows[i]])
  }
}

land_archive <- list(
  coordinates = grid[land_rows, , drop = FALSE],
  years = years,
  members = members,
  units = "Kelvin",
  annual_kelvin = land_corrected,
  original_location_ids = if ("location_id" %in% names(grid))
    grid$location_id[land_rows] else land_rows,
  bias_correction = list(
    reference = "WorldClim v2.1 tavg (1970-2000)",
    applied_to = "direct WorldClim land nodes only; no offshore IDW",
    number_of_original_locations = 266L,
    number_of_corrected_locations = length(land_rows),
    interpolation = "thin-plate spline, k=266, sp=1e-8, applied BEFORE correction"
  )
)
saveRDS(land_archive,
        file.path(out_dir, "lme_land_grid_corrected_annual_1368_1911.rds"),
        compress = "xz")

status <- data.frame(
  original_location_id = if ("location_id" %in% names(grid))
    grid$location_id else seq_len(nrow(grid)),
  site = matched$site,
  longitude = matched$longitude,
  latitude = matched$latitude,
  worldclim_available = matched$worldclim_available,
  correction_method = matched$correction_method,
  offset_C = matched$d_C,
  applied_to_land_product = land,
  stringsAsFactors = FALSE
)
utils::write.csv(status,
                 file.path(out_dir, "worldclim_grid_correction_status.csv"),
                 row.names = FALSE)
message("WorldClim corrections loaded from: ",
        normalizePath(offset_file, winslash = "/"))
message("Saved three-city corrected LME (13 members, 544 years).")
message("Saved land-only corrected grid with ", length(land_rows),
        " directly referenced nodes; ", sum(!land), " unavailable nodes excluded.")
message("Original LME archives are unchanged; no offshore IDW was applied.")
