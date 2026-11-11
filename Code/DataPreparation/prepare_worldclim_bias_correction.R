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
coord_key <- function(lon, lat) {
  sprintf("%.8f|%.8f", as.numeric(lon), as.numeric(lat))
}
key_historical <- coord_key(grid$long, grid$lati)
key_offsets <- coord_key(grid_offsets$longitude, grid_offsets$latitude)
if (anyDuplicated(key_historical) || anyDuplicated(key_offsets) ||
    !setequal(key_historical, key_offsets)) {
  stop("Modern and historical native-grid coordinates do not match.")
}
gi <- match(key_historical, key_offsets)
if (anyNA(gi)) stop("Unable to align grid offsets to native LME grid.")
matched <- grid_offsets[gi, , drop = FALSE]
land <- matched$worldclim_available &
  matched$correction_method == "WorldClim_direct"
if (!any(land)) stop("No directly referenced WorldClim land nodes.")

# Identical full-rank thin-plate setup (k = 266, sp = 1e-8).
train_coords <- data.frame(long = grid$long, lati = grid$lati)
smooth <- mgcv::smoothCon(
  mgcv::s(long, lati, bs = "tp", k = 266),
  data = train_coords, absorb.cons = TRUE, scale.penalty = TRUE
)[[1]]
training <- cbind(intercept = 1, smooth$X)
target <- cbind(intercept = 1,
                mgcv::PredictMat(smooth, train_coords))
penalty <- matrix(0, ncol(training), ncol(training))
penalty[-1, -1] <- smooth$S[[1]]
normal <- crossprod(training) + 1e-8 * penalty
if (!is.finite(rcond(normal)) || rcond(normal) < 1e-14) {
  stop("Historical thin-plate interpolation matrix is singular.")
}
operator <- target %*% solve(normal, t(training))

raw_history <- grid_archive$annual_kelvin
spline_history <- operator %*% matrix(raw_history, nrow = 266L)
dim(spline_history) <- dim(raw_history)
land_rows <- which(land)
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
