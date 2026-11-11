here::i_am("Code/DataPreparation/prepare_worldclim_bias_correction.R")

# WorldClim-based climatological bias correction for CESM-LME.
# Run AFTER DataPreparation/prepare_lme_annual.R, which performs
# the historical LME thin-plate spline interpolation at the cities.
#
# Default mode (reuse): use a previously computed 3-city WorldClim
# offset table. This mode NEVER reads the large NetCDF or GeoTIFF files.
#
# Optional mode (recompute): derive offsets from the 13 original
# modern CESM-LME NetCDF files and 12 WorldClim monthly GeoTIFFs.
# Set WORLDCLIM_OFFSET_MODE=recompute to enable that mode.
#
# Historical city data are always corrected AFTER the existing spline;
# the original LME files are never overwritten. No grid/offshore IDW
# correction is applied to historical outputs by this script.

# ------------------------------------------------------------
# 1. Configuration and historical-input checks
# ------------------------------------------------------------

offset_mode <- tolower(trimws(Sys.getenv("WORLDCLIM_OFFSET_MODE", "reuse")))
if (!offset_mode %in% c("reuse", "recompute")) {
  stop("WORLDCLIM_OFFSET_MODE must be 'reuse' or 'recompute'.")
}

# In a full upstream recalculation, an offshore IDW field can optionally
# be created for diagnostics. It is NEVER applied to the historical grid.
compute_offshore_idw_diagnostic <-
  tolower(Sys.getenv("WORLDCLIM_OFFSHORE_IDW_DIAGNOSTIC", "false")) %in%
  c("true", "1", "yes")

lme_dir <- here::here("Output", "Intermediate", "LME")
offset_dir <- here::here("Output", "LME_Bias_Correction")
corrected_dir <- here::here("Output", "Intermediate", "LME", "WorldClim_Corrected")
reference_offset_file <- here::here(
  "Data", "LME data", "worldclim_offsets_city3.csv"
)
generated_offset_file <- file.path(offset_dir, "worldclim_offsets_city3.csv")

dir.create(offset_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(corrected_dir, recursive = TRUE, showWarnings = FALSE)

city_rds <- file.path(lme_dir, "lme_city3_annual_1368_1911.rds")
if (!file.exists(city_rds)) {
  stop("Missing spline-interpolated historical city archive: ", city_rds,
       "\nRun Code/DataPreparation/prepare_lme_annual.R first.")
}
archive <- readRDS(city_rds)
required_archive <- c("coordinates", "years", "members", "units",
                      "annual_kelvin", "basis_dimension",
                      "fixed_smoothing_parameter")
if (!all(required_archive %in% names(archive))) {
  stop("Historical city archive is missing fields: ",
       paste(setdiff(required_archive, names(archive)), collapse = ", "))
}
if (!identical(as.character(archive$units), "Kelvin") ||
    !isTRUE(all.equal(as.numeric(archive$basis_dimension), 266)) ||
    !isTRUE(all.equal(as.numeric(archive$fixed_smoothing_parameter), 1e-8))) {
  stop("Unexpected units or thin-plate spline settings in historical city archive.")
}

cities <- as.data.frame(archive$coordinates)
if (!all(c("city", "long", "lati") %in% names(cities)) ||
    nrow(cities) != 3L || anyDuplicated(cities$city)) {
  stop("Historical city archive must contain three unique cities with coordinates.")
}
years_hist <- as.integer(archive$years)
members <- as.character(archive$members)
if (!identical(years_hist, 1368:1911) ||
    length(members) != 13L || anyDuplicated(members)) {
  stop("Historical city archive must contain 1368-1911 and 13 unique members.")
}
original_kelvin <- archive$annual_kelvin
if (!identical(as.integer(dim(original_kelvin)),
               c(nrow(cities), length(years_hist), length(members))) ||
    any(!is.finite(original_kelvin))) {
  stop("Malformed or nonfinite historical city LME array.")
}

# ------------------------------------------------------------
# 2. Optional upstream recomputation (large external files)
# ------------------------------------------------------------

if (offset_mode == "recompute") {
  required_packages <- c("ncdf4", "mgcv", "terra")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Install these packages to recompute offsets: ",
         paste(missing_packages, collapse = ", "))
  }

  nc_dir <- path.expand(Sys.getenv(
    "LME_NC_DIR", "~/Downloads/LME nc data"
  ))
  wc_dir <- path.expand(Sys.getenv(
    "WORLDCLIM_DIR", "~/Downloads/wc2.1_30s_tavg"
  ))

  modern_years <- 1970:2000
  lon_idx <- 40:58
  lat_idx <- 58:71
  nlon <- length(lon_idx)
  nlat <- length(lat_idx)
  monthly_start <- (1970L - 1850L) * 12L + 1L
  monthly_count <- length(modern_years) * 12L

  file_for <- function(j) file.path(
    nc_dir,
    sprintf("b.e11.BLMTRC5CN.f19_g16.%03d.cam.h0.TREFHT.185001-200512.nc", j)
  )

  f0 <- file_for(1L)
  if (!file.exists(f0)) stop("Missing modern LME NetCDF: ", f0)
  nc0 <- ncdf4::nc_open(f0)
  lon <- as.numeric(ncdf4::ncvar_get(nc0, "lon"))[lon_idx]
  lat <- as.numeric(ncdf4::ncvar_get(nc0, "lat"))[lat_idx]
  ncdf4::nc_close(nc0)
  coords <- data.frame(
    longitude = rep(lon, times = nlat),
    latitude = rep(lat, each = nlon)
  )
  if (nrow(coords) != 266L) stop("Expected exactly 266 native LME cells.")

  # Both periods must use the same native grid (order may differ).
  grid_rds <- file.path(lme_dir, "lme_annual_1368_1911.rds")
  if (!file.exists(grid_rds)) stop("Missing historical LME grid: ", grid_rds)
  historical_grid <- readRDS(grid_rds)
  if (!"coordinates" %in% names(historical_grid)) {
    stop("Historical LME grid has no coordinate metadata.")
  }
  old_coords <- historical_grid$coordinates
  if (!all(c("long", "lati") %in% names(old_coords))) {
    stop("Historical LME grid needs long and lati coordinates.")
  }
  coord_key <- function(lon, lat) sprintf("%.8f|%.8f", lon, lat)
  if (!setequal(coord_key(coords$longitude, coords$latitude),
                coord_key(old_coords$long, old_coords$lati))) {
    stop("Modern and historical native 266-cell LME grids do not match.")
  }

  grid <- data.frame(
    site = sprintf("grid_%03d", seq_len(nrow(coords))),
    longitude = coords$longitude,
    latitude = coords$latitude,
    type = "grid"
  )
  stations <- data.frame(
    site = as.character(cities$city),
    longitude = as.numeric(cities$long),
    latitude = as.numeric(cities$lati),
    type = "station"
  )
  locations <- rbind(grid, stations)
  if (anyDuplicated(locations$site)) stop("Duplicate location IDs.")
  utils::write.csv(locations, file.path(offset_dir, "locations_all_269.csv"),
                   row.names = FALSE)

  wc_files <- file.path(wc_dir, sprintf("wc2.1_30s_tavg_%02d.tif", 1:12))
  if (!all(file.exists(wc_files))) {
    stop("Missing WorldClim rasters: ",
         paste(wc_files[!file.exists(wc_files)], collapse = ", "))
  }
  pts <- terra::vect(locations, geom = c("longitude", "latitude"),
                     crs = "EPSG:4326")
  wc_monthly <- vapply(wc_files, function(f) {
    raster <- terra::rast(f)
    as.numeric(terra::extract(raster, pts, method = "simple")[, 2])
  }, numeric(nrow(locations)))
  wc_valid <- rowSums(is.finite(wc_monthly)) == 12L
  if (any(!wc_valid & locations$type == "station")) {
    stop("WorldClim is missing at a validation city.")
  }
  if (any(wc_monthly[wc_valid, ] < -80 | wc_monthly[wc_valid, ] > 60)) {
    stop("Unexpected WorldClim units or values.")
  }
  colnames(wc_monthly) <- sprintf("WorldClim_month_%02d_C", 1:12)
  locations$WorldClim_1970_2000_C <- rowMeans(wc_monthly)
  locations$worldclim_available <- wc_valid
  utils::write.csv(locations[!wc_valid, ],
                   file.path(offset_dir, "worldclim_masked_locations.csv"),
                   row.names = FALSE)
  utils::write.csv(cbind(locations, wc_monthly),
                   file.path(offset_dir, "worldclim_monthly_all_269.csv"),
                   row.names = FALSE)

  # Reusable thin-plate operator, matching historical bs='tp',
  # k=266, sp=1e-8 and the same native input coordinate support.
  train_coords <- data.frame(long = coords$longitude, lati = coords$latitude)
  smooth <- mgcv::smoothCon(
    mgcv::s(long, lati, bs = "tp", k = 266),
    data = train_coords, absorb.cons = TRUE, scale.penalty = TRUE
  )[[1]]
  train <- cbind(intercept = 1, smooth$X)
  target <- cbind(
    intercept = 1,
    mgcv::PredictMat(smooth, data.frame(
      long = locations$longitude, lati = locations$latitude
    ))
  )
  penalty <- matrix(0, ncol(train), ncol(train))
  penalty[-1, -1] <- smooth$S[[1]]
  normal <- crossprod(train) + 1e-8 * penalty
  if (!is.finite(rcond(normal)) || rcond(normal) < 1e-14) {
    stop("Modern thin-plate interpolation operator is ill-conditioned.")
  }
  operator <- target %*% solve(normal, t(train))

  lme_rows <- vector("list", 13L)
  for (j in seq_len(13L)) {
    f <- file_for(j)
    if (!file.exists(f)) stop("Missing modern LME NetCDF: ", f)
    message("Reading modern CESM-LME member a", j)
    nc <- ncdf4::nc_open(f)
    lon_j <- as.numeric(ncdf4::ncvar_get(nc, "lon"))[lon_idx]
    lat_j <- as.numeric(ncdf4::ncvar_get(nc, "lat"))[lat_idx]
    if (!isTRUE(all.equal(lon, lon_j, tolerance = 1e-10)) ||
        !isTRUE(all.equal(lat, lat_j, tolerance = 1e-10))) {
      ncdf4::nc_close(nc)
      stop("Modern member has a different coordinate grid: ", f)
    }
    unit_info <- ncdf4::ncatt_get(nc, "TREFHT", "units")
    if (!isTRUE(unit_info$hasatt)) {
      ncdf4::nc_close(nc)
      stop("Missing TREFHT units in: ", f)
    }
    arr <- ncdf4::ncvar_get(
      nc, "TREFHT",
      start = c(min(lon_idx), min(lat_idx), monthly_start),
      count = c(nlon, nlat, monthly_count),
      collapse_degen = FALSE
    )
    ncdf4::nc_close(nc)
    if (!identical(as.integer(dim(arr)), c(nlon, nlat, monthly_count))) {
      stop("Unexpected TREFHT dimensions in: ", f)
    }
    unit <- tolower(trimws(as.character(unit_info$value)))
    if (unit %in% c("k", "kelvin")) {
      arr <- arr - 273.15
    } else if (!unit %in% c("c", "degc", "degree_celsius", "degrees_celsius")) {
      stop("Unsupported TREFHT units: ", unit)
    }
    if (any(!is.finite(arr))) stop("Nonfinite modern TREFHT: ", f)
    native <- matrix(arr, nrow = nrow(coords), ncol = monthly_count)
    annual <- vapply(seq_along(modern_years), function(y) {
      months <- ((y - 1L) * 12L + 1L):(y * 12L)
      rowMeans(native[, months, drop = FALSE])
    }, numeric(nrow(coords)))
    target_annual <- operator %*% annual
    lme_rows[[j]] <- data.frame(
      site = rep(locations$site, times = length(modern_years)),
      member = paste0("a", j),
      year = rep(modern_years, each = nrow(locations)),
      temp_C = as.vector(target_annual)
    )
    rm(arr, native, annual, target_annual)
    gc(verbose = FALSE)
  }
  lme <- do.call(rbind, lme_rows)
  if (nrow(lme) != 13L * 31L * 269L || any(!is.finite(lme$temp_C))) {
    stop("Modern LME output is incomplete.")
  }
  utils::write.csv(lme,
                   file.path(offset_dir, "lme_annual_1970_2000_all_269.csv"),
                   row.names = FALSE)
  climatology <- stats::aggregate(temp_C ~ site, lme, mean)
  names(climatology)[2] <- "LME_1970_2000_C"
  offsets <- merge(locations, climatology, by = "site", sort = FALSE)
  offsets$d_C <- offsets$WorldClim_1970_2000_C - offsets$LME_1970_2000_C
  offsets$correction_method <- ifelse(
    offsets$worldclim_available, "WorldClim_direct", "unavailable"
  )
  offsets$distance_to_land_km <- ifelse(offsets$worldclim_available,
                                        0, NA_real_)

  if (compute_offshore_idw_diagnostic) {
    # Optional EXPLORATORY offset extrapolation. These results are never
    # used for the official three-city correction or spatial assimilation.
    source_idx <- which(offsets$worldclim_available & offsets$type == "grid")
    target_idx <- which(!offsets$worldclim_available & offsets$type == "grid")
    if (length(source_idx) < 4L) stop("Too few valid land grid cells for IDW.")
    for (i in target_idx) {
      lon1 <- offsets$longitude[i] * pi / 180
      lat1 <- offsets$latitude[i] * pi / 180
      lon2 <- offsets$longitude[source_idx] * pi / 180
      lat2 <- offsets$latitude[source_idx] * pi / 180
      a <- sin((lat2 - lat1) / 2)^2 +
        cos(lat1) * cos(lat2) * sin((lon2 - lon1) / 2)^2
      km <- 6371 * 2 * asin(pmin(1, sqrt(pmax(0, a))))
      nearest <- order(km)[seq_len(4L)]
      distance <- km[nearest]
      w <- 1 / pmax(distance, 1)^2
      offsets$d_C[i] <- sum(w * offsets$d_C[source_idx[nearest]]) / sum(w)
      offsets$correction_method[i] <- "IDW_diagnostic_only"
      offsets$distance_to_land_km[i] <- min(km)
    }
  }

  utils::write.csv(offsets,
                   file.path(offset_dir, "worldclim_offsets_all_269.csv"),
                   row.names = FALSE)
  utils::write.csv(offsets[offsets$type == "grid", ],
                   file.path(offset_dir, "worldclim_offsets_grid266.csv"),
                   row.names = FALSE)
  utils::write.csv(offsets[offsets$type == "station", ],
                   generated_offset_file, row.names = FALSE)
  offset_file <- generated_offset_file
  message("Recomputed WorldClim offsets; city table: ", offset_file)

} else {
  # REUSE mode: never load NetCDF or GeoTIFF data. Prefer a locally
  # regenerated offset CSV; otherwise use the small tracked reference.
  offset_override <- Sys.getenv("WORLDCLIM_OFFSETS_FILE", unset = "")
  if (nzchar(offset_override)) {
    offset_file <- path.expand(offset_override)
  } else if (file.exists(generated_offset_file)) {
    offset_file <- generated_offset_file
  } else if (file.exists(reference_offset_file)) {
    offset_file <- reference_offset_file
  } else {
    stop(
      "No cached WorldClim city offsets were found. Searched:\n  ",
      generated_offset_file, "\n  ", reference_offset_file,
      "\nCopy your previously computed city offset CSV to one of these paths ",
      "or set WORLDCLIM_OFFSET_MODE=recompute to regenerate it."
    )
  }
}

# ------------------------------------------------------------
# 3. Apply the three directly observed city offsets
# ------------------------------------------------------------

if (!file.exists(offset_file)) stop("Missing WorldClim offset CSV: ", offset_file)
offsets <- utils::read.csv(offset_file, stringsAsFactors = FALSE)
required_offsets <- c("site", "longitude", "latitude", "d_C",
                      "worldclim_available", "correction_method")
if (!all(required_offsets %in% names(offsets))) {
  stop("WorldClim offset CSV lacks columns: ",
       paste(setdiff(required_offsets, names(offsets)), collapse = ", "))
}
if (anyDuplicated(offsets$site)) stop("Duplicate sites in offset CSV.")
idx <- match(as.character(cities$city), offsets$site)
if (anyNA(idx)) stop("Not all three cities have WorldClim offsets.")
selected <- offsets[idx, , drop = FALSE]
if (!all(tolower(as.character(selected$worldclim_available)) == "true") ||
    !all(as.character(selected$correction_method) == "WorldClim_direct")) {
  stop("City offsets must be direct WorldClim values, not offshore IDW.")
}
if (any(!is.finite(selected$d_C)) ||
    any(abs(as.numeric(cities$long) - selected$longitude) > 1e-6) ||
    any(abs(as.numeric(cities$lati) - selected$latitude) > 1e-6)) {
  stop("City offsets have invalid values or mismatched coordinates.")
}

# A Celsius temperature DIFFERENCE has the same numerical size in Kelvin.
# Keep the original Kelvin interface so downstream code still subtracts 273.15.
corrected_kelvin <- sweep(original_kelvin, 1L, selected$d_C, FUN = "+")

summary_rows <- vector("list", nrow(cities))
for (j in seq_len(nrow(cities))) {
  before <- as.numeric(original_kelvin[j, , ])
  after <- as.numeric(corrected_kelvin[j, , ])
  max_error <- max(abs((after - before) - selected$d_C[j]))
  variance_error <- abs(stats::var(after) - stats::var(before))
  if (!is.finite(max_error) || max_error > 1e-10 ||
      !is.finite(variance_error) || variance_error > 1e-9) {
    stop("Additive correction validation failed at ", cities$city[j])
  }
  summary_rows[[j]] <- data.frame(
    city = as.character(cities$city[j]),
    offset_C = selected$d_C[j],
    original_mean_C = mean(before) - 273.15,
    corrected_mean_C = mean(after) - 273.15,
    max_addition_error = max_error,
    variance_difference = variance_error,
    correction_method = as.character(selected$correction_method[j]),
    stringsAsFactors = FALSE
  )
}

# expand.grid's first factor (city) varies fastest, matching the R array.
out <- expand.grid(
  city = as.character(cities$city),
  year = years_hist,
  member = members,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
pos <- match(out$city, cities$city)
out$long <- as.numeric(cities$long)[pos]
out$lati <- as.numeric(cities$lati)[pos]
out$temperature_kelvin <- as.vector(corrected_kelvin)
out <- out[c("city", "long", "lati", "member", "year", "temperature_kelvin")]
if (nrow(out) != 3L * length(years_hist) * length(members) ||
    anyDuplicated(out[c("city", "member", "year")]) ||
    any(!is.finite(out$temperature_kelvin))) {
  stop("Invalid corrected city CSV.")
}

corrected_archive <- archive
corrected_archive$annual_kelvin <- corrected_kelvin
corrected_archive$bias_correction <- list(
  method = "additive WorldClim-v2.1 reference level adjustment",
  reference_years = 1970:2000,
  reference_ensemble_members = 13L,
  offset_mode = offset_mode,
  offset_source = normalizePath(offset_file, mustWork = TRUE),
  city_offsets_C = stats::setNames(as.numeric(selected$d_C), as.character(cities$city)),
  original_historical_interpolation = "thin-plate regression spline; not repeated",
  applied_to = "three validation stations only; no offshore IDW"
)

utils::write.csv(out,
                 file.path(corrected_dir, "lme_city3_annual_1368_1911.csv"),
                 row.names = FALSE)
saveRDS(corrected_archive,
        file.path(corrected_dir, "lme_city3_annual_1368_1911.rds"),
        compress = "xz")
utils::write.csv(do.call(rbind, summary_rows),
                 file.path(corrected_dir, "worldclim_correction_diagnostics.csv"),
                 row.names = FALSE)
message("WorldClim mode: ", offset_mode,
        "; offset table: ", normalizePath(offset_file, mustWork = TRUE))
message("Corrected historical city LME written to: ", corrected_dir)
message("Original LME files unchanged; no offshore-grid corrections applied.")
print(do.call(rbind, summary_rows), row.names = FALSE)
