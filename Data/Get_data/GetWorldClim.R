here::i_am("Data/Get_data/GetWorldClim.R")

# Optional UPSTREAM preparation (not part of run_all.R).
# Inputs (download separately):
#   - 13 CESM-LME monthly TREFHT NetCDF files covering 1850-2005
#   - WorldClim v2.1 tavg 30-second climatology (12 GeoTIFFs; 1970-2000)
# Outputs:
#   - Data/WorldClim/worldclim_offsets_all_269.csv (small, version-controlled)
#   - Output/LME_Bias_Correction/ (local diagnostics)
#
# Spline and reference settings match prepare_lme_annual.R:
# bs = "tp", k = 266, fixed sp = 1e-8, 12-month arithmetic annual mean.
# No offshore IDW is used for the official correction table.

suppressPackageStartupMessages({
  library(ncdf4)
  library(mgcv)
  library(terra)
})
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

nc_dir <- path.expand(Sys.getenv("LME_NC_DIR", "~/Downloads/LME nc data"))
wc_dir <- path.expand(Sys.getenv("WORLDCLIM_DIR", "~/Downloads/wc2.1_30s_tavg"))
prepared_dir <- here::here("Data", "WorldClim")
diagnostic_dir <- here::here("Output", "LME_Bias_Correction")
dir.create(prepared_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

members <- paste0("a", 1:13)
reference_years <- 1970:2000
lon_idx <- 40:58
lat_idx <- 58:71
n_lon <- length(lon_idx)
n_lat <- length(lat_idx)
stopifnot(n_lon * n_lat == 266L)
monthly_start <- (1970L - 1850L) * 12L + 1L
monthly_count <- length(reference_years) * 12L

netcdf_file <- function(member_number) {
  file.path(
    nc_dir,
    sprintf("b.e11.BLMTRC5CN.f19_g16.%03d.cam.h0.TREFHT.185001-200512.nc",
            member_number)
  )
}

# Read native-grid coordinates, with longitude varying fastest.
first_file <- netcdf_file(1L)
if (!file.exists(first_file)) stop("Missing NetCDF: ", first_file)
nc <- ncdf4::nc_open(first_file)
lon <- as.numeric(ncdf4::ncvar_get(nc, "lon"))[lon_idx]
lat <- as.numeric(ncdf4::ncvar_get(nc, "lat"))[lat_idx]
ncdf4::nc_close(nc)

native_grid <- data.frame(
  site = sprintf("grid_%03d", 1:266),
  longitude = rep(lon, times = n_lat),
  latitude = rep(lat, each = n_lon),
  type = "grid",
  stringsAsFactors = FALSE
)
stations <- data.frame(
  site = c("HongKong", "Shanghai", "Beijing"),
  longitude = c(114.167, 121.433, 116.283),
  latitude = c(22.333, 31.167, 39.933),
  type = "station",
  stringsAsFactors = FALSE
)
locations <- rbind(native_grid, stations)
stopifnot(nrow(locations) == 269L, !anyDuplicated(locations$site))

# Validate that the grid is the one used by the historical LME archive, if present.
historical_file <- here::here("Output", "Intermediate", "LME",
                             "lme_annual_1368_1911.rds")
if (file.exists(historical_file)) {
  historical <- readRDS(historical_file)
  hc <- historical$coordinates
  if (is.null(hc) || nrow(hc) != 266L ||
      !all(c("long", "lati") %in% names(hc))) {
    stop("Historical LME coordinate metadata are not as expected.")
  }
  key <- function(x, y) sprintf("%.8f|%.8f", as.numeric(x), as.numeric(y))
  if (!setequal(key(hc$long, hc$lati),
                key(native_grid$longitude, native_grid$latitude))) {
    stop("Modern and historical native LME coordinates do not match.")
  }
  if (!identical(
        lme_coordinate_key(hc$long, hc$lati),
        lme_coordinate_key(native_grid$longitude, native_grid$latitude))) {
    stop("Historical and modern native LME coordinates match as sets but ",
         "not in row order. Align the training order before using the ",
         "same spline operator in both periods.")
  }
  message("Modern/historical native LME coordinates and order match (266 locations).")
} else {
  message("Historical LME archive not found; coordinate cross-check skipped.")
}

# WorldClim tavg (1970-2000) climatology: 12 equal-weight monthly means.
wc_files <- file.path(wc_dir, sprintf("wc2.1_30s_tavg_%02d.tif", 1:12))
missing_wc <- wc_files[!file.exists(wc_files)]
if (length(missing_wc)) stop("Missing WorldClim GeoTIFF(s): ",
                            paste(missing_wc, collapse = ", "))
points <- terra::vect(locations, geom = c("longitude", "latitude"),
                      crs = "EPSG:4326")
wc_monthly <- vapply(wc_files, function(f) {
  raster <- terra::rast(f)
  as.numeric(terra::extract(raster, points, method = "simple")[[2]])
}, numeric(nrow(locations)))
colnames(wc_monthly) <- sprintf("WorldClim_month_%02d_C", 1:12)
wc_available <- rowSums(is.finite(wc_monthly)) == 12L
if (any(!wc_available[locations$type == "station"])) {
  stop("WorldClim is missing at a validation station.")
}
if (any(wc_monthly[wc_available, , drop = FALSE] < -80 |
        wc_monthly[wc_available, , drop = FALSE] > 60)) {
  stop("Inspect WorldClim temperature units and GeoTIFF scaling.")
}
locations$WorldClim_1970_2000_C <- rowMeans(wc_monthly)
locations$worldclim_available <- wc_available
# For the spline reconstruction-error audit, score only native GRID nodes
# with a direct land reference. Keep ALL 266 nodes in the spline fit.
land_diagnostic_rows <- which(wc_available[seq_len(nrow(native_grid))])
if (length(land_diagnostic_rows) == 0L) {
  stop("No native-grid land locations have valid WorldClim reference values.")
}
message("Spline fitting uses 266 native nodes; native reconstruction ",
        "errors will be scored only at ", length(land_diagnostic_rows),
        " directly referenced land nodes.")
utils::write.csv(cbind(locations, wc_monthly),
                 file.path(diagnostic_dir, "worldclim_monthly_all_269.csv"),
                 row.names = FALSE)

# Same operator constructor is also used to audit historical fields.
train_coords <- data.frame(long = native_grid$longitude,
                           lati = native_grid$latitude)
target_coords <- data.frame(long = locations$longitude,
                            lati = locations$latitude)
spline_fit <- build_lme_spline_operator(train_coords, target_coords)
operator <- spline_fit$operator

# Modern CESM-LME annual values: all 13 members, each year 1970-2000.
modern_records <- vector("list", 13L)
modern_native_errors <- vector("list", 13L)
for (j in 1:13) {
  path <- netcdf_file(j)
  if (!file.exists(path)) stop("Missing NetCDF: ", path)
  message("Reading reference-period TREFHT: ", basename(path))
  nc <- ncdf4::nc_open(path)
  lon_j <- as.numeric(ncdf4::ncvar_get(nc, "lon"))[lon_idx]
  lat_j <- as.numeric(ncdf4::ncvar_get(nc, "lat"))[lat_idx]
  if (!isTRUE(all.equal(lon, lon_j, tolerance = 1e-10)) ||
      !isTRUE(all.equal(lat, lat_j, tolerance = 1e-10))) {
    ncdf4::nc_close(nc)
    stop("Modern LME member coordinates differ: ", path)
  }
  u <- ncdf4::ncatt_get(nc, "TREFHT", "units")
  if (!isTRUE(u$hasatt)) {
    ncdf4::nc_close(nc)
    stop("Missing TREFHT units metadata: ", path)
  }
  arr <- ncdf4::ncvar_get(
    nc, "TREFHT",
    start = c(min(lon_idx), min(lat_idx), monthly_start),
    count = c(n_lon, n_lat, monthly_count),
    collapse_degen = FALSE
  )
  ncdf4::nc_close(nc)
  if (!identical(as.integer(dim(arr)), c(n_lon, n_lat, monthly_count))) {
    stop("Unexpected TREFHT dimensions: ", path)
  }
  unit <- tolower(trimws(as.character(u$value)))
  if (unit %in% c("k", "kelvin")) {
    arr <- arr - 273.15
  } else if (!unit %in% c("c", "degc", "degree_celsius", "degrees_celsius")) {
    stop("Unexpected TREFHT units: ", unit)
  }
  if (any(!is.finite(arr))) stop("Invalid TREFHT values in: ", path)
  native_monthly <- matrix(arr, nrow = nrow(native_grid),
                           ncol = monthly_count)
  annual <- vapply(seq_along(reference_years), function(y) {
    cols <- ((y - 1L) * 12L + 1L):(y * 12L)
    rowMeans(native_monthly[, cols, drop = FALSE])
  }, numeric(nrow(native_grid)))
  annual_locations <- operator %*% annual
  # At the first 266 targets, compare spline predictions against the
  # original annual native-grid temperatures, not against corrected values.
  modern_native_errors[[j]] <- lme_spline_native_grid_errors(
    predicted_native = annual_locations[seq_len(266L), , drop = FALSE],
    original_native = annual,
    years = reference_years,
    members = members[j],
    period = "modern_1970_2000",
    diagnostic_rows = land_diagnostic_rows
  )
  modern_records[[j]] <- data.frame(
    site = rep(locations$site, times = length(reference_years)),
    member = members[j],
    year = rep(reference_years, each = nrow(locations)),
    temp_C = as.vector(annual_locations),
    stringsAsFactors = FALSE
  )
  rm(arr, native_monthly, annual, annual_locations)
}
modern <- do.call(rbind, modern_records)
modern_errors <- do.call(rbind, modern_native_errors)
utils::write.csv(
  modern_errors,
  file.path(diagnostic_dir, "lme_spline_native_grid_errors_modern.csv"),
  row.names = FALSE
)
message("Modern LAND-only native-grid spline RMSE (pooled, ",
        length(land_diagnostic_rows), " scored land nodes): ",
        signif(sqrt(mean(modern_errors$rmse_C^2)), 5), " C; max abs: ",
        signif(max(modern_errors$max_abs_error_C), 5), " C")
stopifnot(nrow(modern) == 269L * 13L * 31L,
          all(is.finite(modern$temp_C)))
utils::write.csv(modern,
                 file.path(diagnostic_dir, "lme_annual_1970_2000_all_269.csv"),
                 row.names = FALSE)

clim <- stats::aggregate(temp_C ~ site, data = modern, FUN = mean)
names(clim)[2] <- "LME_1970_2000_C"
idx <- match(locations$site, clim$site)
if (anyNA(idx)) stop("Reference-period LME mean missing for a site.")
locations$LME_1970_2000_C <- clim$LME_1970_2000_C[idx]
locations$d_C <- locations$WorldClim_1970_2000_C - locations$LME_1970_2000_C

# No fabricated WorldClim ocean values. The 122 (data-dependent) masked
# sites remain NA for official corrections. No offshore IDW is applied.
locations$correction_method <- ifelse(
  locations$worldclim_available, "WorldClim_direct", "unavailable"
)
locations$d_C[!locations$worldclim_available] <- NA_real_
if (!all(is.finite(locations$d_C[locations$worldclim_available]))) {
  stop("Nonfinite direct WorldClim correction.")
}

# Compare regenerated offsets with the previously tracked prepared input.
# If they differ, the newly computed values replace the old file, and the
# change is recorded for review before committing the new input to GitHub.
prepared_file <- file.path(prepared_dir, "worldclim_offsets_all_269.csv")
if (file.exists(prepared_file)) {
  old <- utils::read.csv(prepared_file, stringsAsFactors = FALSE)
  if (all(c("site", "d_C", "worldclim_available") %in% names(old)) &&
      !anyDuplicated(old$site)) {
    index <- match(locations$site, old$site)
    if (!anyNA(index)) {
      difference <- locations$d_C - old$d_C[index]
      comparison <- data.frame(
        site = locations$site,
        old_offset_C = old$d_C[index],
        regenerated_offset_C = locations$d_C,
        difference_C = difference,
        old_worldclim_available = old$worldclim_available[index],
        new_worldclim_available = locations$worldclim_available
      )
      utils::write.csv(
        comparison,
        file.path(diagnostic_dir, "lme_spline_offset_before_after.csv"),
        row.names = FALSE
      )
      comparable <- is.finite(difference)
      if (any(comparable) && max(abs(difference[comparable])) > 1e-5) {
        warning("Recomputed WorldClim offsets differ from the previous ",
                "prepared input. Review lme_spline_offset_before_after.csv ",
                "and commit the regenerated CSV after checking the changes.")
      }
      if (any(as.logical(old$worldclim_available[index]) !=
              locations$worldclim_available)) {
        warning("WorldClim availability masks have changed since the ",
                "previous prepared offsets.")
      }
    }
  }
}
# Small prepared input intended for version control.
utils::write.csv(locations, prepared_file, row.names = FALSE)
# Mark the modern spline audit complete only after all 13 NetCDF members
# have been processed and the regenerated offset CSV has been saved.
modern_specification <- lme_spline_specification(
  spline_fit, period = "modern_1970_2000"
)
utils::write.csv(
  modern_specification,
  file.path(diagnostic_dir, "lme_spline_specification_modern.csv"),
  row.names = FALSE
)
utils::write.csv(
  modern_specification,
  file.path(prepared_dir, "lme_spline_specification_modern.csv"),
  row.names = FALSE
)
utils::write.csv(locations[locations$type == "grid", ],
                 file.path(diagnostic_dir, "worldclim_offsets_grid266.csv"),
                 row.names = FALSE)
utils::write.csv(locations[locations$type == "station", ],
                 file.path(diagnostic_dir, "worldclim_offsets_city3.csv"),
                 row.names = FALSE)
utils::write.csv(locations[!locations$worldclim_available, ],
                 file.path(diagnostic_dir, "worldclim_masked_locations.csv"),
                 row.names = FALSE)
message("Prepared offsets: ", prepared_file)
message("Direct WorldClim locations: ", sum(locations$worldclim_available), "/269")
message("Uncorrected locations: ", sum(!locations$worldclim_available))
print(locations[locations$type == "station",
                c("site", "WorldClim_1970_2000_C", "LME_1970_2000_C", "d_C")],
      row.names = FALSE)
