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
  message("Modern/historical native LME coordinates match (266 locations).")
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
utils::write.csv(cbind(locations, wc_monthly),
                 file.path(diagnostic_dir, "worldclim_monthly_all_269.csv"),
                 row.names = FALSE)

# Same thin-plate regression-spline operator used for historical city LME.
train_coords <- data.frame(long = native_grid$longitude,
                           lati = native_grid$latitude)
smooth <- mgcv::smoothCon(
  mgcv::s(long, lati, bs = "tp", k = 266),
  data = train_coords, absorb.cons = TRUE, scale.penalty = TRUE
)[[1]]
train_design <- cbind(intercept = 1, smooth$X)
target_design <- cbind(
  intercept = 1,
  mgcv::PredictMat(smooth, data.frame(long = locations$longitude,
                                      lati = locations$latitude))
)
penalty <- matrix(0, ncol(train_design), ncol(train_design))
penalty[-1, -1] <- smooth$S[[1]]
fixed_sp <- 1e-8
normal <- crossprod(train_design) + fixed_sp * penalty
if (!is.finite(rcond(normal)) || rcond(normal) < 1e-14) {
  stop("Thin-plate interpolation operator is numerically singular.")
}
operator <- target_design %*% solve(normal, t(train_design))

# Modern CESM-LME annual values: all 13 members, each year 1970-2000.
modern_records <- vector("list", 13L)
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

# Small prepared input intended for version control.
prepared_file <- file.path(prepared_dir, "worldclim_offsets_all_269.csv")
utils::write.csv(locations, prepared_file, row.names = FALSE)
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
