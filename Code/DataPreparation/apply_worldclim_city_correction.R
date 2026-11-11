here::i_am("Code/DataPreparation/apply_worldclim_city_correction.R")

# Apply location-specific WorldClim corrections to historical LME city series.
# Retains the original thin-plate spline results and never overwrites them.
# This is an intermediate validation product; downstream models are unchanged.

output_dir <- here::here("Output", "LME_Bias_Correction")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

offsets_path <- Sys.getenv(
  "WORLDCLIM_OFFSETS_FILE",
  here::here("Output", "L1", "L1_worldclim_offsets_all_269.csv")
)
if (!file.exists(offsets_path)) {
  stop("Cannot find WorldClim offsets: ", offsets_path,
       "\nSet WORLDCLIM_OFFSETS_FILE to your actual offsets CSV.")
}

offsets <- read.csv(offsets_path, stringsAsFactors = FALSE)
required_offset_cols <- c("site", "type", "longitude", "latitude", "d_C",
                          "worldclim_available", "correction_method")
stopifnot(all(required_offset_cols %in% names(offsets)))
stopifnot(!anyDuplicated(offsets$site))

grid_offsets <- offsets[offsets$type == "grid", ]
stopifnot(nrow(grid_offsets) == 266L)

historical_grid <- readRDS(here::here(
  "Output", "Intermediate", "LME", "lme_annual_1368_1911.rds"
))
stopifnot(all(c("coordinates", "annual_kelvin") %in% names(historical_grid)))
stopifnot(nrow(historical_grid$coordinates) == 266L)

# Compare grids as coordinate sets, rather than presuming identical row order.
coord_keys <- function(lon, lat) {
  sprintf("%.8f|%.8f", as.numeric(lon), as.numeric(lat))
}
modern_key <- coord_keys(grid_offsets$longitude, grid_offsets$latitude)
historical_key <- coord_keys(historical_grid$coordinates$long,
                             historical_grid$coordinates$lati)
if (!setequal(modern_key, historical_key)) {
  stop("Modern and historical native LME grids differ. Check spatial coordinates.")
}
message("Native grid coordinates verified: 266 locations match.")

city_archive <- readRDS(here::here(
  "Output", "Intermediate", "LME", "lme_city3_annual_1368_1911.rds"
))
required_city_fields <- c("coordinates", "years", "members", "annual_kelvin")
stopifnot(all(required_city_fields %in% names(city_archive)))
stopifnot(identical(dim(city_archive$annual_kelvin),
                    c(3L, 544L, 13L)))
stopifnot(identical(as.integer(city_archive$years), 1368:1911))

cities <- city_archive$coordinates
stopifnot(all(c("city", "long", "lati") %in% names(cities)))
idx <- match(cities$city, offsets$site)
if (anyNA(idx)) stop("One or more city names missing from offsets CSV")
city_offsets <- offsets[idx, , drop = FALSE]
if (!all(city_offsets$worldclim_available) ||
    !all(city_offsets$correction_method == "WorldClim_direct")) {
  stop("All three city offsets must come directly from WorldClim")
}
if (!identical(coord_keys(cities$long, cities$lati),
               coord_keys(city_offsets$longitude, city_offsets$latitude))) {
  stop("City coordinates in historical LME and WorldClim offsets do not match")
}
if (any(!is.finite(city_offsets$d_C))) stop("Nonfinite city corrections")

# The historical city archive already contains thin-plate interpolated values.
# Add one city-specific shift to every year/member without another spline.
uncorrected <- city_archive$annual_kelvin - 273.15
corrected <- sweep(uncorrected, MARGIN = 1,
                   STATS = city_offsets$d_C, FUN = "+")

# Structural checks: the only change is the location-specific additive shift.
for (i in seq_len(nrow(cities))) {
  delta <- as.numeric(corrected[i, , ] - uncorrected[i, , ])
  stopifnot(max(abs(delta - city_offsets$d_C[i])) < 1e-10)
  stopifnot(abs(stats::var(as.numeric(corrected[i, , ])) -
                stats::var(as.numeric(uncorrected[i, , ]))) < 1e-9)
}

out <- expand.grid(
  city = cities$city,
  year = as.integer(city_archive$years),
  member = as.character(city_archive$members),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
out$uncorrected_C <- as.vector(uncorrected)
out$offset_C <- city_offsets$d_C[match(out$city, cities$city)]
out$corrected_C <- as.vector(corrected)
out$longitude <- cities$long[match(out$city, cities$city)]
out$latitude <- cities$lati[match(out$city, cities$city)]
out <- out[c("city", "longitude", "latitude", "member", "year",
             "uncorrected_C", "offset_C", "corrected_C")]

out_file <- file.path(output_dir, "lme_city3_corrected_annual_1368_1911.csv")
write.csv(out, out_file, row.names = FALSE)
write.csv(city_offsets[c("site", "longitude", "latitude", "d_C",
                         "correction_method")],
          file.path(output_dir, "city_correction_diagnostics.csv"),
          row.names = FALSE)
message("Saved corrected city LME series: ", out_file)
print(city_offsets[c("site", "d_C", "correction_method")], row.names = FALSE)
