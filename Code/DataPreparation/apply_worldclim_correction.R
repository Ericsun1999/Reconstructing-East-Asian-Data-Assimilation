here::i_am("Code/DataPreparation/apply_worldclim_correction.R")

# Apply WorldClim-derived additive offsets AFTER the existing historical
# LME thin-plate spline. Original files are left unchanged.
# No IDW/ocean-grid corrections are applied by this city-only module.

original_dir <- here::here("Output", "Intermediate", "LME")
corrected_dir <- here::here("Output", "Intermediate", "LME", "WorldClim_Corrected")
dir.create(corrected_dir, recursive = TRUE, showWarnings = FALSE)

city_rds <- file.path(original_dir, "lme_city3_annual_1368_1911.rds")
if (!file.exists(city_rds)) {
  stop("Missing original spline-interpolated city LME archive: ", city_rds,
       "\nRun Code/DataPreparation/prepare_lme_annual.R first.")
}

# Default input is the small, version-controlled three-station offset table.
# Use WORLDCLIM_OFFSETS_FILE to test a newly regenerated offset file.
offset_override <- Sys.getenv("WORLDCLIM_OFFSETS_FILE", unset = "")
offset_file <- if (nzchar(offset_override)) {
  path.expand(offset_override)
} else {
  offset_override <- Sys.getenv(
  "WORLDCLIM_OFFSETS_FILE",
  unset = ""
)

offset_file <- if (nzchar(offset_override)) {
  path.expand(offset_override)
} else {
  here::here(
    "Output",
    "LME_Bias_Correction",
    "worldclim_offsets_city3.csv"
  )
}
}
if (!file.exists(offset_file)) {
  stop("Missing WorldClim offset file: ", offset_file)
}

archive <- readRDS(city_rds)
required_archive <- c("coordinates", "years", "members", "units", "annual_kelvin")
if (!all(required_archive %in% names(archive))) {
  stop("City LME archive missing: ",
       paste(setdiff(required_archive, names(archive)), collapse = ", "))
}
if (!identical(as.character(archive$units), "Kelvin")) {
  stop("Original city archive must contain Kelvin values.")
}

cities <- as.data.frame(archive$coordinates)
if (!all(c("city", "long", "lati") %in% names(cities))) {
  stop("City archive coordinates require city, long, and lati.")
}
if (nrow(cities) != 3L || anyDuplicated(cities$city)) {
  stop("Expected exactly three unique validation stations.")
}
years <- as.integer(archive$years)
members <- as.character(archive$members)
if (!identical(years, 1368:1911) || length(members) != 13L ||
    anyDuplicated(members)) {
  stop("Unexpected years or members in the historical LME archive.")
}
original_kelvin <- archive$annual_kelvin
if (!identical(as.integer(dim(original_kelvin)),
               c(nrow(cities), length(years), length(members))) ||
    any(!is.finite(original_kelvin))) {
  stop("Malformed/nonfinite Kelvin array in original city archive.")
}

offsets <- utils::read.csv(offset_file, stringsAsFactors = FALSE)
required_offsets <- c("site", "longitude", "latitude", "d_C",
                      "worldclim_available", "correction_method")
if (!all(required_offsets %in% names(offsets))) {
  stop("Offset CSV missing: ",
       paste(setdiff(required_offsets, names(offsets)), collapse = ", "))
}
if (anyDuplicated(offsets$site)) stop("Duplicate sites in offsets file.")
i <- match(as.character(cities$city), offsets$site)
if (anyNA(i)) stop("Offsets missing one or more validation stations.")
selected <- offsets[i, , drop = FALSE]
if (!all(as.character(selected$correction_method) == "WorldClim_direct") ||
    !all(tolower(as.character(selected$worldclim_available)) == "true")) {
  stop("All city offsets must be directly observed WorldClim values, not IDW.")
}
if (any(!is.finite(selected$d_C))) stop("Non-finite city offsets.")
if (any(abs(as.numeric(cities$long) - selected$longitude) > 1e-6) ||
    any(abs(as.numeric(cities$lati) - selected$latitude) > 1e-6)) {
  stop("WorldClim and LME city coordinates differ; do not apply offsets.")
}

# Add each city's d(s) exactly once. A °C temperature DIFFERENCE has the
# same numeric size in kelvin; preserving the Kelvin downstream interface.
corrected_kelvin <- sweep(original_kelvin, 1L, selected$d_C, "+")

# Validate a pure translation for each complete member/year panel.
summary_rows <- vector("list", nrow(cities))
for (j in seq_len(nrow(cities))) {
  before <- as.numeric(original_kelvin[j, , ])
  after <- as.numeric(corrected_kelvin[j, , ])
  max_error <- max(abs((after - before) - selected$d_C[j]))
  variance_error <- abs(stats::var(after) - stats::var(before))
  if (!is.finite(max_error) || max_error > 1e-10 ||
      !is.finite(variance_error) || variance_error > 1e-9) {
    stop("Translation check failed for ", cities$city[j])
  }
  summary_rows[[j]] <- data.frame(
    city = as.character(cities$city[j]),
    offset_C = selected$d_C[j],
    uncorrected_mean_C = mean(before) - 273.15,
    corrected_mean_C = mean(after) - 273.15,
    max_addition_error = max_error,
    variance_difference = variance_error,
    correction_method = as.character(selected$correction_method[j]),
    stringsAsFactors = FALSE
  )
}

# Match the ORIGINAL file's exact columns and Kelvin unit; no downstream
# change to any 273.15 conversion is necessary.
out <- expand.grid(
  city = as.character(cities$city),
  year = years,
  member = members,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
city_pos <- match(out$city, cities$city)
out$long <- as.numeric(cities$long)[city_pos]
out$lati <- as.numeric(cities$lati)[city_pos]
out$temperature_kelvin <- as.vector(corrected_kelvin)
out <- out[c("city", "long", "lati", "member", "year", "temperature_kelvin")]
if (nrow(out) != 3L * length(years) * length(members) ||
    anyDuplicated(out[c("city", "member", "year")]) ||
    any(!is.finite(out$temperature_kelvin))) {
  stop("Invalid corrected city CSV.")
}

corrected_archive <- archive
corrected_archive$annual_kelvin <- corrected_kelvin
corrected_archive$bias_correction <- list(
  type = "additive_climatological_level_adjustment",
  reference = "WorldClim v2.1 monthly tavg, 1970-2000",
  ensemble_reference = "13 CESM-LME members, 1970-2000",
  offset_source = normalizePath(offset_file, mustWork = TRUE),
  cities = as.character(cities$city),
  offsets_C = as.numeric(selected$d_C),
  method = "direct WorldClim, no IDW",
  interpolation = "existing historical thin-plate spline; not repeated"
)

csv_out <- file.path(corrected_dir, "lme_city3_annual_1368_1911.csv")
rds_out <- file.path(corrected_dir, "lme_city3_annual_1368_1911.rds")
diag_out <- file.path(corrected_dir, "worldclim_correction_diagnostics.csv")
utils::write.csv(out, csv_out, row.names = FALSE)
saveRDS(corrected_archive, rds_out, compress = "xz")
utils::write.csv(do.call(rbind, summary_rows), diag_out, row.names = FALSE)
message("Created corrected, Kelvin-compatible three-city LME:")
message("  ", csv_out)
message("  ", rds_out)
message("  ", diag_out)
print(do.call(rbind, summary_rows), row.names = FALSE)
message("Original city LME files were NOT modified. No grid/ocean IDW was applied.")
