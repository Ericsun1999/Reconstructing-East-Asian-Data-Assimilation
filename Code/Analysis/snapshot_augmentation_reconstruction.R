here::i_am("Code/Analysis/snapshot_augmentation_reconstruction.R")

# Archive a verified generated reconstruction BEFORE switching to the
# other documentary input. No model is estimated or altered here.
#
# Rscript --vanilla Code/Analysis/snapshot_augmentation_reconstruction.R original city
# Rscript --vanilla Code/Analysis/snapshot_augmentation_reconstruction.R augmented city
# Replace city with grid only after the corresponding spatial assimilation
# has completed; a grid snapshot also includes all city diagnostics.

library(here)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L ||
    !args[1L] %in% c("original", "augmented") ||
    !args[2L] %in% c("city", "grid")) {
  stop("Usage: snapshot_augmentation_reconstruction.R ",
       "<original|augmented> <city|grid>")
}
input_mode <- args[1L]
stage <- args[2L]
root <- here::here()
out_root <- file.path(root, "Output", "Intermediate")
reaches_dir <- file.path(out_root, "REACHES")
city_dir <- file.path(out_root, "Assimilation")
source_file <- if (input_mode == "original") {
  here::here("Data", "temperature index value.v1.xlsx")
} else {
  file.path(reaches_dir, "temperature_site_year_augmented.csv")
}
calibration_file <- file.path(out_root, "calibration_parameters.rds")
kriging_meta_file <- file.path(
  reaches_dir, "reaches_kriging_metadata.csv"
)
city_meta_file <- file.path(city_dir, "assimilation_metadata.csv")
city_names <- c("Beijing", "Shanghai", "HongKong")
city_files <- file.path(
  city_dir, paste0("assimilation_", city_names, ".csv")
)
inputs <- c(
  calibration_file,
  kriging_meta_file,
  file.path(reaches_dir, "reaches_kriging_city3_mean.csv"),
  file.path(reaches_dir, "reaches_kriging_city3_sd.csv"),
  city_files, city_meta_file,
  file.path(city_dir, "assimilation_metrics.csv")
)
if (stage == "grid") {
  inputs <- c(
    inputs,
    file.path(city_dir, paste0(
      "assimilated_posterior_lme_grid_",
      c("mean.csv", "variance.csv", "metadata.csv")
    )),
    file.path(reaches_dir, c(
      "reaches_kriging_lme_grid_mean.csv",
      "reaches_kriging_lme_grid_variance.csv"
    ))
  )
}
missing <- unique(c(
  source_file, inputs
))[!file.exists(unique(c(source_file, inputs)))]
if (length(missing) > 0L) {
  stop("Missing required generated inputs:\n",
       paste(missing, collapse = "\n"),
       "\nRun the corresponding model stage before snapshotting.")
}

calibration <- readRDS(calibration_file)
kriging_meta <- readr::read_csv(
  kriging_meta_file, show_col_types = FALSE
)
city_meta <- readr::read_csv(
  city_meta_file, show_col_types = FALSE
)
required_kriging <- c(
  "temperature_input_mode", "temperature_source_md5",
  "calibration_file_md5", "process_variance"
)
required_city <- c(
  "city", "measurement_projection_method",
  "extra_observation_variance", "measurement_mc_size",
  "process_variance"
)
if (nrow(kriging_meta) != 1L ||
    !all(required_kriging %in% names(kriging_meta)) ||
    !all(required_city %in% names(city_meta)) ||
    nrow(city_meta) != length(city_names)) {
  stop("Outdated or incomplete model metadata. Regenerate the relevant ",
       "kriging and city-assimilation products.")
}
source_md5 <- unname(tools::md5sum(source_file))
calibration_md5 <- unname(tools::md5sum(calibration_file))
if (!identical(calibration$temperature_input_mode, input_mode) ||
    !identical(as.character(kriging_meta$temperature_input_mode[1L]),
               input_mode) ||
    !identical(as.character(kriging_meta$temperature_source_md5[1L]),
               source_md5) ||
    !identical(as.character(kriging_meta$calibration_file_md5[1L]),
               calibration_md5) ||
    !identical(as.character(calibration$temperature_source_md5),
               source_md5)) {
  stop("The calibration, kriging, and source documentary input do not ",
       "match. Do not compare cross-version outputs.")
}
if (!setequal(as.character(city_meta$city), city_names) ||
    length(unique(city_meta$measurement_projection_method)) != 1L ||
    anyNA(city_meta$measurement_projection_method) ||
    length(unique(city_meta$extra_observation_variance)) != 1L ||
    any(!is.finite(city_meta$extra_observation_variance)) ||
    length(unique(city_meta$measurement_mc_size)) != 1L ||
    any(!is.finite(city_meta$measurement_mc_size)) ||
    any(abs(city_meta$process_variance -
              kriging_meta$process_variance[1L]) > 1e-8)) {
  stop("City metadata disagree on method, variance, Monte Carlo size, ",
       "or the associated kriging calibration.")
}

for (f in city_files) {
  z <- readr::read_csv(f, show_col_types = FALSE)
  columns <- c(
    "year", "smoothed_temperature_celsius", "smoothed_variance",
    "reaches_index_mean", "reaches_mspe"
  )
  if (!all(columns %in% names(z)) ||
      nrow(z) != 544L ||
      !identical(as.integer(z$year), 1368:1911) ||
      any(!is.finite(z$smoothed_temperature_celsius)) ||
      any(!is.finite(z$smoothed_variance)) ||
      any(z$smoothed_variance < 0)) {
    stop("Incomplete or invalid annual city reconstruction: ", f)
  }
}

if (stage == "grid") {
  if ("kriging_scope" %in% names(kriging_meta) &&
      identical(as.character(kriging_meta$kriging_scope[1L]),
                "city_only")) {
    stop("A city-only kriging run cannot be snapshotted as a spatial grid.")
  }
  grid_meta <- readr::read_csv(
    file.path(
      city_dir, "assimilated_posterior_lme_grid_metadata.csv"
    ), show_col_types = FALSE
  )
  if (nrow(grid_meta) != 1L ||
      !all(c(
        "process_variance", "measurement_mc_size",
        "measurement_projection_method", "extra_observation_variance",
        "reaches_mean_md5", "reaches_variance_md5",
        "reaches_metadata_md5"
      ) %in% names(grid_meta)) ||
      !isTRUE(all.equal(
        as.numeric(grid_meta$process_variance[1L]),
        as.numeric(kriging_meta$process_variance[1L]),
        tolerance = 1e-8
      )) ||
      !isTRUE(all.equal(
        as.numeric(grid_meta$extra_observation_variance[1L]),
        as.numeric(city_meta$extra_observation_variance[1L]),
        tolerance = 1e-8
      )) ||
      !identical(
        as.character(grid_meta$measurement_projection_method[1L]),
        as.character(city_meta$measurement_projection_method[1L])
      ) ||
      as.integer(grid_meta$measurement_mc_size[1L]) !=
        as.integer(city_meta$measurement_mc_size[1L]) ||
      !identical(
        as.character(grid_meta$reaches_metadata_md5[1L]),
        unname(tools::md5sum(kriging_meta_file))
      ) ||
      !identical(
        as.character(grid_meta$reaches_mean_md5[1L]),
        unname(tools::md5sum(file.path(
          reaches_dir, "reaches_kriging_lme_grid_mean.csv"
        )))
      ) ||
      !identical(
        as.character(grid_meta$reaches_variance_md5[1L]),
        unname(tools::md5sum(file.path(
          reaches_dir, "reaches_kriging_lme_grid_variance.csv"
        )))
      )) {
    stop("Spatial assimilation metadata do not match city/kriging run.")
  }
}

destination <- file.path(
  out_root, "AugmentationSensitivity", input_mode, stage
)
if (dir.exists(destination)) {
  stop("A snapshot already exists at ", destination,
       ". Preserve it; manually archive or remove it before ",
       "capturing a new comparable reconstruction.")
}
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
if (!all(file.copy(inputs, destination, overwrite = FALSE))) {
  stop("Could not archive all requested reconstruction files.")
}

manifest <- data.frame(
  input_mode = input_mode,
  stage = stage,
  source_md5 = source_md5,
  calibration_md5 = calibration_md5,
  kriging_metadata_md5 = unname(tools::md5sum(kriging_meta_file)),
  projection_method =
    as.character(city_meta$measurement_projection_method[1L]),
  extra_observation_variance =
    as.numeric(city_meta$extra_observation_variance[1L]),
  measurement_mc_size = as.integer(city_meta$measurement_mc_size[1L]),
  copied_files = length(inputs),
  captured_at = as.character(Sys.time())
)
readr::write_csv(
  manifest, file.path(destination, "snapshot_manifest.csv")
)
print(manifest)
message("Verified ", input_mode, " snapshot saved to ", destination)
