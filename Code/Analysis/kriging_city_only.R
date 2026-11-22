here::i_am("Code/Analysis/kriging_city_only.R")

# Fast three-city sensitivity run. Uses the identical interval-censored
# kriging functions and calibrated covariance from Get_tempe_all_data.R
# without computing the 53x49 or 266-location spatial grids.
# The main full-grid pipeline remains unchanged.

library(here)
library(sp)
library(readr)

calibration_file <- here::here(
  "Output", "Intermediate", "calibration_parameters.rds"
)
output_dir <- here::here("Output", "Intermediate", "REACHES")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(calibration_file)) {
  stop("Missing calibration_parameters.rds; run prepare_calibration.R.")
}
calibration <- readRDS(calibration_file)
required <- c("var_fit2", "vario_fit2", "temp2",
              "temperature_input_mode", "temperature_source_md5",
              "temperature_source_file")
if (!all(required %in% names(calibration)) ||
    !calibration$temperature_input_mode %in% c("original", "augmented")) {
  stop("Calibration RDS is incomplete or has an invalid input mode.")
}
input_file <- calibration$temperature_source_file
source_md5 <- unname(tools::md5sum(input_file))
if (!file.exists(input_file) ||
    !identical(as.character(calibration$temperature_source_md5),
               source_md5)) {
  stop("The original documentary input does not match this calibration.")
}
sigma_Y2 <- as.numeric(calibration$vario_fit2$psill2)
sigma_E2 <- as.numeric(calibration$vario_fit2$psill1)
alpha_km <- as.numeric(calibration$vario_fit2$range)
if (any(!is.finite(c(sigma_Y2, sigma_E2, alpha_km))) ||
    sigma_Y2 <= 0 || sigma_E2 < 0 || alpha_km <= 0) {
  stop("The calibrated covariance parameters are invalid.")
}
y2 <- calibration$var_fit2$y2
if (!methods::is(y2, "SpatialPointsDataFrame") ||
    !all(c("level", "year") %in% names(y2))) {
  stop("The fitted observation data are not a usable spatial dataset.")
}
years <- sort(unique(as.integer(calibration$temp2$year)))
if (!identical(years, sort(unique(as.integer(y2$year)))) ||
    anyNA(years) || any(!years %in% 1368:1911)) {
  stop("Inconsistent calibrated historical years.")
}

# Reuse only the actual definitions, not the heavyweight script body.
load_functions <- function(path, required_names) {
  exprs <- parse(file = path)
  loaded <- character()
  for (x in exprs) {
    if (is.call(x) &&
        identical(x[[1L]], as.name("<-")) &&
        is.symbol(x[[2L]]) &&
        is.call(x[[3L]]) &&
        identical(x[[3L]][[1L]], as.name("function"))) {
      name <- as.character(x[[2L]])
      if (name %in% required_names) {
        assign(name, eval(x[[3L]], envir = .GlobalEnv),
               envir = .GlobalEnv)
        loaded <- c(loaded, name)
      }
    }
  }
  if (!setequal(required_names, loaded)) {
    stop("Missing expected kriging functions: ",
         paste(setdiff(required_names, loaded), collapse = ", "))
  }
}
load_functions(
  here::here("Code", "Get_tempe_all_data.R"),
  c("Ez_discrete", "Varz_discrete", "EZstar_h",
    "ordinal_covariance", "build_covariance_lookup",
    "solve_symmetric_system", "predict_grid_one_year")
)

city_locations <- data.frame(
  city = c("HongKong", "Shanghai", "Beijing"),
  long = c(114.167, 121.433, 116.283),
  lat = c(22.333, 31.167, 39.933)
)
coords <- as.matrix(city_locations[c("long", "lat")])
total_variance <- sigma_Y2 + sigma_E2
sigma <- sqrt(total_variance)
cuts <- c(-Inf, -1.5, -0.5, 0.5, Inf)
vals <- c(-2, -1, 0, 1)
mean_z <- Ez_discrete(sigma, cuts, vals)
var_z <- Varz_discrete(sigma, cuts, vals, mean_z)
ezstar_h <- EZstar_h(sigma, cuts, vals)
rho_lookup_size <- 4001L
message("Building the original ordinal covariance lookup ...")
lookup <- build_covariance_lookup(
  sigma_Y2 = sigma_Y2, sigma_E2 = sigma_E2, sigma = sigma,
  number_of_points = rho_lookup_size
)

n <- length(years)
mean_matrix <- matrix(
  NA_real_, nrow = 3L, ncol = n,
  dimnames = list(city_locations$city, paste0("x", years))
)
variance_matrix <- mean_matrix
diagnostics <- data.frame(
  year = years, number_of_observations = NA_integer_,
  numerical_jitter = NA_real_
)
start_time <- Sys.time()
for (i in seq_along(years)) {
  pred <- predict_grid_one_year(
    current_year = years[i], observation_data = y2,
    prediction_coordinates = coords,
    sigma_Y2 = sigma_Y2, sigma_E2 = sigma_E2,
    alpha_km = alpha_km, total_variance = total_variance,
    mean_z = mean_z, variance_z = var_z, ezstar_h = ezstar_h,
    covariance_lookup_function = lookup$function_object
  )
  mean_matrix[, i] <- as.numeric(pred$mean)
  variance_matrix[, i] <- as.numeric(pred$variance)
  diagnostics$number_of_observations[i] <-
    as.integer(pred$number_of_observations)
  diagnostics$numerical_jitter[i] <- as.numeric(pred$jitter)
  if (i %% 50L == 0L || i == n) {
    message("City kriging ", i, "/", n, " years")
  }
}
if (any(!is.finite(mean_matrix)) ||
    any(!is.finite(variance_matrix)) ||
    any(variance_matrix < 0)) {
  stop("The three-city kriging results have invalid values.")
}
elapsed <- as.numeric(difftime(
  Sys.time(), start_time, units = "mins"
))
city_mean <- cbind(city_locations, as.data.frame(
  mean_matrix, check.names = FALSE
))
city_sd <- cbind(city_locations, as.data.frame(
  sqrt(pmax(variance_matrix, 0)), check.names = FALSE
))
readr::write_csv(
  city_mean, file.path(output_dir, "reaches_kriging_city3_mean.csv")
)
readr::write_csv(
  city_sd, file.path(output_dir, "reaches_kriging_city3_sd.csv")
)
readr::write_csv(
  transform(
    city_locations,
    coordinate_source = "GHCN nominal station coordinate",
    prediction_method = "direct interval-censored kriging",
    number_of_event_years = n, elapsed_minutes = elapsed
  ), file.path(output_dir, "reaches_kriging_city3_metadata.csv")
)

metadata <- data.frame(
  temperature_input_mode = calibration$temperature_input_mode,
  temperature_source_md5 = source_md5,
  calibration_file_md5 = unname(tools::md5sum(calibration_file)),
  kriging_scope = "city_only",
  number_of_event_years = n,
  number_of_grid_locations = NA_integer_,
  number_of_lme_grid_locations = NA_integer_,
  number_of_city_locations = nrow(city_locations),
  city_coordinate_source = "GHCN nominal station coordinates",
  city_prediction_method = "direct interval-censored kriging",
  process_variance = sigma_Y2, nugget_variance = sigma_E2,
  total_variance = total_variance, range_km = alpha_km,
  covariance_lookup_size = rho_lookup_size,
  elapsed_minutes = elapsed,
  city_kriging_elapsed_minutes = elapsed
)
readr::write_csv(
  metadata, file.path(output_dir, "reaches_kriging_metadata.csv")
)
readr::write_csv(
  diagnostics, file.path(output_dir, "reaches_kriging_city3_year_diagnostics.csv")
)
message("Completed three-city-only kriging. Spatial-grid results were ",
        "not recomputed; rerun Get_tempe_all_data.R before grid assimilation.")
