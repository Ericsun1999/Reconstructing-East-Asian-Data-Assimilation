# Shared measurement-discrepancy variance utilities.
# The fitted addition is in squared degrees Celsius and is added once
# to the Monte Carlo measurement residual variance (not to the
# Kalman smoother variance after filtering).

read_extra_observation_variance <- function(
    reaches_metadata_file,
    variance_file = Sys.getenv("REACHES_EXTRA_VARIANCE_FILE", unset = ""),
    projection_method = "empirical_reference_joint_mc_v1") {

  variance_file <- path.expand(trimws(variance_file))
  # An explicitly fixed common u is useful for documentary-input
  # sensitivity comparisons. It is NOT an independently fitted u
  # for the augmented dataset, and is labeled accordingly.
  fixed <- trimws(Sys.getenv(
    "REACHES_EXTRA_VARIANCE_FIXED", unset = ""
  ))
  if (nzchar(fixed)) {
    if (nzchar(variance_file)) {
      stop("Set either REACHES_EXTRA_VARIANCE_FIXED or ",
           "REACHES_EXTRA_VARIANCE_FILE, not both.")
    }
    value <- suppressWarnings(as.numeric(fixed))
    if (length(value) != 1L || !is.finite(value) || value < 0) {
      stop("REACHES_EXTRA_VARIANCE_FIXED must be a nonnegative number.")
    }
    return(list(
      u = value, file = NA_character_, file_md5 = NA_character_,
      mode = "fixed_input_sensitivity"
    ))
  }
  if (!nzchar(variance_file)) {
    return(list(
      u = 0,
      file = NA_character_,
      file_md5 = NA_character_,
      mode = "legacy_zero"
    ))
  }
  if (!file.exists(variance_file)) {
    stop("Requested measurement variance file does not exist: ",
         variance_file)
  }
  if (!file.exists(reaches_metadata_file)) {
    stop("REACHES kriging metadata is missing: ", reaches_metadata_file)
  }

  fitted <- readr::read_csv(variance_file, show_col_types = FALSE)
  required <- c(
    "selected_u", "temperature_input_mode", "temperature_source_md5",
    "number_of_folds", "variance_units", "measurement_projection_method"
  )
  if (nrow(fitted) != 1L || !all(required %in% names(fitted))) {
    stop("Measurement variance file must contain a single valid fit: ",
         variance_file)
  }

  kriging <- readr::read_csv(reaches_metadata_file, show_col_types = FALSE)
  kriging_required <- c("temperature_input_mode",
                        "temperature_source_md5")
  if (nrow(kriging) != 1L ||
      !all(kriging_required %in% names(kriging))) {
    stop("Kriging metadata lacks source identity. Rerun ",
         "Code/Get_tempe_all_data.R before applying the extra variance.")
  }

  u <- as.numeric(fitted$selected_u[1L])
  folds <- as.integer(fitted$number_of_folds[1L])
  if (length(u) != 1L || !is.finite(u) || u < 0 ||
      length(folds) != 1L || is.na(folds) || folds != 5L ||
      !identical(as.character(fitted$variance_units[1L]), "celsius_squared")) {
    stop("Selected discrepancy variance is invalid or incomplete.")
  }
  fitted_method <- as.character(fitted$measurement_projection_method[1L])
  if (length(fitted_method) != 1L || is.na(fitted_method) ||
      !identical(fitted_method, projection_method)) {
    stop(
      "The variance file was fitted with a different measurement ",
      "projection method (or predates versioned projection fitting). ",
      "Use REACHES_MEASUREMENT_PROJECTION_METHOD=",
      fitted_method, " and rerun blocked CV if needed."
    )
  }
  if (!identical(as.character(fitted$temperature_input_mode[1L]),
                 as.character(kriging$temperature_input_mode[1L])) ||
      !identical(as.character(fitted$temperature_source_md5[1L]),
                 as.character(kriging$temperature_source_md5[1L]))) {
    stop(
      "The selected discrepancy variance was fitted to different ",
      "REACHES input data. Rerun blocked CV for this input version."
    )
  }

  list(
    u = u,
    file = normalizePath(variance_file),
    file_md5 = unname(tools::md5sum(variance_file)),
    mode = "blocked_cv"
  )
}

annual_unconditional_variance <- function(M, r2) {
  if (length(M) != length(r2) - 1L || length(r2) < 2L ||
      any(!is.finite(c(M, r2))) || any(r2 <= 0)) {
    stop("Invalid annual state-space prior parameters.")
  }
  V <- numeric(length(r2))
  V[1L] <- r2[1L]
  for (t in 2:length(r2)) {
    V[t] <- M[t - 1L]^2 * V[t - 1L] + r2[t]
  }
  if (any(!is.finite(V)) || any(V <= 0)) {
    stop("Unconditional state variance is non-positive or non-finite.")
  }
  V
}

predictive_log_score <- function(residual, base_variance, u) {
  if (length(residual) != length(base_variance) ||
      length(residual) < 1L ||
      any(!is.finite(residual)) ||
      any(!is.finite(base_variance)) ||
      any(base_variance <= 0) ||
      length(u) != 1L || !is.finite(u) || u < 0) {
    stop("Invalid held-out score inputs.")
  }
  S <- base_variance + u
  mean(0.5 * (log(2 * pi * S) + residual^2 / S))
}
