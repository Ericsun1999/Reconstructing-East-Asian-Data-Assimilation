here::i_am("Code/Analysis/validate_observation_variance_cv.R")

# Independent arithmetic and partition checks for a completed five-fold
# predictive-log-score fit. This script does not rerun model fitting.
library(here)
library(readr)

source(here::here("Code", "Analysis", "extra_observation_variance.R"))

mode <- Sys.getenv("REACHES_TEMPERATURE_INPUT", unset = "original")
if (!mode %in% c("original", "augmented")) {
  stop("REACHES_TEMPERATURE_INPUT must be original or augmented.")
}
directory <- here::here(
  "Output", "Intermediate", "ObservationVariance", mode
)
files <- c(
  selected = "selected_u.csv",
  components = "heldout_score_components.csv",
  covariance = "fold_covariance.csv",
  grid = "cv_score_curve.csv"
)
paths <- file.path(directory, files)
if (any(!file.exists(paths))) {
  stop("Missing completed five-fold CV outputs: ",
       paste(paths[!file.exists(paths)], collapse = ", "))
}

selected <- readr::read_csv(paths[["selected"]], show_col_types = FALSE)
components <- readr::read_csv(paths[["components"]], show_col_types = FALSE)
covariance <- readr::read_csv(paths[["covariance"]], show_col_types = FALSE)
grid <- readr::read_csv(paths[["grid"]], show_col_types = FALSE)

required_components <- c(
  "fold", "city", "year", "Xstar", "mu_prior", "V_prior",
  "alpha", "beta", "vdelta", "predictive_mean",
  "baseline_predictive_variance", "predictive_residual"
)
if (nrow(selected) != 1L ||
    !all(required_components %in% names(components)) ||
    !all(c("u", "predictive_log_score") %in% names(grid))) {
  stop("CV outputs are incomplete or malformed.")
}
if (nrow(components) == 0L || nrow(covariance) != 5L) {
  stop("Five-fold CV has incomplete observations or fold metadata.")
}

fold_bounds <- data.frame(
  fold = sprintf("fold_%d", 1:5),
  start = c(1368L, 1501L, 1601L, 1701L, 1801L),
  end = c(1500L, 1600L, 1700L, 1800L, 1911L)
)
if (!setequal(as.character(covariance$fold), fold_bounds$fold) ||
    anyDuplicated(covariance$fold) ||
    !setequal(as.character(components$fold), fold_bounds$fold)) {
  stop("Unexpected, repeated, or absent historical blocks.")
}
index <- match(as.character(components$fold), fold_bounds$fold)
if (anyNA(index) ||
    any(components$year < fold_bounds$start[index]) ||
    any(components$year > fold_bounds$end[index]) ||
    anyDuplicated(components[c("year", "city")])) {
  stop("Duplicate or incorrectly partitioned held-out city/year records.")
}
if (!setequal(as.character(components$city),
              c("HongKong", "Shanghai", "Beijing"))) {
  stop("Unexpected set of city scores.")
}

numeric_columns <- setdiff(required_components, c("fold", "city"))
if (any(!is.finite(as.matrix(components[numeric_columns]))) ||
    any(components$V_prior <= 0) ||
    any(components$vdelta <= 0)) {
  stop("Non-finite or non-positive predictive components.")
}
close_enough <- function(a, b, scale = 1) {
  isTRUE(all.equal(
    as.numeric(a), as.numeric(b),
    tolerance = 1e-8, scale = scale,
    check.attributes = FALSE
  ))
}
expected_mean <- with(
  components, alpha + beta * mu_prior
)
expected_base <- with(
  components, beta^2 * V_prior + vdelta
)
expected_residual <- components$Xstar - expected_mean
if (!close_enough(expected_mean, components$predictive_mean) ||
    !close_enough(expected_base,
                  components$baseline_predictive_variance) ||
    !close_enough(expected_residual,
                  components$predictive_residual)) {
  stop("A saved predictive-score component fails recomputation.")
}

input_file <- if (mode == "original") {
  here::here("Data", "temperature index value.v1.xlsx")
} else {
  here::here(
    "Output", "Intermediate", "REACHES",
    "temperature_site_year_augmented.csv"
  )
}
if (!file.exists(input_file)) stop("Missing documentary input.")
current_md5 <- unname(tools::md5sum(input_file))
if (!identical(as.character(selected$temperature_input_mode[1L]), mode) ||
    !identical(as.character(selected$temperature_source_md5[1L]),
               current_md5)) {
  stop("CV estimate is not associated with the current input dataset.")
}

u <- as.numeric(selected$selected_u[1L])
if (!is.finite(u) || u < 0 ||
    selected$number_of_folds[1L] != 5L ||
    selected$number_of_heldout_city_years[1L] != nrow(components)) {
  stop("The selected variance or sample size is inconsistent.")
}
loss <- function(v) predictive_log_score(
  expected_residual, expected_base, v
)
at_zero <- loss(0)
at_fit <- loss(u)
if (!close_enough(at_zero, selected$cv_score_at_zero[1L]) ||
    !close_enough(at_fit, selected$cv_score_selected[1L]) ||
    at_fit > at_zero + 1e-8 ||
    !any(abs(as.numeric(grid$u) - u) <= 1e-10) ||
    abs(min(grid$predictive_log_score) - at_fit) > 1e-7) {
  stop("Selected CV variance or objective is inconsistent.")
}

summary <- data.frame(
  temperature_input_mode = mode,
  input_md5 = current_md5,
  heldout_city_years = nrow(components),
  folds = length(unique(components$fold)),
  selected_extra_variance_celsius_squared = u,
  score_at_zero = at_zero,
  score_at_selected = at_fit,
  score_improvement = at_zero - at_fit,
  checked_at = as.character(Sys.time())
)
readr::write_csv(
  summary, file.path(directory, "cv_validation_summary.csv")
)
print(summary)
message("Five-fold observation-variance diagnostics verified.")
