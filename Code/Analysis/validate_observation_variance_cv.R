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
# file.path() does not reliably retain vector names. Keep the
# filename keys for the explicit paths[["selected"]] lookups below.
paths <- stats::setNames(
  file.path(directory, unname(files)),
  names(files)
)
stopifnot(identical(names(paths), names(files)))
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
# Inspect why a few held-out years may dominate the baseline
# predictive log score when the modeled measurement variance is tiny.
components$variance_without_u <- expected_base
components$standardized_residual_without_u <-
  expected_residual / sqrt(expected_base)
components$log_score_without_u <-
  0.5 * (log(2 * pi * expected_base) +
         expected_residual^2 / expected_base)
components$log_score_with_u <-
  0.5 * (log(2 * pi * (expected_base + u)) +
         expected_residual^2 / (expected_base + u))

fold_diagnostics <- do.call(rbind, lapply(fold_bounds$fold, function(fold) {
  z <- components[components$fold == fold, , drop = FALSE]
  data.frame(
    fold = fold,
    heldout_city_years = nrow(z),
    minimum_baseline_variance = min(z$variance_without_u),
    p01_baseline_variance = unname(
      stats::quantile(z$variance_without_u, 0.01)
    ),
    median_baseline_variance = stats::median(z$variance_without_u),
    n_baseline_variance_below_1e_4 = sum(
      z$variance_without_u < 1e-4
    ),
    maximum_absolute_standardized_residual_at_zero = max(
      abs(z$standardized_residual_without_u)
    ),
    log_score_at_zero = mean(z$log_score_without_u),
    log_score_at_selected = mean(z$log_score_with_u),
    improvement = mean(z$log_score_without_u -
                         z$log_score_with_u)
  )
}))
readr::write_csv(
  fold_diagnostics,
  file.path(directory, "cv_fold_variance_diagnostics.csv")
)

extreme <- components[
  order(components$log_score_without_u, decreasing = TRUE),
  c("fold", "city", "year", "Xstar", "predictive_mean",
    "predictive_residual", "variance_without_u",
    "standardized_residual_without_u", "log_score_without_u",
    "log_score_with_u")
]
extreme <- utils::head(extreme, 30L)
readr::write_csv(
  extreme, file.path(directory, "cv_largest_zero_variance_scores.csv")
)

# The full curve remains in cv_score_curve.png. A second plot focuses
# on the region of the optimum rather than being compressed by the
# large scores at u = 0. Producing it does not rerun the CV fit.
positive_u <- as.numeric(grid$u[grid$u > 0])
zoom_lower <- max(min(positive_u), u / 8, 1e-6)
zoom_upper <- min(max(grid$u), max(u * 5, 0.1))
if (is.finite(zoom_lower) && is.finite(zoom_upper) &&
    zoom_lower < zoom_upper) {
  zoom_data <- grid[
    grid$u >= zoom_lower & grid$u <= zoom_upper, ,
    drop = FALSE
  ]
  if (nrow(zoom_data) > 2L) {
    p <- ggplot2::ggplot(
      zoom_data,
      ggplot2::aes(x = u, y = predictive_log_score)
    ) +
      ggplot2::geom_line(linewidth = 0.8, colour = "black") +
      ggplot2::geom_vline(
        xintercept = u, linetype = "dashed", colour = "red"
      ) +
      ggplot2::labs(
        title = "Held-out predictive log score near the optimum",
        subtitle = paste0("Selected u = ", signif(u, 5),
                          " Celsius squared"),
        x = expression("Additional observation variance " * u *
                         " (" * degree * "C)"^2),
        y = "Mean predictive negative log score"
      ) +
      ggplot2::theme_gray(base_size = 12)
    ggplot2::ggsave(
      file.path(directory, "cv_score_curve_zoom.png"),
      p, width = 7, height = 4.5, dpi = 300, bg = "white"
    )
  }
}

print(summary)
print(fold_diagnostics)
message("Five-fold observation-variance diagnostics verified.")
message(
  "Additional diagnostics saved: cv_fold_variance_diagnostics.csv, ",
  "cv_largest_zero_variance_scores.csv, cv_score_curve_zoom.png"
)
