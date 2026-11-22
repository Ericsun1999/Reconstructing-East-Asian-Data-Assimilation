here::i_am("Code/SupportingAnalyses/validate_observation_variance_cv.R")

# Independent arithmetic and partition checks for a completed five-fold
# predictive-log-score fit. This script does not rerun model fitting.
library(here)
library(readr)

# Import only the reusable Figure 9(d) mathematical definitions.
# Do not source Figure9d.R: that would execute the full figure pipeline.
required_figure9_helpers <- c("predictive_log_score")
figure9_expressions <- parse(file = here::here("Code", "Figure9d.R"))
figure9_helpers <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
    is.symbol(expr[[2L]]) &&
    as.character(expr[[2L]]) %in% required_figure9_helpers &&
    is.call(expr[[3L]]) &&
    identical(expr[[3L]][[1L]], as.name("function"))
}, as.list(figure9_expressions))
found_figure9_helpers <- vapply(
  figure9_helpers, function(expr) as.character(expr[[2L]]), character(1L)
)
if (!setequal(found_figure9_helpers, required_figure9_helpers) ||
    anyDuplicated(found_figure9_helpers)) {
  stop("Missing or repeated mathematical definitions in Figure9d.R.")
}
for (expr in figure9_helpers) eval(expr, envir = .GlobalEnv)
rm(figure9_expressions, figure9_helpers, found_figure9_helpers,
   required_figure9_helpers)

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

# CV results fitted using the legacy reference-mean projection must
# not be silently reused after correcting the zero-mean model.
required_projection <- "empirical_reference_zero_centered_joint_mc_v2"
if (!"measurement_projection_method" %in% names(selected) ||
    nrow(selected) != 1L ||
    is.na(selected$measurement_projection_method[1L]) ||
    !identical(as.character(selected$measurement_projection_method[1L]),
               required_projection)) {
  stop(
    "The selected_u.csv file predates the zero-centered measurement ",
    "projection. Rerun Code/Analysis/cv_observation_variance.R to ",
    "refit u with the corrected observation model."
  )
}

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

# Recover exact held-out kriging information without rerunning models.
kriging_files <- file.path(
  directory, fold_bounds$fold, "city_kriging.csv"
)
if (!all(file.exists(kriging_files))) {
  stop("Per-fold city kriging diagnostics are unavailable: ",
       paste(kriging_files[!file.exists(kriging_files)],
             collapse = ", "))
}
kriging <- do.call(rbind, lapply(seq_along(kriging_files), function(i) {
  x <- readr::read_csv(kriging_files[i], show_col_types = FALSE)
  if (!all(c("year", "city", "yhat", "mspe") %in% names(x))) {
    stop("Incomplete city kriging diagnostics for ",
         fold_bounds$fold[i])
  }
  data.frame(
    fold = fold_bounds$fold[i],
    city = as.character(x$city),
    year = as.integer(x$year),
    kriging_index_mean = as.numeric(x$yhat),
    kriging_mspe = as.numeric(x$mspe)
  )
}))
key_string <- function(x) paste(x$fold, x$city, x$year, sep = "|")
if (anyDuplicated(key_string(kriging))) {
  stop("Duplicate fold/city/year keys in kriging diagnostics.")
}
kriging_idx <- match(key_string(components), key_string(kriging))
if (anyNA(kriging_idx)) {
  stop("Unable to match all held-out scores to kriging diagnostics.")
}
components$kriging_index_mean <- kriging$kriging_index_mean[kriging_idx]
components$kriging_mspe <- kriging$kriging_mspe[kriging_idx]
fold_idx <- match(components$fold, covariance$fold)
sigmaY2 <- as.numeric(covariance$process_variance[fold_idx])
if (anyNA(sigmaY2) || any(!is.finite(sigmaY2)) ||
    any(sigmaY2 <= 0) ||
    any(!is.finite(components$kriging_mspe))) {
  stop("Invalid fold-specific process variances or kriging MSPE.")
}
components$information_fraction <- pmin(
  pmax(1 - components$kriging_mspe / sigmaY2, 0), 1
)
components$proxy_minus_projection_intercept <-
  components$Xstar - components$alpha

# The highest baseline log-score points can now be inspected for
# weak kriging information, near-zero beta, and tiny residual variance.
display_cols <- c(
  "fold", "city", "year", "Xstar", "kriging_index_mean",
  "kriging_mspe", "information_fraction", "alpha", "beta",
  "vdelta", "mu_prior", "V_prior", "predictive_mean",
  "predictive_residual", "proxy_minus_projection_intercept",
  "variance_without_u", "standardized_residual_without_u",
  "log_score_without_u", "log_score_with_u"
)
extreme <- components[
  order(components$log_score_without_u, decreasing = TRUE),
  display_cols
]
extreme <- utils::head(extreme, 30L)
readr::write_csv(
  extreme, file.path(directory, "cv_largest_zero_variance_scores.csv")
)

low_information <- components[
  components$information_fraction < 0.01,
  display_cols,
  drop = FALSE
]
readr::write_csv(
  low_information,
  file.path(directory, "cv_low_information_city_years.csv")
)

# Quantify how much one chronological block influences the selected u.
fold_profile_file <- file.path(directory, "cv_score_by_fold.csv")
if (file.exists(fold_profile_file)) {
  fold_profile <- readr::read_csv(
    fold_profile_file, show_col_types = FALSE
  )
  required_profile <- c(
    "fold", "u", "n", "predictive_log_score"
  )
  if (!all(required_profile %in% names(fold_profile)) ||
      anyDuplicated(fold_profile[c("fold", "u")])) {
    stop("The per-fold CV profile is missing or inconsistent.")
  }
  leave_one_out <- do.call(rbind, lapply(
    c("none", fold_bounds$fold),
    function(excluded) {
      z <- if (excluded == "none") fold_profile else
        fold_profile[fold_profile$fold != excluded, , drop = FALSE]
      weights <- tapply(z$n, z$fold, function(v) unique(v))
      if (any(lengths(weights) != 1L)) {
        stop("Inconsistent held-out counts in CV fold profile.")
      }
      total_n <- sum(unlist(weights, use.names = FALSE))
      weighted <- aggregate(
        z$predictive_log_score * z$n,
        by = list(u = z$u),
        FUN = sum
      )
      weighted$score <- weighted$x / total_n
      best <- which.min(weighted$score)
      data.frame(
        excluded_fold = excluded,
        included_city_years = total_n,
        u_minimum_on_saved_grid = weighted$u[best],
        cv_score_at_grid_minimum = weighted$score[best],
        cv_score_at_zero = weighted$score[which.min(abs(weighted$u))]
      )
    }
  ))
  readr::write_csv(
    leave_one_out,
    file.path(directory, "cv_leave_one_fold_out.csv")
  )
}

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
  "cv_largest_zero_variance_scores.csv, cv_low_information_city_years.csv, ",
  "cv_leave_one_fold_out.csv, cv_score_curve_zoom.png"
)
