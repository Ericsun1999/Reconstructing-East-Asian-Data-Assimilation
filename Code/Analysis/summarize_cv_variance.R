here::i_am("Code/Analysis/summarize_cv_variance.R")

# Generate a reviewer-facing five-fold summary from EXISTING CV outputs.
# No REACHES calibration, kriging, or prior fitting is performed here.
# Run from the Original reconstruction workspace after the full CV.

library(here)
library(readr)
library(dplyr)

output_dir <- here::here(
  "Output", "Intermediate", "ObservationVariance", "original"
)
required_files <- file.path(
  output_dir,
  c("fold_covariance.csv", "heldout_score_components.csv", "selected_u.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Existing full CV files are missing: ",
       paste(missing_files, collapse = "; "))
}

covariance <- readr::read_csv(required_files[1L], show_col_types = FALSE)
scores <- readr::read_csv(required_files[2L], show_col_types = FALSE)
selection <- readr::read_csv(required_files[3L], show_col_types = FALSE)

required_covariance <- c(
  "fold", "holdout_start", "holdout_end", "n_train", "n_all",
  "process_variance", "nugget_variance", "range_km"
)
required_scores <- c(
  "fold", "year", "city", "predictive_residual",
  "baseline_predictive_variance"
)
if (!all(required_covariance %in% names(covariance)) ||
    !all(required_scores %in% names(scores)) ||
    nrow(selection) != 1L ||
    !"selected_u" %in% names(selection)) {
  stop("The existing CV outputs are incomplete or incompatible.")
}
if (nrow(covariance) != 5L ||
    !setequal(covariance$fold, paste0("fold_", 1:5)) ||
    anyDuplicated(covariance$fold) ||
    !setequal(unique(scores$fold), covariance$fold) ||
    anyDuplicated(scores[c("year", "city")]) ||
    any(!is.finite(scores$predictive_residual)) ||
    any(!is.finite(scores$baseline_predictive_variance)) ||
    any(scores$baseline_predictive_variance <= 0)) {
  stop("The five-fold CV output contains missing, repeated or invalid scores.")
}
u <- as.numeric(selection$selected_u[1L])
if (length(u) != 1L || !is.finite(u) || u < 0) {
  stop("The selected additional observation variance is invalid.")
}
if ("temperature_input_mode" %in% names(selection) &&
    as.character(selection$temperature_input_mode[1L]) != "original") {
  stop("The shared observation variance must be fitted on Original input.")
}

predictive_log_score <- function(residual, base_variance, u) {
  S <- base_variance + u
  mean(0.5 * (log(2 * pi * S) + residual^2 / S))
}

document_years <- sort(unique(as.integer(scores$year)))
if (any(!document_years %in% 1368:1911)) {
  stop("Held-out documentary years fall outside 1368–1911.")
}

summary_rows <- lapply(seq_len(nrow(covariance)), function(i) {
  row <- covariance[i, , drop = FALSE]
  fold <- as.character(row$fold)
  held <- scores %>% filter(.data$fold == fold)
  a <- as.integer(row$holdout_start)
  b <- as.integer(row$holdout_end)
  if (nrow(held) < 2L || a > b ||
      any(held$year < a | held$year > b) ||
      any(row$n_train < 0L) || any(row$n_all <= row$n_train)) {
    stop("Invalid training or validation partition for ", fold)
  }
  train_ranges <- character()
  if (a > 1368L) train_ranges <- c(train_ranges, paste0("1368-", a - 1L))
  if (b < 1911L) train_ranges <- c(train_ranges, paste0(b + 1L, "-1911"))

  data.frame(
    fold = fold,
    training_ranges = paste(train_ranges, collapse = "; "),
    validation_start = a,
    validation_end = b,
    n_training_observations = as.integer(row$n_train),
    n_validation_observations = as.integer(row$n_all - row$n_train),
    n_training_document_years = sum(!between(document_years, a, b)),
    n_validation_document_years = n_distinct(held$year),
    n_scored_city_years = nrow(held),
    process_variance = as.numeric(row$process_variance),
    nugget_variance = as.numeric(row$nugget_variance),
    range_km = as.numeric(row$range_km),
    pooled_selected_u_celsius2 = u,
    predictive_log_score_u0 = predictive_log_score(
      held$predictive_residual, held$baseline_predictive_variance, 0
    ),
    predictive_log_score_selected_u = predictive_log_score(
      held$predictive_residual, held$baseline_predictive_variance, u
    )
  )
})
summary_table <- dplyr::bind_rows(summary_rows) %>%
  arrange(validation_start)
readr::write_csv(
  summary_table, file.path(output_dir, "cv_fold_summary.csv")
)
print(summary_table)
message("CV audit saved without refitting: ",
        file.path(output_dir, "cv_fold_summary.csv"))
