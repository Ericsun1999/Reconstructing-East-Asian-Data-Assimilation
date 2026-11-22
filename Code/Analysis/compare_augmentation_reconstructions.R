here::i_am("Code/Analysis/compare_augmentation_reconstructions.R")

# Compare archival Original and Augmented reconstructions using an
# identical observation model and an identical fixed extra variance.
# Each input must be snapshotted with
# Code/Analysis/snapshot_augmentation_reconstruction.R first.
#
# Rscript --vanilla Code/Analysis/compare_augmentation_reconstructions.R city
# Rscript --vanilla Code/Analysis/compare_augmentation_reconstructions.R grid

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[1L] %in% c("city", "grid")) {
  stop("Usage: compare_augmentation_reconstructions.R <city|grid>")
}
stage <- args[1L]
version_root <- trimws(Sys.getenv(
  "REACHES_L4_VERSION_ROOT", unset = ""
))
if (nzchar(version_root)) {
  version_root <- normalizePath(
    path.expand(version_root), winslash = "/", mustWork = TRUE
  )
  snapshots <- c(
    original = file.path(
      version_root, "original", "Output", "Intermediate",
      "AugmentationSensitivity", "original", stage
    ),
    augmented = file.path(
      version_root, "inferred-zero", "Output", "Intermediate",
      "AugmentationSensitivity", "augmented", stage
    )
  )
  result_dir <- file.path(version_root, "Comparison", stage)
} else {
  root <- here::here(
    "Output", "Intermediate", "AugmentationSensitivity"
  )
  snapshots <- setNames(
    file.path(root, c("original", "augmented"), stage),
    c("original", "augmented")
  )
  result_dir <- file.path(root, "Comparison", stage)
}
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

get_file <- function(mode, filename) {
  p <- file.path(snapshots[[mode]], filename)
  if (!file.exists(p)) stop("Missing snapshot file: ", p)
  p
}
get_data <- function(mode, filename) {
  readr::read_csv(get_file(mode, filename), show_col_types = FALSE)
}
manifests <- lapply(names(snapshots), function(mode) {
  m <- get_data(mode, "snapshot_manifest.csv")
  if (nrow(m) != 1L || !identical(as.character(m$input_mode), mode) ||
      !identical(as.character(m$stage), stage)) {
    stop("Snapshot mode/stage mismatch: ", mode)
  }
  m
})
names(manifests) <- names(snapshots)
original <- manifests$original
augmented <- manifests$augmented
if (!identical(
      as.character(original$projection_method[1L]),
      as.character(augmented$projection_method[1L])
    ) ||
    as.character(original$projection_method[1L]) !=
      "empirical_reference_joint_mc_v1") {
  stop("Both reconstructions must use the same original ",
       "measurement projection (v1) for this comparison.")
}
if (!isTRUE(all.equal(
      as.numeric(original$extra_observation_variance[1L]),
      as.numeric(augmented$extra_observation_variance[1L]),
      tolerance = 1e-10
    )) ||
    as.integer(original$measurement_mc_size[1L]) !=
      as.integer(augmented$measurement_mc_size[1L])) {
  stop("Original and Augmented must use the same fixed u and ",
       "Monte Carlo draw count. Do not confound dataset changes ",
       "with changes in the measurement model.")
}
if (identical(as.character(original$source_md5[1L]),
              as.character(augmented$source_md5[1L]))) {
  stop("The documentary source inputs are identical; expected ",
       "Original and Augmented versions.")
}

block <- function(year) {
  cut(year,
      breaks = c(1367, 1500, 1600, 1700, 1800, 1911),
      labels = c("1368-1500", "1501-1600", "1601-1700",
                 "1701-1800", "1801-1911"))
}
city_names <- c("Beijing", "Shanghai", "HongKong")
city_result <- vector("list", length(city_names))
for (i in seq_along(city_names)) {
  name <- city_names[i]
  filename <- paste0("assimilation_", name, ".csv")
  a <- get_data("original", filename)
  b <- get_data("augmented", filename)
  required <- c(
    "year", "smoothed_temperature_celsius", "smoothed_variance",
    "prior_mean", "lme_ensemble_mean_celsius",
    "reaches_index_mean", "reaches_mspe"
  )
  if (!all(required %in% names(a)) ||
      !all(required %in% names(b)) ||
      nrow(a) != 544L || nrow(b) != 544L ||
      anyDuplicated(a$year) || anyDuplicated(b$year) ||
      !setequal(a$year, 1368:1911) ||
      !setequal(b$year, 1368:1911)) {
    stop("Unaligned city model outputs for ", name)
  }
  x <- inner_join(
    a %>% select(all_of(required)),
    b %>% select(all_of(required)),
    by = "year", suffix = c("_original", "_augmented")
  ) %>% arrange(year)
  for (column in c("prior_mean", "lme_ensemble_mean_celsius")) {
    one <- x[[paste0(column, "_original")]]
    two <- x[[paste0(column, "_augmented")]]
    if (anyNA(one) || anyNA(two) ||
        max(abs(one - two)) > 1e-6) {
      stop("LME/prior inputs differ across data versions for ", name,
           ": ", column)
    }
  }
  va <- x$smoothed_variance_original
  vb <- x$smoothed_variance_augmented
  if (any(!is.finite(va)) || any(!is.finite(vb)) ||
      any(va < 0) || any(vb < 0)) {
    stop("Invalid posterior variances for ", name)
  }
  x$city <- name
  x$block <- block(x$year)
  x$delta_mean <- x$smoothed_temperature_celsius_augmented -
    x$smoothed_temperature_celsius_original
  x$sd_original <- sqrt(va)
  x$sd_augmented <- sqrt(vb)
  x$delta_sd <- x$sd_augmented - x$sd_original
  x$sd_ratio <- x$sd_augmented / pmax(x$sd_original, 1e-12)
  x$shift_over_original_sd <- abs(x$delta_mean) /
    pmax(x$sd_original, 1e-12)
  x$delta_reaches_index <- x$reaches_index_mean_augmented -
    x$reaches_index_mean_original
  x$delta_reaches_mspe <- x$reaches_mspe_augmented -
    x$reaches_mspe_original
  city_result[[i]] <- x
}
city_years <- bind_rows(city_result)
readr::write_csv(
  city_years, file.path(result_dir, "city_year_differences.csv")
)

summarize_differences <- function(z) {
  tibble::tibble(
    n = nrow(z),
    mean_signed_shift_celsius = mean(z$delta_mean),
    mean_absolute_shift_celsius = mean(abs(z$delta_mean)),
    rmse_shift_celsius = sqrt(mean(z$delta_mean^2)),
    maximum_absolute_shift_celsius = max(abs(z$delta_mean)),
    mean_original_posterior_sd = mean(z$sd_original),
    mean_augmented_posterior_sd = mean(z$sd_augmented),
    mean_sd_change_celsius = mean(z$delta_sd),
    median_absolute_shift_over_original_sd =
      median(z$shift_over_original_sd),
    fraction_shift_larger_than_original_sd =
      mean(z$shift_over_original_sd > 1),
    fraction_posterior_sd_increased = mean(z$delta_sd > 0)
  )
}
city_summary <- city_years %>%
  group_by(city) %>%
  group_modify(~ summarize_differences(.x)) %>%
  ungroup()
city_block_summary <- city_years %>%
  group_by(city, block) %>%
  group_modify(~ summarize_differences(.x)) %>%
  ungroup()
readr::write_csv(
  city_summary, file.path(result_dir, "city_summary.csv")
)
readr::write_csv(
  city_block_summary, file.path(result_dir, "city_by_block.csv")
)

# Both variations share the same LME prior; the band is the
# baseline posterior SD. It is a scale reference, not a confidence
# interval for the difference of two highly correlated posteriors.
p_shift <- ggplot(city_years, aes(x = year)) +
  geom_ribbon(
    aes(ymin = -sd_original, ymax = sd_original),
    fill = "grey80", alpha = 0.6
  ) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_line(aes(y = delta_mean), colour = "#377EB8",
            linewidth = 0.4) +
  facet_wrap(~city, ncol = 1) +
  labs(
    title = "Temperature reconstruction sensitivity to inferred zeros",
    subtitle = "Augmented minus Original; shaded: Original posterior SD",
    x = "Year", y = "Posterior mean difference (°C)"
  ) +
  theme_gray(base_size = 12)
ggsave(file.path(result_dir, "city_posterior_mean_shifts.png"),
       p_shift, width = 11, height = 8, dpi = 300, bg = "white")

p_sd <- ggplot(city_years, aes(x = year, y = sd_ratio)) +
  geom_hline(yintercept = 1, colour = "grey30", linewidth = 0.3) +
  geom_line(colour = "#377EB8", linewidth = 0.4) +
  facet_wrap(~city, ncol = 1) +
  labs(
    title = "Sensitivity of posterior uncertainty",
    subtitle = "Ratio greater than 1: larger uncertainty with inferred zeros",
    x = "Year", y = "Augmented / Original posterior SD"
  ) +
  theme_gray(base_size = 12)
ggsave(file.path(result_dir, "city_posterior_sd_ratios.png"),
       p_sd, width = 11, height = 8, dpi = 300, bg = "white")

covariance <- bind_rows(lapply(names(snapshots), function(mode) {
  fit <- readRDS(get_file(mode, "calibration_parameters.rds"))
  z <- fit$vario_fit2
  if (!all(c("psill1", "psill2", "range", "total_variance") %in%
           names(z))) stop("Missing calibrated covariance parameters.")
  tibble::tibble(
    mode = mode,
    nugget_variance = as.numeric(z$psill1),
    process_variance = as.numeric(z$psill2),
    spatial_range_km = as.numeric(z$range),
    total_variance = as.numeric(z$total_variance)
  )
}))
readr::write_csv(covariance, file.path(result_dir, "covariance_comparison.csv"))

if (stage == "grid") {
  yearly_columns <- paste0("x", 1368:1911)
  keys <- c("location_id", "lati", "long", "figure6_cluster")
  load_grid <- function(mode, kind) {
    d <- get_data(
      mode,
      paste0("assimilated_posterior_lme_grid_", kind, ".csv")
    )
    if (!all(c(keys, yearly_columns) %in% names(d)) ||
        anyDuplicated(d$location_id)) {
      stop("The spatial grid has invalid location/year structure.")
    }
    d %>%
      pivot_longer(
        cols = all_of(yearly_columns),
        names_to = "year", values_to = kind
      ) %>%
      mutate(year = as.integer(sub("^x", "", year)),
             mode = mode)
  }
  grid_model <- bind_rows(lapply(names(snapshots), function(mode) {
    means <- load_grid(mode, "mean")
    vars <- load_grid(mode, "variance")
    inner_join(
      means, vars,
      by = c(keys, "year", "mode")
    )
  }))
  keys2 <- c(keys, "year")
  grid_original <- grid_model %>%
    filter(mode == "original") %>% select(-mode)
  grid_augmented <- grid_model %>%
    filter(mode == "augmented") %>% select(-mode)
  if (nrow(grid_original) != nrow(grid_augmented)) {
    stop("Original and Augmented have different spatial grid sizes.")
  }
  grid_pairs <- inner_join(
    grid_original, grid_augmented,
    by = keys2, suffix = c("_original", "_augmented")
  )
  if (nrow(grid_pairs) != nrow(grid_original) ||
      anyNA(grid_pairs$mean_original) ||
      anyNA(grid_pairs$mean_augmented) ||
      anyNA(grid_pairs$variance_original) ||
      anyNA(grid_pairs$variance_augmented) ||
      any(grid_pairs$variance_original < 0) ||
      any(grid_pairs$variance_augmented < 0)) {
    stop("Incomplete or invalid spatial posterior comparison.")
  }
  grid_pairs <- grid_pairs %>%
    mutate(
      block = block(year),
      delta_mean = mean_augmented - mean_original,
      sd_original = sqrt(variance_original),
      sd_augmented = sqrt(variance_augmented),
      delta_sd = sd_augmented - sd_original,
      shift_over_original_sd = abs(delta_mean) /
        pmax(sd_original, 1e-12)
    )
  readr::write_csv(
    grid_pairs, file.path(result_dir, "grid_year_differences.csv")
  )
  grid_summary <- grid_pairs %>%
    group_by(location_id, lati, long, figure6_cluster) %>%
    summarise(
      mean_signed_shift_celsius = mean(delta_mean),
      mean_absolute_shift_celsius = mean(abs(delta_mean)),
      mean_sd_change_celsius = mean(delta_sd),
      fraction_shift_larger_than_original_sd =
        mean(shift_over_original_sd > 1),
      .groups = "drop"
    )
  readr::write_csv(
    grid_summary, file.path(result_dir, "grid_location_summary.csv")
  )
  grid_block_summary <- grid_pairs %>%
    group_by(block) %>%
    group_modify(~ summarize_differences(
      tibble::tibble(
        delta_mean = .x$delta_mean,
        sd_original = .x$sd_original,
        sd_augmented = .x$sd_augmented,
        delta_sd = .x$delta_sd,
        shift_over_original_sd = .x$shift_over_original_sd
      )
    )) %>%
    ungroup()
  readr::write_csv(
    grid_block_summary, file.path(result_dir, "grid_by_block.csv")
  )
  if (requireNamespace("maps", quietly = TRUE) &&
      requireNamespace("mapproj", quietly = TRUE) &&
      requireNamespace("RColorBrewer", quietly = TRUE)) {
    limit <- max(abs(grid_summary$mean_signed_shift_celsius))
    if (!is.finite(limit) || limit <= 0) limit <- 1e-8
    p_map <- ggplot() +
      borders(
        "world", xlim = c(95, 126), ylim = c(19, 45),
        fill = NA, colour = "grey30"
      ) +
      geom_point(
        data = grid_summary,
        aes(x = long, y = lati, colour = mean_signed_shift_celsius),
        size = 2.5
      ) +
      coord_map(xlim = c(98, 124.5), ylim = c(19, 42.5)) +
      scale_colour_gradientn(
        colours = rev(RColorBrewer::brewer.pal(9, "RdBu")),
        limits = c(-limit, limit),
        name = "Mean difference (°C)"
      ) +
      labs(
        title = "Mean reconstruction difference, 1368–1911",
        subtitle = "Augmented minus Original",
        x = "long", y = "lat"
      ) +
      theme_gray(base_size = 12)
    ggsave(
      file.path(result_dir, "grid_mean_difference_map.png"),
      p_map, width = 8, height = 6, dpi = 300, bg = "white"
    )
  }
}
comparison_manifest <- data.frame(
  stage = stage,
  projection_method = as.character(original$projection_method[1L]),
  fixed_extra_variance = as.numeric(
    original$extra_observation_variance[1L]
  ),
  mc_size = as.integer(original$measurement_mc_size[1L]),
  original_source_md5 = as.character(original$source_md5[1L]),
  augmented_source_md5 = as.character(augmented$source_md5[1L]),
  compared_at = as.character(Sys.time())
)
readr::write_csv(
  comparison_manifest,
  file.path(result_dir, "comparison_manifest.csv")
)
print(city_summary)
print(covariance)
message("Verified L4 sensitivity comparison saved to: ", result_dir)
