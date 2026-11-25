here::i_am("Code/cv_observation_variance.R")

# Five-fold temporally blocked predictive-log-score calibration of a
# nonnegative additive measurement discrepancy variance (Celsius^2).
#
# Each fold fits REACHES covariance using training years only, but
# kriges documentary observations in both training and held-out years
# with that fitted covariance. The quantile-mapping reference and joint
# Monte Carlo reference are training-year-only; no held-out proxy is
# used to estimate either. The annual LME state prior is independent
# of REACHES and remains fixed across folds.
#
# Run from the repository root, after generating the LME/Prior inputs:
# Rscript --vanilla Code/cv_observation_variance.R
#
# Results are isolated in Output/Intermediate/ObservationVariance/.
# Nothing in the standard calibration, kriging, or assimilation output
# directories is overwritten.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(sp)
library(np)
library(ggplot2)

# Import only the reusable Figure 9(d) mathematical definitions.
# Do not source Figure9d.R: that would execute the full figure pipeline.
required_figure9_helpers <- c("annual_unconditional_variance", "predictive_log_score")
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

years <- 1368:1911
analysis_years <- years
folds <- data.frame(
  fold = sprintf("fold_%d", 1:5),
  start = c(1368L, 1501L, 1601L, 1701L, 1801L),
  end = c(1500L, 1600L, 1700L, 1800L, 1911L)
)
if (!identical(unlist(lapply(seq_len(nrow(folds)), function(i) {
  seq(folds$start[i], folds$end[i])
}), use.names = FALSE), as.integer(years))) {
  stop("The five historical folds must partition the complete time range.")
}

input_mode <- Sys.getenv("REACHES_TEMPERATURE_INPUT", unset = "original")
if (!input_mode %in% c("original", "augmented")) {
  stop("REACHES_TEMPERATURE_INPUT must be original or augmented.")
}
source_path <- if (input_mode == "original") {
  here::here("Data", "temperature index value.v1.xlsx")
} else {
  here::here(
    "Output", "Intermediate", "REACHES",
    "temperature_site_year_augmented.csv"
  )
}
if (!file.exists(source_path)) {
  stop("REACHES source is missing: ", source_path)
}
source_md5 <- unname(tools::md5sum(source_path))
output_dir <- here::here(
  "Output", "Intermediate", "ObservationVariance", input_mode
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_mc <- suppressWarnings(as.integer(Sys.getenv(
  "REACHES_CV_MC_SIZE", unset = "100000"
)))
if (length(n_mc) != 1L || is.na(n_mc) || n_mc < 1000L) {
  stop("REACHES_CV_MC_SIZE must be an integer >= 1000.")
}
calibration_engine <- Sys.getenv(
  "REACHES_CV_CALIBRATION_ENGINE", unset = "analytic"
)
if (!calibration_engine %in% c("analytic", "monte_carlo")) {
  stop("Invalid REACHES_CV_CALIBRATION_ENGINE.")
}
projection_method <- Sys.getenv(
  "REACHES_CV_PROJECTION_METHOD",
  unset = "empirical_reference_joint_mc_v1"
)
if (!projection_method %in% c(
  "empirical_reference_joint_mc_v1",
  "empirical_reference_zero_centered_joint_mc_v2"
)) stop("Invalid REACHES_CV_PROJECTION_METHOD.")
projection_center <- if (
  projection_method == "empirical_reference_joint_mc_v1"
) "reference_mean" else "zero"
rho_lookup_size <- 4001L
mapping_grid_size <- 2001L
lme_cdf_grid_size <- 5001L
latent_tail_probability <- 1e-14
variance_tolerance <- 1e-10
coordinate_tolerance <- 1e-8

# This option permits a fast single-fold implementation check. A partial
# run writes diagnostics but NEVER publishes a five-fold estimate.
only_fold <- Sys.getenv("REACHES_CV_ONLY_FOLD", unset = "")
active_folds <- folds
if (nzchar(only_fold)) {
  active_folds <- folds[folds$fold == only_fold, , drop = FALSE]
  if (nrow(active_folds) != 1L) {
    stop("REACHES_CV_ONLY_FOLD must be fold_1, ..., fold_5.")
  }
}

city_config <- list(
  HongKong = list(code = "H", long = 114.167, lat = 22.333),
  Shanghai = list(code = "S", long = 121.433, lat = 31.167),
  Beijing = list(code = "B", long = 116.283, lat = 39.933)
)

parameter_names <- c(
  "mtB.csv", "muB.csv", "rtB.csv",
  "mtS.csv", "muS.csv", "rtS.csv",
  "mtH.csv", "muH.csv", "rtH.csv"
)
candidate_inputs <- list(
  generated = list(
    lme = here::here(
      "Output", "Intermediate", "LME", "WorldClim_Corrected",
      "lme_city3_annual_1368_1911.csv"
    ),
    parameters = here::here("Output", "Intermediate", "Prior")
  ),
  precomputed = list(
    lme = here::here(
      "Data", "LME data", "precomputed",
      "lme_city3_annual_1368_1911.csv"
    ),
    parameters = here::here("Data", "par")
  )
)
complete_set <- function(z) {
  all(file.exists(c(z$lme, file.path(z$parameters, parameter_names))))
}
valid_inputs <- Filter(complete_set, candidate_inputs)
if (length(valid_inputs) == 0L) {
  stop(
    "Complete LME city and penalized-prior inputs are missing. ",
    "Run the LME preparation and Prior.R, or provide the complete ",
    "precomputed city input set."
  )
}
input_files <- valid_inputs[[1L]]
message("Selected LME/prior input set: ", names(valid_inputs)[1L])

# Reuse the exact kriging and joint-MC implementations without sourcing
# their top-level pipelines (which generate large full-grid products).
load_function_definitions <- function(path, names_required) {
  expressions <- parse(file = path)
  loaded <- character()
  for (expr in expressions) {
    if (!is.call(expr) ||
        !identical(expr[[1L]], as.name("<-")) ||
        !is.symbol(expr[[2L]]) ||
        !is.call(expr[[3L]]) ||
        !identical(expr[[3L]][[1L]], as.name("function"))) {
      next
    }
    name <- as.character(expr[[2L]])
    if (!name %in% names_required) next
    assign(name, eval(expr[[3L]], envir = .GlobalEnv),
           envir = .GlobalEnv)
    loaded <- c(loaded, name)
  }
  missing <- setdiff(names_required, loaded)
  if (length(missing) > 0L) {
    stop("Expected reusable function definitions missing from ",
         path, ": ", paste(missing, collapse = ", "))
  }
}
load_function_definitions(
  here::here("Code", "Get_tempe_all_data.R"),
  c("Ez_discrete", "Varz_discrete", "EZstar_h",
    "ordinal_covariance", "build_covariance_lookup",
    "solve_symmetric_system", "predict_grid_one_year")
)
load_function_definitions(
  here::here("Code", "Figure9d.R"),
  c("normalize_city_name", "read_lme_city",
    "read_penalized_parameter", "read_city_prior",
    "build_quantile_mapping", "compute_measurement_parameters_mc")
)

lme_city_data <- readr::read_csv(input_files$lme, show_col_types = FALSE)
required_lme <- c("city", "long", "lati", "member", "year",
                  "temperature_kelvin")
if (!all(required_lme %in% names(lme_city_data))) {
  stop("The LME city input is missing required columns.")
}
lme_city_data <- lme_city_data %>%
  transmute(
    city = normalize_city_name(city),
    long = as.numeric(long),
    lat = as.numeric(lati),
    member = as.character(member),
    year = as.integer(year),
    temperature_celsius = as.numeric(temperature_kelvin) - 273.15
  ) %>%
  filter(year %in% years) %>%
  arrange(city, member, year)
if (anyNA(lme_city_data)) stop("Invalid LME city observations.")

# Cache the independent LME sample and annual prior, not the
# documentary calibration, across all REACHES held-out blocks.
city_fixed <- lapply(names(city_config), function(name) {
  config <- city_config[[name]]
  lme <- read_lme_city(name, config, lme_city_data)
  prior <- read_city_prior(config$code, input_files$parameters)
  list(lme = lme, prior = prior,
       state_variance = annual_unconditional_variance(prior$M, prior$r2))
})
names(city_fixed) <- names(city_config)

rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Rscript executable was not found.")

fit_covariance_fold <- function(fold) {
  fold_dir <- file.path(output_dir, fold$fold)
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  model_file <- file.path(fold_dir, "calibration_parameters.rds")
  compatible <- FALSE
  if (file.exists(model_file)) {
    old <- try(readRDS(model_file), silent = TRUE)
    compatible <- !inherits(old, "try-error") &&
      identical(old$temperature_input_mode, input_mode) &&
      identical(old$temperature_source_md5, source_md5) &&
      identical(old$calibration_script_md5, unname(tools::md5sum(
        here::here("Code", "prepare_calibration.R")
      ))) &&
      identical(as.integer(old$cv_holdout), c(fold$start, fold$end)) &&
      identical(old$vario_fit2$calibration_engine, calibration_engine) &&
      !is.null(old$all_site_years)
  }
  if (!compatible) {
    log_file <- file.path(fold_dir, "calibration.log")
    message("Fitting training-only spatial covariance: ", fold$fold)
    env <- c(
      paste0("REACHES_TEMPERATURE_INPUT=", input_mode),
      paste0("REACHES_CALIBRATION_ENGINE=", calibration_engine),
      paste0("REACHES_CV_HOLDOUT_START=", fold$start),
      paste0("REACHES_CV_HOLDOUT_END=", fold$end),
      paste0("REACHES_CV_CALIBRATION_OUTPUT_DIR=", fold_dir)
    )
    status <- system2(
      rscript,
      args = c("--vanilla",
               shQuote(here::here("Code", "prepare_calibration.R"))),
      stdout = log_file, stderr = log_file, env = env
    )
    if (!identical(as.integer(status), 0L)) {
      stop("Fold covariance fit failed: ", fold$fold,
           ". Inspect ", log_file, " (exit code ", status, ").")
    }
  } else {
    message("Reusing validated covariance fit: ", fold$fold)
  }
  fit <- readRDS(model_file)
  if (is.null(fit$all_site_years) ||
      nrow(fit$all_site_years) == 0L ||
      any(between(fit$temp2$year, fold$start, fold$end))) {
    stop("Held-out year leakage detected in ", fold$fold, ".")
  }
  fit
}

cuts <- c(-Inf, -1.5, -0.5, 0.5, Inf)
vals <- c(-2, -1, 0, 1)
fold_observations <- vector("list", nrow(active_folds))
fold_covariance <- vector("list", nrow(active_folds))

for (k in seq_len(nrow(active_folds))) {
  fold <- active_folds[k, , drop = FALSE]
  fit <- fit_covariance_fold(fold)
  variance_process <- as.numeric(fit$vario_fit2$psill2)
  variance_nugget <- as.numeric(fit$vario_fit2$psill1)
  spatial_range <- as.numeric(fit$vario_fit2$range)
  if (any(!is.finite(c(variance_process, variance_nugget,
                        spatial_range))) ||
      variance_process <= 0 || variance_nugget < 0 ||
      spatial_range <= 0) {
    stop("Invalid training-only covariance parameters in ", fold$fold)
  }

  all_spatial <- as.data.frame(fit$all_site_years[
    , c("level", "year", "long", "lat")
  ])
  sp::coordinates(all_spatial) <- ~ long + lat
  sp::proj4string(all_spatial) <- sp::CRS(
    "+proj=longlat +datum=WGS84"
  )
  all_years <- sort(unique(as.integer(all_spatial$year)))
  target_coords <- do.call(rbind, lapply(city_config, function(z) {
    c(z$long, z$lat)
  }))
  colnames(target_coords) <- c("long", "lat")
  sigma_total <- sqrt(variance_process + variance_nugget)
  mean_z <- Ez_discrete(sigma_total, cuts, vals)
  var_z <- Varz_discrete(sigma_total, cuts, vals, mean_z)
  ezstar_h <- EZstar_h(sigma_total, cuts, vals)
  # Reuse earlier city kriging only by explicit opt-in and only after
  # checking source identity, fold covariance, complete year/city keys,
  # and the saved baseline calibration parameters. Changing the
  # measurement projection does not affect the kriging predictions.
  city_file <- file.path(
    output_dir, as.character(fold$fold), "city_kriging.csv"
  )
  old_selection_file <- file.path(output_dir, "selected_u.csv")
  old_covariance_file <- file.path(output_dir, "fold_covariance.csv")
  reuse_city <- identical(
    Sys.getenv("REACHES_CV_REUSE_CITY_KRIGING", unset = "0"), "1"
  ) && all(file.exists(c(
    city_file, old_selection_file, old_covariance_file
  )))
  if (reuse_city) {
    old_selection <- readr::read_csv(
      old_selection_file, show_col_types = FALSE
    )
    old_covariance <- readr::read_csv(
      old_covariance_file, show_col_types = FALSE
    ) %>%
      filter(.data$fold == .env$fold$fold)
    existing_city <- readr::read_csv(
      city_file, show_col_types = FALSE
    )
    required_city <- c("year", "city", "yhat", "mspe", "nu")
    reusable <- nrow(old_selection) == 1L &&
      all(c("temperature_input_mode", "temperature_source_md5",
            "calibration_engine", "calibration_script_md5") %in%
            names(old_selection)) &&
      identical(as.character(old_selection$temperature_input_mode[1L]),
                input_mode) &&
      identical(as.character(old_selection$temperature_source_md5[1L]),
                source_md5) &&
      identical(as.character(old_selection$calibration_engine[1L]),
                calibration_engine) &&
      identical(as.character(old_selection$calibration_script_md5[1L]),
                unname(tools::md5sum(
                  here::here("Code", "prepare_calibration.R")
                ))) &&
      nrow(old_covariance) == 1L &&
      all(c("process_variance", "nugget_variance", "range_km",
            "n_train", "n_all") %in% names(old_covariance)) &&
      isTRUE(all.equal(old_covariance$process_variance[1L],
                       variance_process, tolerance = 1e-10)) &&
      isTRUE(all.equal(old_covariance$nugget_variance[1L],
                       variance_nugget, tolerance = 1e-10)) &&
      isTRUE(all.equal(old_covariance$range_km[1L],
                       spatial_range, tolerance = 1e-10)) &&
      as.integer(old_covariance$n_train[1L]) == nrow(fit$temp2) &&
      as.integer(old_covariance$n_all[1L]) ==
        nrow(fit$all_site_years) &&
      all(required_city %in% names(existing_city))
    if (reusable) {
      reusable <- nrow(existing_city) ==
        length(all_years) * length(city_config) &&
        !anyDuplicated(existing_city[c("year", "city")]) &&
        setequal(as.integer(existing_city$year), all_years) &&
        setequal(as.character(existing_city$city),
                 names(city_config)) &&
        all(is.finite(as.matrix(existing_city[
          c("yhat", "mspe", "nu")
        ]))) &&
        all(existing_city$mspe >= 0) &&
        all(existing_city$mspe <= variance_process + 1e-7) &&
        all(existing_city$nu >= 0) &&
        isTRUE(all.equal(
          existing_city$nu^2, existing_city$mspe,
          tolerance = 1e-6
        ))
    }
    if (!reusable) {
      stop("Existing city kriging cache failed validation for ",
           fold$fold, ". Rerun without REACHES_CV_REUSE_CITY_KRIGING=1.")
    }
    city_predictions <- existing_city
    message("Reusing validated city kriging for ", fold$fold)
  } else {
    message("Preparing shared ordinal covariance lookup: ", fold$fold)
    lookup <- build_covariance_lookup(
      sigma_Y2 = variance_process, sigma_E2 = variance_nugget,
      sigma = sigma_total, number_of_points = rho_lookup_size
    )

    predictions <- vector("list", length(all_years))
    for (j in seq_along(all_years)) {
      pred <- predict_grid_one_year(
        current_year = all_years[j],
        observation_data = all_spatial,
        prediction_coordinates = target_coords,
        sigma_Y2 = variance_process,
        sigma_E2 = variance_nugget,
        alpha_km = spatial_range,
        total_variance = variance_process + variance_nugget,
        mean_z = mean_z,
        variance_z = var_z,
        ezstar_h = ezstar_h,
        covariance_lookup_function = lookup$function_object
      )
      predictions[[j]] <- data.frame(
        year = all_years[j],
        city = names(city_config),
        yhat = as.numeric(pred$mean),
        mspe = as.numeric(pred$variance)
      )
      if (j %% 100L == 0L) {
        message("  ", fold$fold, " city kriging: ", j, "/",
                length(all_years), " years")
      }
    }
    city_predictions <- bind_rows(predictions) %>%
      mutate(nu = sqrt(pmax(mspe, 0)))
    readr::write_csv(city_predictions, city_file)
  }

  city_results <- vector("list", length(city_config))
  for (ci in seq_along(city_config)) {
    name <- names(city_config)[ci]
    city <- city_predictions %>% filter(.data$city == name)
    training <- city %>%
      filter(!between(year, fold$start, fold$end))
    held_out <- city %>%
      filter(between(year, fold$start, fold$end))
    if (nrow(training) < 10L || nrow(held_out) < 2L) {
      stop("Insufficient training/held-out documentary years in ",
           fold$fold, ", city ", name)
    }
    independent <- city_fixed[[name]]
    mapping <- build_quantile_mapping(
      lme_sample = independent$lme$sample,
      yhat = training$yhat,
      nu = training$nu,
      sigmaY2 = variance_process,
      bandwidth_seed = 10L + ci,
      mapping_grid_size = mapping_grid_size,
      lme_cdf_grid_size = lme_cdf_grid_size,
      latent_tail_probability = latent_tail_probability
    )
    implied_vhat <- pmin(
      pmax(variance_process - held_out$mspe, 0),
      variance_process
    )
    measurement <- compute_measurement_parameters_mc(
      sigmaY2 = variance_process,
      vhat = implied_vhat,
      reference_yhat = training$yhat,
      reference_nu = training$nu,
      target_nu = held_out$nu,
      g = mapping$g,
      center_mode = projection_center,
      n_mc = n_mc,
      seed = 1L + ci
    )
    prior_year_idx <- match(held_out$year, years)
    if (anyNA(prior_year_idx)) stop("Invalid held-out year.")
    mu_prior <- independent$prior$mu[prior_year_idx]
    V_prior <- independent$state_variance[prior_year_idx]
    Xstar <- as.numeric(mapping$g(held_out$yhat))
    m <- measurement$alpha + measurement$beta * mu_prior
    base <- measurement$beta^2 * V_prior + measurement$vdelta
    city_results[[ci]] <- data.frame(
      fold = fold$fold,
      year = held_out$year,
      city = name,
      Xstar = Xstar,
      mu_prior = mu_prior,
      V_prior = V_prior,
      alpha = measurement$alpha,
      beta = measurement$beta,
      vdelta = measurement$vdelta,
      predictive_mean = m,
      baseline_predictive_variance = base,
      predictive_residual = Xstar - m
    )
    if (any(!is.finite(unlist(city_results[[ci]][, -c(1, 3)]))) ||
        any(base <= 0)) {
      stop("Non-finite held-out score data for ", fold$fold, " ", name)
    }
  }
  fold_observations[[k]] <- bind_rows(city_results)
  fold_covariance[[k]] <- data.frame(
    fold = fold$fold, holdout_start = fold$start,
    holdout_end = fold$end, n_train = nrow(fit$temp2),
    n_all = nrow(fit$all_site_years),
    process_variance = variance_process,
    nugget_variance = variance_nugget,
    range_km = spatial_range
  )
  readr::write_csv(
    fold_observations[[k]],
    file.path(output_dir, fold$fold, "heldout_score_components.csv")
  )
  invisible(gc())
}

all_scores <- bind_rows(fold_observations)
covariance_table <- bind_rows(fold_covariance)
# A single-fold smoke test must never overwrite the pooled five-fold
# products, including any previously fitted selected_u.csv.
is_partial_run <- nrow(active_folds) != nrow(folds)
results_dir <- if (is_partial_run) {
  file.path(output_dir, "smoke", as.character(active_folds$fold[1L]))
} else {
  output_dir
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  all_scores, file.path(results_dir, "heldout_score_components.csv")
)
readr::write_csv(
  covariance_table, file.path(results_dir, "fold_covariance.csv")
)

if (is_partial_run) {
  message("One-fold diagnostic run complete: ", results_dir,
          ". Existing five-fold results were not modified.")
} else {
  if (anyDuplicated(all_scores[c("year", "city")])) {
    stop("Duplicate held-out city/year evaluations.")
  }
  if (!identical(
    sort(unique(as.character(all_scores$fold))),
    sort(folds$fold)
  )) stop("Some historical fold scores are missing.")

  residual <- all_scores$predictive_residual
  base <- all_scores$baseline_predictive_variance
  loss <- function(u) predictive_log_score(residual, base, u)
  derivative <- function(u) {
    S <- base + u
    mean(0.5 * (1 / S - residual^2 / S^2))
  }
  upper <- max(1, 10 * stats::var(residual), 10 * max(base))
  for (i in seq_len(10L)) {
    if (derivative(upper) > 0) break
    upper <- 4 * upper
  }
  if (!is.finite(upper) || !is.finite(derivative(upper)) ||
      derivative(upper) <= 0) {
    stop("Could not bracket the additional-variance objective.")
  }

  # Heterogeneous base variances need not give a unimodal profile.
  # Scan logarithmically, refine every detected local minimum, and
  # compare with the constrained boundary u = 0.
  coarse <- sort(unique(c(
    0, exp(seq(log(1e-10), log(upper), length.out = 401L))
  )))
  coarse_loss <- vapply(coarse, loss, numeric(1L))
  n_grid <- length(coarse)
  local <- which(
    coarse_loss[2:(n_grid - 1L)] <= coarse_loss[1:(n_grid - 2L)] &
    coarse_loss[2:(n_grid - 1L)] <= coarse_loss[3:n_grid]
  ) + 1L
  candidates <- c(0, upper)
  for (idx in local) {
    refined <- stats::optimize(
      loss, c(coarse[idx - 1L], coarse[idx + 1L]), tol = 1e-10
    )
    candidates <- c(candidates, coarse[idx], refined$minimum)
  }
  u_opt <- candidates[which.min(
    vapply(candidates, loss, numeric(1L))
  )]
  if (u_opt < 1e-8) u_opt <- 0
  u_grid <- sort(unique(c(coarse, u_opt)))
  cv_grid <- data.frame(
    u = u_grid,
    predictive_log_score = vapply(u_grid, loss, numeric(1L))
  )
  score_by_fold <- bind_rows(lapply(seq_len(nrow(folds)), function(i) {
    z <- all_scores %>% filter(.data$fold == folds$fold[i])
    data.frame(
      fold = folds$fold[i], u = u_grid,
      n = nrow(z),
      predictive_log_score = vapply(
        u_grid, function(u) predictive_log_score(
          z$predictive_residual, z$baseline_predictive_variance, u
        ), numeric(1L)
      )
    )
  }))
  score_by_city <- bind_rows(lapply(names(city_config), function(name) {
    z <- all_scores %>% filter(.data$city == name)
    data.frame(
      city = name, u = u_grid, n = nrow(z),
      predictive_log_score = vapply(
        u_grid, function(u) predictive_log_score(
          z$predictive_residual, z$baseline_predictive_variance, u
        ), numeric(1L)
      )
    )
  }))
  selected <- data.frame(
    selected_u = u_opt,
    variance_units = "celsius_squared",
    temperature_input_mode = input_mode,
    temperature_source_md5 = source_md5,
    number_of_folds = nrow(folds),
    number_of_heldout_city_years = nrow(all_scores),
    cv_score_at_zero = loss(0),
    cv_score_selected = loss(u_opt),
    source_file = normalizePath(source_path),
    calibration_engine = calibration_engine,
    calibration_script_md5 = unname(tools::md5sum(
      here::here("Code", "prepare_calibration.R")
    )),
    mc_size = n_mc,
    measurement_projection_method = projection_method,
    score_type = "pooled_city_year_predictive_log_score",
    lme_prior_source = names(valid_inputs)[1L],
    cv_script_md5 = unname(tools::md5sum(
      here::here("Code", "cv_observation_variance.R")
    )),
    shared_manifest_md5 = Sys.getenv(
      "REACHES_CV_SHARED_MANIFEST_MD5", unset = ""
    )
  )
  readr::write_csv(cv_grid, file.path(output_dir, "cv_score_curve.csv"))
  readr::write_csv(score_by_fold,
                   file.path(output_dir, "cv_score_by_fold.csv"))
  readr::write_csv(score_by_city,
                   file.path(output_dir, "cv_score_by_city.csv"))
  readr::write_csv(selected, file.path(output_dir, "selected_u.csv"))

  plot <- ggplot(cv_grid, aes(x = u, y = predictive_log_score)) +
    geom_line(linewidth = 0.7) +
    geom_vline(xintercept = u_opt, linetype = "dashed") +
    labs(x = expression("Additional variance, " * u * " (" *
                        degree * "C)"^2),
         y = "Held-out predictive log score") +
    theme_minimal(base_size = 12)
  ggsave(
    file.path(output_dir, "cv_score_curve.png"), plot,
    width = 7, height = 4.5, dpi = 300
  )
  message(
    "Selected extra measurement variance u = ",
    signif(u_opt, 7), " Celsius^2; log score: ",
    signif(loss(0), 7), " (u=0) -> ",
    signif(loss(u_opt), 7), " (fitted)."
  )
  message("Fit saved to: ", file.path(output_dir, "selected_u.csv"))
}
