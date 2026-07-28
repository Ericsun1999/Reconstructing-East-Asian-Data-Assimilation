here::i_am("Code/Prior.R")

library(here)
library(readr)
library(tidyr)
library(dplyr)

# ============================================================
# Configuration
# ============================================================

analysis_years <- 1368:1911

# The three smoothing parameters are selected by leave-one-
# ensemble-member-out cross-validation over a logarithmic grid.
#
# A direct closed-form minimizer is not available because every
# CV evaluation requires refitting the nonlinear penalized model.
# This coarse grid contains 7^3 = 343 candidate triples and is
# used to keep the reproducibility run computationally feasible.
#
# For a finer search, increase length.out (for example to 13 or
# 25) and/or expand the exponent range. If a selected lambda lies
# on 10^0 or 10^3, the boundary warning below indicates that the
# grid should be expanded or refined.


lambda1_grid <- 10^seq(2, 3, length.out = 3)
lambda2_grid <- 10^seq(1, 3, length.out = 5)
lambda3_grid <- 10^seq(2, 3, length.out = 3)

max_iter <- 50
tol <- 1e-8
verbose <- TRUE

input_mode <- "auto"

allowed_input_modes <- c(
  "auto",
  "generated",
  "precomputed"
)

output_dir <- here::here(
  "Output",
  "Intermediate",
  "Prior"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!input_mode %in% allowed_input_modes) {
  stop(
    "input_mode must be one of: ",
    paste(
      allowed_input_modes,
      collapse = ", "
    )
  )
}

lme_input_files <- c(
  generated = here::here(
    "Output",
    "Intermediate",
    "LME",
    "lme_city3_annual_1368_1911.csv"
  ),
  precomputed = here::here(
    "Data",
    "LME data",
    "precomputed",
    "lme_city3_annual_1368_1911.csv"
  )
)

select_lme_input_file <- function(
    input_mode,
    input_files) {

  generated_file <- unname(
    input_files["generated"]
  )

  precomputed_file <- unname(
    input_files["precomputed"]
  )

  generated_available <- file.exists(
    generated_file
  )

  precomputed_available <- file.exists(
    precomputed_file
  )

  if (input_mode == "generated") {
    if (!generated_available) {
      stop(
        "Generated LME city input was not found: ",
        generated_file,
        "\nRun Code/DataPreparation/prepare_lme_annual.R first."
      )
    }

    return(
      generated_file
    )
  }

  if (input_mode == "precomputed") {
    if (!precomputed_available) {
      stop(
        "Precomputed LME city input was not found: ",
        precomputed_file
      )
    }

    return(
      precomputed_file
    )
  }

  if (generated_available) {
    return(
      generated_file
    )
  }

  if (precomputed_available) {
    return(
      precomputed_file
    )
  }

  stop(
    "Neither generated nor precomputed LME city input was found."
  )
}

lme_input_file <- select_lme_input_file(
  input_mode = input_mode,
  input_files = lme_input_files
)

message(
  "Prior input selected: ",
  lme_input_file
)

lme_city_data <- readr::read_csv(
  lme_input_file,
  show_col_types = FALSE
)

required_lme_columns <- c(
  "city",
  "member",
  "year",
  "temperature_kelvin"
)

missing_lme_columns <- setdiff(
  required_lme_columns,
  names(
    lme_city_data
  )
)

if (length(missing_lme_columns) > 0L) {
  stop(
    "The LME city input is missing columns: ",
    paste(
      missing_lme_columns,
      collapse = ", "
    )
  )
}

normalize_city_name <- function(city) {

  city <- tolower(
    gsub(
      "[[:space:]_-]",
      "",
      as.character(
        city
      )
    )
  )

  dplyr::recode(
    city,
    hongkong = "HongKong",
    shanghai = "Shanghai",
    beijing = "Beijing",
    .default = NA_character_
  )
}

lme_city_data <- lme_city_data %>%
  transmute(
    city = normalize_city_name(
      city
    ),
    member = as.character(
      member
    ),
    year = as.integer(
      year
    ),
    temperature_celsius =
      as.numeric(
        temperature_kelvin
      ) -
      273.15
  ) %>%
  filter(
    year %in%
      analysis_years
  ) %>%
  arrange(
    city,
    member,
    year
  )

if (
  anyNA(
    lme_city_data
  ) ||
    any(
      !is.finite(
        lme_city_data$temperature_celsius
      )
    )
) {
  stop(
    "The LME city input contains invalid or missing values."
  )
}

city_config <- list(
  HongKong = list(
    code = "H"
  ),
  Shanghai = list(
    code = "S"
  ),
  Beijing = list(
    code = "B"
  )
)

# ============================================================
# Model functions
# ============================================================

build_D <- function(n) {

  if (n < 2L) {
    return(
      matrix(
        numeric(0),
        nrow = 0,
        ncol = n
      )
    )
  }

  D <- matrix(
    0,
    nrow = n - 1,
    ncol = n
  )

  for (i in seq_len(n - 1)) {
    D[i, i] <- -1
    D[i, i + 1] <- 1
  }

  D
}


ell_one_series <- function(
    xj,
    mu,
    M,
    r2) {

  N <- length(mu)

  if (
    length(xj) != N ||
      length(M) != N - 1L ||
      length(r2) != N ||
      any(!is.finite(r2)) ||
      any(r2 <= 0)
  ) {
    return(Inf)
  }

  value <- 0.5 * sum(
    log(r2)
  )

  value <- value +
    (xj[1] - mu[1])^2 /
      (2 * r2[1])

  if (N >= 2L) {

    residual <- xj[2:N] -
      mu[2:N] -
      M[1:(N - 1)] *
      (
        xj[1:(N - 1)] -
          mu[1:(N - 1)]
      )

    value <- value +
      sum(
        residual^2 /
          (2 * r2[2:N])
      )
  }

  value
}


compute_St <- function(
    X,
    mu,
    M) {

  N <- nrow(X)
  S <- numeric(N)

  S[1] <- mean(
    (
      X[1, ] -
        mu[1]
    )^2
  )

  if (N >= 2L) {
    for (t in 2:N) {

      residual_t <- X[t, ] -
        mu[t] -
        M[t - 1] *
        (
          X[t - 1, ] -
            mu[t - 1]
        )

      S[t] <- mean(
        residual_t^2
      )
    }
  }

  S
}


solve_r2_cubic <- function(
    S,
    neighbor_values,
    lambda3,
    J,
    x_init = NULL,
    eps = 1e-10) {

  d <- length(neighbor_values)
  q <- sum(neighbor_values)

  if (lambda3 == 0 || d == 0L) {
    return(
      max(S, eps)
    )
  }

  roots <- polyroot(
    c(
      -J * S,
      J,
      -4 * lambda3 * q,
      4 * lambda3 * d
    )
  )

  candidates <- Re(
    roots[
      abs(Im(roots)) < 1e-8
    ]
  )

  candidates <- candidates[
    is.finite(candidates) &
      candidates > eps
  ]

  coordinate_objective <- function(x) {

    likelihood_part <- (
      J / 2
    ) * (
      log(x) +
        S / x
    )

    penalty_part <- lambda3 *
      sum(
        (
          x -
            neighbor_values
        )^2
      )

    likelihood_part +
      penalty_part
  }

  if (length(candidates) > 0L) {

    values <- vapply(
      candidates,
      coordinate_objective,
      numeric(1)
    )

    return(
      candidates[
        which.min(values)
      ]
    )
  }

  center <- max(
    c(
      S,
      neighbor_values,
      x_init,
      eps
    ),
    na.rm = TRUE
  )

  optimization <- optimize(
    function(z) {
      coordinate_objective(
        exp(z)
      )
    },
    interval = c(
      log(center) - 20,
      log(center) + 20
    )
  )

  max(
    exp(optimization$minimum),
    eps
  )
}


update_M <- function(
    X,
    mu,
    r2,
    lambda1) {

  N <- nrow(X)
  D_M <- build_D(N - 1)

  A_diag <- numeric(N - 1)
  a <- numeric(N - 1)

  for (t in seq_len(N - 1)) {

    A_diag[t] <- sum(
      (
        X[t, ] -
          mu[t]
      )^2
    ) / r2[t + 1]

    a[t] <- sum(
      (
        X[t, ] -
          mu[t]
      ) *
        (
          X[t + 1, ] -
            mu[t + 1]
        )
    ) / r2[t + 1]
  }

  A <- diag(
    A_diag,
    nrow = N - 1
  )

  K <- A +
    2 * lambda1 *
    crossprod(D_M)

  as.vector(
    solve(K, a)
  )
}


update_mu <- function(
    X,
    M,
    r2,
    lambda2) {

  N <- nrow(X)
  J <- ncol(X)
  D_mu <- build_D(N)

  B <- matrix(
    0,
    nrow = N,
    ncol = N
  )

  b <- numeric(N)

  B[1, 1] <- B[1, 1] + J / r2[1]

  b[1] <- b[1] + sum(X[1, ]) / r2[1]

  for (t in 2:N) {

    weight <- J /
      r2[t]

    m <- M[t - 1]

    B[t, t] <- B[t, t] +
      weight

    B[t - 1, t - 1] <-
      B[t - 1, t - 1] +
      weight * m^2

    B[t, t - 1] <-
      B[t, t - 1] -
      weight * m

    B[t - 1, t] <-
      B[t - 1, t] -
      weight * m

    summed_term <- sum(
      X[t, ]
    ) -
      m *
      sum(
        X[t - 1, ]
      )

    b[t] <- b[t] +
      summed_term /
      r2[t]

    b[t - 1] <- b[t - 1] -
      m *
      summed_term /
      r2[t]
  }

  K <- B +
    2 * lambda2 *
    crossprod(D_mu)

  as.vector(
    solve(K, b)
  )
}


update_r2 <- function(
    X,
    mu,
    M,
    r2,
    lambda3) {

  N <- nrow(X)
  J <- ncol(X)

  S <- compute_St(
    X,
    mu,
    M
  )

  r2_new <- r2

  for (t in seq_len(N)) {

    neighbor_values <- numeric(0)

    if (t > 1L) {
      neighbor_values <- c(
        neighbor_values,
        r2_new[t - 1]
      )
    }

    if (t < N) {
      neighbor_values <- c(
        neighbor_values,
        r2[t + 1]
      )
    }

    r2_new[t] <- solve_r2_cubic(
      S = S[t],
      neighbor_values = neighbor_values,
      lambda3 = lambda3,
      J = J,
      x_init = r2[t]
    )

    if (
      !is.finite(r2_new[t]) ||
        r2_new[t] <= 0
    ) {
      r2_new[t] <- max(
        S[t],
        1e-10
      )
    }
  }

  r2_new
}


fit_theta_once <- function(
    X,
    lambda1,
    lambda2,
    lambda3,
    max_iter = 100,
    tol = 1e-6,
    verbose = FALSE) {

  stopifnot(
    is.matrix(X),
    all(is.finite(X))
  )

  N <- nrow(X)
  J <- ncol(X)

  mu <- as.vector(
    rowMeans(X)
  )

  r2 <- as.vector(
    apply(
      X,
      1,
      function(z) {
        mean(
          (
            z -
              mean(z)
          )^2
        )
      }
    )
  )

  r2 <- pmax(
    r2,
    1e-10
  )

  M <- rep(
    0,
    N - 1
  )

  objective <- function() {

    likelihood_sum <- 0

    for (j in seq_len(J)) {
      likelihood_sum <- likelihood_sum +
        ell_one_series(
          X[, j],
          mu,
          M,
          r2
        )
    }

    D_M <- build_D(N - 1)
    D_mu <- build_D(N)
    D_r2 <- build_D(N)

    penalty <- lambda1 *
      sum(
        (
          D_M %*% M
        )^2
      ) +
      lambda2 *
      sum(
        (
          D_mu %*% mu
        )^2
      ) +
      lambda3 *
      sum(
        (
          D_r2 %*% r2
        )^2
      )

    likelihood_sum +
      penalty
  }

  previous_objective <- objective()
  converged <- FALSE
  iteration <- 0L

  for (iteration in seq_len(max_iter)) {

    M <- update_M(
      X,
      mu,
      r2,
      lambda1
    )

    mu <- update_mu(
      X,
      M,
      r2,
      lambda2
    )

    r2 <- update_r2(
      X,
      mu,
      M,
      r2,
      lambda3
    )

    current_objective <- objective()

    if (verbose) {
      cat(
        sprintf(
          "Iteration %d: objective = %.6f\n",
          iteration,
          current_objective
        )
      )
    }

    if (
      abs(
        previous_objective -
          current_objective
      ) <
        tol *
        (
          1 +
            abs(previous_objective)
        )
    ) {
      converged <- TRUE
      previous_objective <- current_objective
      break
    }

    previous_objective <- current_objective
  }

  list(
    M = M,
    mu = mu,
    r2 = r2,
    obj = previous_objective,
    iters = iteration,
    converged = converged
  )
}


cv_select_lambdas <- function(
    X,
    lambda1_grid,
    lambda2_grid,
    lambda3_grid,
    max_iter,
    tol,
    verbose = FALSE) {

  J <- ncol(X)

  combinations <- expand.grid(
    lambda1 = lambda1_grid,
    lambda2 = lambda2_grid,
    lambda3 = lambda3_grid,
    KEEP.OUT.ATTRS = FALSE
  )

  combinations$CV <- NA_real_

  best <- list(
    score = Inf,
    lambda1 = NA_real_,
    lambda2 = NA_real_,
    lambda3 = NA_real_
  )

  for (i in seq_len(nrow(combinations))) {

    lambda1 <- combinations$lambda1[i]
    lambda2 <- combinations$lambda2[i]
    lambda3 <- combinations$lambda3[i]

    scores <- numeric(J)

    for (fold in seq_len(J)) {

      training_columns <- setdiff(
        seq_len(J),
        fold
      )

      X_train <- X[
        ,
        training_columns,
        drop = FALSE
      ]

      X_test <- X[
        ,
        fold
      ]

      fit <- fit_theta_once(
        X = X_train,
        lambda1 = lambda1,
        lambda2 = lambda2,
        lambda3 = lambda3,
        max_iter = max_iter,
        tol = tol,
        verbose = FALSE
      )

      scores[fold] <- ell_one_series(
        X_test,
        fit$mu,
        fit$M,
        fit$r2
      )
    }

    cv_score <- mean(scores)
    combinations$CV[i] <- cv_score

    if (verbose) {
      cat(
        sprintf(
          paste0(
            "Combination %d/%d: ",
            "lambda1=%.6g, ",
            "lambda2=%.6g, ",
            "lambda3=%.6g, ",
            "CV=%.6f\n"
          ),
          i,
          nrow(combinations),
          lambda1,
          lambda2,
          lambda3,
          cv_score
        )
      )
    }

    if (cv_score < best$score) {
      best <- list(
        score = cv_score,
        lambda1 = lambda1,
        lambda2 = lambda2,
        lambda3 = lambda3
      )
    }
  }

  list(
    best = best,
    table = combinations
  )
}


fit_LME_fusedridge <- function(
    X,
    lambda1_grid,
    lambda2_grid,
    lambda3_grid,
    max_iter,
    tol,
    verbose = FALSE) {

  stopifnot(
    is.matrix(X),
    ncol(X) == 13L
  )

  cv_result <- cv_select_lambdas(
    X = X,
    lambda1_grid = lambda1_grid,
    lambda2_grid = lambda2_grid,
    lambda3_grid = lambda3_grid,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose
  )

  best <- cv_result$best

  if (verbose) {
    cat(
      sprintf(
        paste0(
          "Best lambdas: ",
          "lambda1=%.6g, ",
          "lambda2=%.6g, ",
          "lambda3=%.6g; ",
          "CV=%.6f\n"
        ),
        best$lambda1,
        best$lambda2,
        best$lambda3,
        best$score
      )
    )
  }

  fit <- fit_theta_once(
    X = X,
    lambda1 = best$lambda1,
    lambda2 = best$lambda2,
    lambda3 = best$lambda3,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose
  )

  list(
    lambdas = c(
      lambda1 = best$lambda1,
      lambda2 = best$lambda2,
      lambda3 = best$lambda3
    ),
    CV = best$score,
    theta = list(
      M = fit$M,
      mu = fit$mu,
      r2 = fit$r2
    ),
    optimization = list(
      obj = fit$obj,
      iters = fit$iters,
      converged = fit$converged
    ),
    cv_table = cv_result$table
  )
}

# ============================================================
# Input and output functions
# ============================================================

read_city_lme <- function(
    lme_city_data,
    city_name) {

  city_data <- lme_city_data %>%
    filter(
      city == city_name
    ) %>%
    dplyr::select(
      member,
      year,
      temperature_celsius
    )

  if (nrow(city_data) == 0L) {
    stop(
      "No LME data were found for ",
      city_name,
      "."
    )
  }

  duplicate_rows <- city_data %>%
    count(
      member,
      year,
      name = "number_of_rows"
    ) %>%
    filter(
      number_of_rows != 1L
    )

  if (nrow(duplicate_rows) > 0L) {
    stop(
      "Each member-year combination must occur exactly once ",
      "for ",
      city_name,
      "."
    )
  }

  member_names <- sort(
    unique(
      city_data$member
    )
  )

  if (length(member_names) != 13L) {
    stop(
      city_name,
      " should have 13 ensemble members, but has ",
      length(member_names),
      "."
    )
  }

  missing_years <- setdiff(
    analysis_years,
    unique(
      city_data$year
    )
  )

  if (length(missing_years) > 0L) {
    stop(
      city_name,
      " is missing ",
      length(missing_years),
      " analysis years. First missing year: ",
      missing_years[1],
      "."
    )
  }

  wide_data <- city_data %>%
    filter(
      year %in%
        analysis_years
    ) %>%
    tidyr::pivot_wider(
      names_from = member,
      values_from = temperature_celsius
    ) %>%
    arrange(
      year
    )

  if (!identical(
    wide_data$year,
    as.integer(
      analysis_years
    )
  )) {
    stop(
      "The LME years are not aligned correctly for ",
      city_name,
      "."
    )
  }

  X <- as.matrix(
    wide_data[
      ,
      member_names,
      drop = FALSE
    ]
  )

  storage.mode(
    X
  ) <- "double"

  if (
    nrow(X) !=
      length(
        analysis_years
      ) ||
      ncol(X) != 13L ||
      any(!is.finite(X))
  ) {
    stop(
      "The LME matrix has invalid dimensions or values for ",
      city_name,
      "."
    )
  }

  rownames(
    X
  ) <- as.character(
    analysis_years
  )

  colnames(
    X
  ) <- member_names

  X
}


make_parameter_table <- function(
    years,
    penalized,
    unpenalized) {

  wide_table <- data.frame(
    year = years,
    `Penalized ML` = penalized,
    ML = unpenalized,
    check.names = FALSE
  )

  tidyr::pivot_longer(
    wide_table,
    cols = -year,
    names_to = "coefficient",
    values_to = "value"
  )
}


estimate_city <- function(
    city_name,
    lme_city_data,
    city_code) {

  message(
    "\nEstimating prior parameters for ",
    city_name,
    "..."
  )

  X <- read_city_lme(
    lme_city_data = lme_city_data,
    city_name = city_name
  )

  penalized_fit <- fit_LME_fusedridge(
    X = X,
    lambda1_grid = lambda1_grid,
    lambda2_grid = lambda2_grid,
    lambda3_grid = lambda3_grid,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose
  )

  unpenalized_fit <- fit_theta_once(
    X = X,
    lambda1 = 0,
    lambda2 = 0,
    lambda3 = 0,
    max_iter = max_iter,
    tol = tol,
    verbose = FALSE
  )

  mt_table <- make_parameter_table(
    years = analysis_years[
      -length(analysis_years)
    ],
    penalized = penalized_fit$theta$M,
    unpenalized = unpenalized_fit$M
  )

  mu_table <- make_parameter_table(
    years = analysis_years,
    penalized = penalized_fit$theta$mu,
    unpenalized = unpenalized_fit$mu
  )

  rt_table <- make_parameter_table(
    years = analysis_years,
    penalized = penalized_fit$theta$r2,
    unpenalized = unpenalized_fit$r2
  )

  readr::write_excel_csv(
    mt_table,
    file.path(
      output_dir,
      paste0("mt", city_code, ".csv")
    )
  )

  readr::write_excel_csv(
    mu_table,
    file.path(
      output_dir,
      paste0("mu", city_code, ".csv")
    )
  )

  readr::write_excel_csv(
    rt_table,
    file.path(
      output_dir,
      paste0("rt", city_code, ".csv")
    )
  )

  city_cv_table <- penalized_fit$cv_table %>%
    mutate(
      city = city_name,
      .before = 1
    ) %>%
    arrange(
      CV
    )

  readr::write_csv(
    city_cv_table,
    file.path(
      output_dir,
      paste0(
        "prior_cv_",
        city_name,
        ".csv"
      )
    )
  )

  saveRDS(
    list(
      city = city_name,
      input_file = lme_input_file,
      years = analysis_years,
      penalized = penalized_fit,
      unpenalized = unpenalized_fit
    ),
    file.path(
      output_dir,
      paste0(
        "prior_",
        city_name,
        ".rds"
      )
    )
  )

  selected_lambdas <- unname(
    penalized_fit$lambdas[
      c(
        "lambda1",
        "lambda2",
        "lambda3"
      )
    ]
  )

  lambda_on_boundary <- c(
    selected_lambdas[1] %in%
      range(
        lambda1_grid
      ),
    selected_lambdas[2] %in%
      range(
        lambda2_grid
      ),
    selected_lambdas[3] %in%
      range(
        lambda3_grid
      )
  )

  best_lambdas <- data.frame(
    city = city_name,
    lambda1 = selected_lambdas[1],
    lambda2 = selected_lambdas[2],
    lambda3 = selected_lambdas[3],
    lambda1_on_grid_boundary =
      lambda_on_boundary[1],
    lambda2_on_grid_boundary =
      lambda_on_boundary[2],
    lambda3_on_grid_boundary =
      lambda_on_boundary[3],
    CV = penalized_fit$CV,
    iterations =
      penalized_fit$optimization$iters,
    converged =
      penalized_fit$optimization$converged
  )

  if (any(
    lambda_on_boundary
  )) {
    warning(
      city_name,
      ": at least one selected lambda lies on the edge of ",
      "the search grid [",
      min(
        lambda1_grid
      ),
      ", ",
      max(
        lambda1_grid
      ),
      "]. Consider expanding or refining that grid before the ",
      "final manuscript run."
    )
  }

  message(
    "Completed ",
    city_name,
    "."
  )

  list(
    mt = mt_table,
    mu = mu_table,
    rt = rt_table,
    penalized = penalized_fit,
    unpenalized = unpenalized_fit,
    best_lambdas = best_lambdas
  )
}

# ============================================================
# Run Hong Kong, Shanghai, and Beijing automatically
# ============================================================

start_time <- Sys.time()

prior_results <- lapply(
  names(city_config),
  function(city_name) {

    config <- city_config[[city_name]]

    result <- estimate_city(
      city_name = city_name,
      lme_city_data = lme_city_data,
      city_code = config$code
    )

    gc()

    result
  }
)

names(prior_results) <- names(
  city_config
)

best_lambda_summary <- do.call(
  rbind,
  lapply(
    prior_results,
    function(x) {
      x$best_lambdas
    }
  )
)

rownames(best_lambda_summary) <- NULL

readr::write_csv(
  best_lambda_summary,
  file.path(
    output_dir,
    "prior_best_lambdas.csv"
  )
)

saveRDS(
  prior_results,
  file.path(
    output_dir,
    "prior_results_all_cities.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "Prior_sessionInfo.txt"
  )
)

elapsed_minutes <- as.numeric(
  difftime(
    Sys.time(),
    start_time,
    units = "mins"
  )
)

message(
  "\nPrior estimation completed for Hong Kong, Shanghai, and Beijing."
)

message(
  "Outputs saved to: ",
  output_dir
)

message(
  "Elapsed time: ",
  round(
    elapsed_minutes,
    2
  ),
  " minutes."
)
