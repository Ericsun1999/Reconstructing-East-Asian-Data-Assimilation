here::i_am("Code/Assimilation_grid.R")

# ============================================================
# Generate the spatial assimilated posterior required by
# Supplementary Figure S5.
#
# Default spatial support:
#   The exact LME-grid locations retained and clustered by
#   Code/Figure6.R. Set location_scope = "all_lme" to run all
#   native LME grid locations instead.
#
# Spatial-prior smoothing parameters:
#   lambda1 = 1000
#   lambda2 = 10^1.5
#   lambda3 = 1000
#
# These fixed values reproduce the original spatial analysis;
# location-specific cross-validation is intentionally not run
# because of the number of grid locations.
#
# Required generated inputs:
#   Output/Figure6/Figure6_cluster_assignments.csv
#   Output/Intermediate/LME/lme_annual_1368_1911.rds
#   Output/Intermediate/REACHES/
#     reaches_kriging_lme_grid_mean.csv
#     reaches_kriging_lme_grid_variance.csv
#     reaches_kriging_metadata.csv
#
# Generated outputs:
#   Output/Intermediate/Assimilation/
#     assimilated_posterior_lme_grid_mean.csv
#     assimilated_posterior_lme_grid_variance.csv
#     assimilated_posterior_lme_grid_diagnostics.csv
#     assimilated_posterior_lme_grid_metadata.csv
#     Assimilation_grid_sessionInfo.txt
#
# Optional repository validation copies:
#   Data/Valid/
#     assimilated_posterior_lme_grid_mean.csv
#     assimilated_posterior_lme_grid_variance.csv
#
# The filter and RTS smoother run on every year from 1368 to
# 1911. Years with no REACHES documentary information receive a
# prediction-only Kalman step.
# ============================================================

library(dplyr)
library(readr)

# ------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------

analysis_years <- 1368:1911

input_mode <- "auto"
location_scope <- "figure6"

allowed_input_modes <- c(
  "auto",
  "generated",
  "precomputed"
)

allowed_location_scopes <- c(
  "figure6",
  "all_lme"
)

if (!input_mode %in% allowed_input_modes) {
  stop(
    "input_mode must be one of: ",
    paste(allowed_input_modes, collapse = ", ")
  )
}

if (!location_scope %in% allowed_location_scopes) {
  stop(
    "location_scope must be one of: ",
    paste(allowed_location_scopes, collapse = ", ")
  )
}

fixed_lambda1 <- 1000
fixed_lambda2 <- 10^1.5
fixed_lambda3 <- 1000

prior_max_iter <- 100L
prior_tol <- 1e-6
prior_verbose <- FALSE

measurement_mc_size <- 10000L
measurement_mc_seed <- 1L
mapping_grid_size <- 2001L
lme_cdf_grid_size <- 5001L
latent_tail_probability <- 1e-14

variance_tolerance <- 1e-10
beta_tolerance <- 1e-10
identity_tolerance <- 1e-8
coordinate_tolerance <- 1e-8

resume_from_checkpoint <- TRUE
checkpoint_every <- 1L
remove_checkpoint_on_success <- TRUE
save_data_valid_copy <- TRUE

output_dir <- here::here(
  "Output",
  "Intermediate",
  "Assimilation"
)

validation_dir <- here::here(
  "Data",
  "Valid"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  validation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_files <- list(
  mean = file.path(
    output_dir,
    "assimilated_posterior_lme_grid_mean.csv"
  ),
  variance = file.path(
    output_dir,
    "assimilated_posterior_lme_grid_variance.csv"
  ),
  diagnostics = file.path(
    output_dir,
    "assimilated_posterior_lme_grid_diagnostics.csv"
  ),
  metadata = file.path(
    output_dir,
    "assimilated_posterior_lme_grid_metadata.csv"
  ),
  checkpoint = file.path(
    output_dir,
    "assimilated_posterior_lme_grid_checkpoint.rds"
  ),
  session_info = file.path(
    output_dir,
    "Assimilation_grid_sessionInfo.txt"
  )
)

validation_files <- list(
  mean = file.path(
    validation_dir,
    "assimilated_posterior_lme_grid_mean.csv"
  ),
  variance = file.path(
    validation_dir,
    "assimilated_posterior_lme_grid_variance.csv"
  )
)

# ------------------------------------------------------------
# 2. Input selection
# ------------------------------------------------------------

generated_inputs <- list(
  clusters = here::here(
    "Output",
    "Figure6",
    "Figure6_cluster_assignments.csv"
  ),
  lme_archive = here::here(
    "Output",
    "Intermediate",
    "LME",
    "lme_annual_1368_1911.rds"
  ),
  reaches_mean = here::here(
    "Output",
    "Intermediate",
    "REACHES",
    "reaches_kriging_lme_grid_mean.csv"
  ),
  reaches_variance = here::here(
    "Output",
    "Intermediate",
    "REACHES",
    "reaches_kriging_lme_grid_variance.csv"
  ),
  reaches_metadata = here::here(
    "Output",
    "Intermediate",
    "REACHES",
    "reaches_kriging_metadata.csv"
  )
)

precomputed_inputs <- list(
  clusters = here::here(
    "Data",
    "REACHES",
    "precomputed",
    "Figure6_cluster_assignments.csv"
  ),
  lme_archive = here::here(
    "Data",
    "LME data",
    "precomputed",
    "lme_annual_1368_1911.rds"
  ),
  reaches_mean = here::here(
    "Data",
    "REACHES",
    "precomputed",
    "reaches_kriging_lme_grid_mean.csv"
  ),
  reaches_variance = here::here(
    "Data",
    "REACHES",
    "precomputed",
    "reaches_kriging_lme_grid_variance.csv"
  ),
  reaches_metadata = here::here(
    "Data",
    "REACHES",
    "precomputed",
    "reaches_kriging_metadata.csv"
  )
)

required_input_names <- function(location_scope) {
  required <- c(
    "lme_archive",
    "reaches_mean",
    "reaches_variance",
    "reaches_metadata"
  )

  if (location_scope == "figure6") {
    required <- c("clusters", required)
  }

  required
}

input_set_is_complete <- function(
    input_set,
    required_names) {

  all(
    file.exists(
      unlist(
        input_set[required_names],
        use.names = FALSE
      )
    )
  )
}

missing_input_files <- function(
    input_set,
    required_names) {

  paths <- unlist(
    input_set[required_names],
    use.names = FALSE
  )

  paths[!file.exists(paths)]
}

select_input_set <- function(
    input_mode,
    location_scope,
    generated_inputs,
    precomputed_inputs) {

  required_names <- required_input_names(
    location_scope
  )

  generated_complete <- input_set_is_complete(
    generated_inputs,
    required_names
  )

  precomputed_complete <- input_set_is_complete(
    precomputed_inputs,
    required_names
  )

  if (input_mode == "generated") {
    if (!generated_complete) {
      stop(
        "The generated spatial-assimilation input set is incomplete. Missing:\n  ",
        paste(
          missing_input_files(
            generated_inputs,
            required_names
          ),
          collapse = "\n  "
        )
      )
    }

    return(generated_inputs)
  }

  if (input_mode == "precomputed") {
    if (!precomputed_complete) {
      stop(
        "The precomputed spatial-assimilation input set is incomplete. Missing:\n  ",
        paste(
          missing_input_files(
            precomputed_inputs,
            required_names
          ),
          collapse = "\n  "
        )
      )
    }

    return(precomputed_inputs)
  }

  if (generated_complete) {
    return(generated_inputs)
  }

  if (precomputed_complete) {
    return(precomputed_inputs)
  }

  stop(
    "Neither a complete generated nor a complete precomputed ",
    "spatial-assimilation input set was found."
  )
}

input_files <- select_input_set(
  input_mode = input_mode,
  location_scope = location_scope,
  generated_inputs = generated_inputs,
  precomputed_inputs = precomputed_inputs
)

message(
  "Spatial assimilation location scope: ",
  location_scope
)

message(
  "LME archive selected: ",
  input_files$lme_archive
)

# ------------------------------------------------------------
# 3. Input standardization and alignment
# ------------------------------------------------------------

coordinate_key <- function(
    long,
    lat) {

  paste(
    sprintf("%.8f", as.numeric(long)),
    sprintf("%.8f", as.numeric(lat)),
    sep = "_"
  )
}

standardize_grid_file <- function(
    input_file,
    object_name) {

  data <- readr::read_csv(
    input_file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  if (
    !"lati" %in% names(data) &&
      "lat" %in% names(data)
  ) {
    data <- data %>%
      rename(lati = lat)
  }

  required_coordinates <- c(
    "long",
    "lati"
  )

  missing_coordinates <- setdiff(
    required_coordinates,
    names(data)
  )

  if (length(missing_coordinates) > 0L) {
    stop(
      object_name,
      " is missing coordinate columns: ",
      paste(missing_coordinates, collapse = ", ")
    )
  }

  raw_year_columns <- grep(
    "^[Xx]?[0-9]{4}$",
    names(data),
    value = TRUE
  )

  if (length(raw_year_columns) == 0L) {
    stop(
      object_name,
      " contains no annual columns."
    )
  }

  raw_year_values <- as.integer(
    sub("^[Xx]", "", raw_year_columns)
  )

  names(data)[match(raw_year_columns, names(data))] <- paste0("x", raw_year_values)

  data <- data %>%
    mutate(
      location_id = if (
        "location_id" %in% names(data)
      ) {
        as.integer(location_id)
      } else {
        NA_integer_
      },
      long = as.numeric(long),
      lati = as.numeric(lati),
      cell_id = coordinate_key(long, lati)
    )

  if (
    any(!is.finite(data$long)) ||
      any(!is.finite(data$lati)) ||
      anyDuplicated(data$cell_id)
  ) {
    stop(
      object_name,
      " contains invalid or duplicated coordinates."
    )
  }

  data
}

extract_grid_years <- function(data) {
  columns <- grep(
    "^x[0-9]{4}$",
    names(data),
    value = TRUE
  )

  sort(
    as.integer(
      sub("^x", "", columns)
    )
  )
}

lme_archive <- readRDS(
  input_files$lme_archive
)

required_archive_elements <- c(
  "coordinates",
  "years",
  "members",
  "annual_kelvin"
)

missing_archive_elements <- setdiff(
  required_archive_elements,
  names(lme_archive)
)

if (length(missing_archive_elements) > 0L) {
  stop(
    "The LME annual archive is missing: ",
    paste(missing_archive_elements, collapse = ", ")
  )
}

lme_coordinates <- as.data.frame(
  lme_archive$coordinates
)

if (
  !"lati" %in% names(lme_coordinates) &&
    "lat" %in% names(lme_coordinates)
) {
  lme_coordinates <- lme_coordinates %>%
    rename(lati = lat)
}

required_lme_coordinate_columns <- c(
  "location_id",
  "long",
  "lati"
)

if (!all(
  required_lme_coordinate_columns %in%
    names(lme_coordinates)
)) {
  stop(
    "The LME archive coordinates must contain location_id, ",
    "long, and lati."
  )
}

lme_coordinates <- lme_coordinates %>%
  transmute(
    location_id = as.integer(location_id),
    long = as.numeric(long),
    lati = as.numeric(lati),
    cell_id = coordinate_key(long, lati)
  )

if (
  anyNA(lme_coordinates) ||
    anyDuplicated(lme_coordinates$location_id) ||
    anyDuplicated(lme_coordinates$cell_id)
) {
  stop(
    "The LME archive contains invalid or duplicated coordinates."
  )
}

lme_years <- as.integer(
  lme_archive$years
)

analysis_year_indices <- match(
  analysis_years,
  lme_years
)

if (anyNA(analysis_year_indices)) {
  stop(
    "The LME archive does not contain every analysis year from ",
    min(analysis_years),
    " to ",
    max(analysis_years),
    "."
  )
}

annual_kelvin <- lme_archive$annual_kelvin

if (
  length(dim(annual_kelvin)) != 3L ||
    dim(annual_kelvin)[1] != nrow(lme_coordinates) ||
    dim(annual_kelvin)[2] != length(lme_years) ||
    dim(annual_kelvin)[3] != 13L
) {
  stop(
    "The LME annual array must have dimensions locations x ",
    "years x 13 members."
  )
}

reaches_mean <- standardize_grid_file(
  input_files$reaches_mean,
  "REACHES mean grid"
)

reaches_variance <- standardize_grid_file(
  input_files$reaches_variance,
  "REACHES variance grid"
)

mean_event_years <- extract_grid_years(
  reaches_mean
)

variance_event_years <- extract_grid_years(
  reaches_variance
)

event_years <- sort(
  Reduce(
    intersect,
    list(
      mean_event_years,
      variance_event_years,
      analysis_years
    )
  )
)

if (length(event_years) < 2L) {
  stop(
    "Fewer than two common REACHES event years were found."
  )
}

event_columns <- paste0(
  "x",
  event_years
)

event_year_indices <- match(
  event_years,
  analysis_years
)

reaches_metadata <- readr::read_csv(
  input_files$reaches_metadata,
  show_col_types = FALSE
)

if (
  nrow(reaches_metadata) != 1L ||
    !"process_variance" %in% names(reaches_metadata)
) {
  stop(
    "The REACHES metadata must contain one process_variance value."
  )
}

sigmaY2 <- as.numeric(
  reaches_metadata$process_variance[1]
)

if (
  length(sigmaY2) != 1L ||
    !is.finite(sigmaY2) ||
    sigmaY2 <= 0
) {
  stop(
    "The REACHES process variance is invalid."
  )
}

read_figure6_locations <- function(input_file) {

  assignments <- readr::read_csv(
    input_file,
    show_col_types = FALSE
  )

  if (
    !"lati" %in% names(assignments) &&
      "lat" %in% names(assignments)
  ) {
    assignments <- assignments %>%
      rename(lati = lat)
  }

  required <- c(
    "long",
    "lati",
    "figure6_cluster"
  )

  missing <- setdiff(
    required,
    names(assignments)
  )

  if (length(missing) > 0L) {
    stop(
      "Figure6_cluster_assignments.csv is missing: ",
      paste(missing, collapse = ", ")
    )
  }

  assignments %>%
    transmute(
      location_id = if (
        "location_id" %in% names(assignments)
      ) {
        as.integer(location_id)
      } else {
        NA_integer_
      },
      long = as.numeric(long),
      lati = as.numeric(lati),
      figure6_cluster = as.integer(figure6_cluster),
      cell_id = coordinate_key(long, lati)
    ) %>%
    arrange(
      figure6_cluster,
      location_id,
      lati,
      long
    )
}

if (location_scope == "figure6") {
  selected_locations <- read_figure6_locations(
    input_files$clusters
  )
} else {
  selected_locations <- lme_coordinates %>%
    mutate(
      figure6_cluster = NA_integer_
    ) %>%
    dplyr::select(
      location_id,
      long,
      lati,
      figure6_cluster,
      cell_id
    )
}

if (
  any(!is.finite(selected_locations$long)) ||
    any(!is.finite(selected_locations$lati)) ||
    anyDuplicated(selected_locations$cell_id)
) {
  stop(
    "The selected spatial locations are invalid or duplicated."
  )
}

match_locations <- function(
    selected_locations,
    grid_data,
    object_name) {

  use_location_id <-
    all(!is.na(selected_locations$location_id)) &&
    all(!is.na(grid_data$location_id)) &&
    !anyDuplicated(selected_locations$location_id) &&
    !anyDuplicated(grid_data$location_id)

  if (use_location_id) {
    matched_rows <- match(
      selected_locations$location_id,
      grid_data$location_id
    )
  } else {
    matched_rows <- match(
      selected_locations$cell_id,
      grid_data$cell_id
    )
  }

  if (anyNA(matched_rows)) {
    stop(
      sum(is.na(matched_rows)),
      " selected locations are absent from ",
      object_name,
      "."
    )
  }

  aligned <- grid_data[
    matched_rows,
    ,
    drop = FALSE
  ]

  coordinate_difference <- pmax(
    abs(aligned$long - selected_locations$long),
    abs(aligned$lati - selected_locations$lati)
  )

  if (any(
    coordinate_difference > coordinate_tolerance
  )) {
    stop(
      object_name,
      " does not align exactly with the selected coordinates."
    )
  }

  aligned
}

lme_rows <- match(
  selected_locations$location_id,
  lme_coordinates$location_id
)

if (anyNA(lme_rows)) {
  # Coordinate fallback for a precomputed assignment file that
  # lacks usable location identifiers.
  lme_rows <- match(
    selected_locations$cell_id,
    lme_coordinates$cell_id
  )
}

if (anyNA(lme_rows)) {
  stop(
    "At least one selected location is absent from the LME archive."
  )
}

selected_locations$location_id <-
  lme_coordinates$location_id[lme_rows]

reaches_mean_aligned <- match_locations(
  selected_locations,
  reaches_mean,
  "the REACHES mean grid"
)

reaches_variance_aligned <- match_locations(
  selected_locations,
  reaches_variance,
  "the REACHES variance grid"
)

message(
  "Number of spatial locations to assimilate: ",
  nrow(selected_locations)
)

message(
  "Number of REACHES event years: ",
  length(event_years)
)

message(
  "Fixed lambdas: lambda1 = ",
  fixed_lambda1,
  ", lambda2 = ",
  signif(fixed_lambda2, 8),
  ", lambda3 = ",
  fixed_lambda3
)

# ------------------------------------------------------------
# 4. Penalized annual-prior model functions
# ------------------------------------------------------------

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

  B[1, 1] <- B[1, 1] +
    J / r2[1]

  b[1] <- b[1] +
    sum(X[1, ]) /
      r2[1]

  for (t in 2:N) {

    weight <- J /
      r2[t]

    m <- M[
      t - 1
    ]

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


# ------------------------------------------------------------
# 5. Quantile mapping, measurement model, and Kalman smoother
# ------------------------------------------------------------

build_quantile_mapping <- function(
    lme_sample,
    yhat,
    nu,
    sigmaY2,
    bandwidth_seed,
    mapping_grid_size,
    lme_cdf_grid_size,
    latent_tail_probability) {

  lme_sample <- as.numeric(
    lme_sample
  )

  yhat <- as.numeric(
    yhat
  )

  nu <- as.numeric(
    nu
  )

  if (
    length(
      yhat
    ) !=
      length(
        nu
      ) ||
      any(!is.finite(
        lme_sample
      )) ||
      any(!is.finite(
        yhat
      )) ||
      any(!is.finite(
        nu
      )) ||
      any(
        nu < 0
      )
  ) {
    stop(
      "Invalid inputs were supplied to the quantile-mapping ",
      "construction."
    )
  }

  lme_sd <- stats::sd(
    lme_sample
  )

  if (
    !is.finite(
      lme_sd
    ) ||
      lme_sd <= 0
  ) {
    stop(
      "The LME calibration sample has zero or invalid ",
      "variation."
    )
  }

  set.seed(
    bandwidth_seed
  )

  fx_bandwidth <- np::npudistbw(
    dat = data.frame(
      x = lme_sample
    )
  )

  lme_grid <- seq(
    min(
      lme_sample
    ) -
      8 *
        lme_sd,
    max(
      lme_sample
    ) +
      8 *
        lme_sd,
    length.out = lme_cdf_grid_size
  )

  fx_grid <- as.numeric(
    fitted(
      np::npudist(
        bws = fx_bandwidth,
        edat = data.frame(
          x = lme_grid
        )
      )
    )
  )

  # Numerical monotonicity protection for the estimated CDF.
  fx_grid <- cummax(
    pmin(
      pmax(
        fx_grid,
        0
      ),
      1
    )
  )

  keep_unique_probability <- !duplicated(
    fx_grid
  )

  fx_probability_unique <- fx_grid[
    keep_unique_probability
  ]

  lme_grid_unique <- lme_grid[
    keep_unique_probability
  ]

  if (length(
    fx_probability_unique
  ) < 2L) {
    stop(
      "The estimated LME CDF could not be inverted."
    )
  }

  Fx_inv <- stats::approxfun(
    x = fx_probability_unique,
    y = lme_grid_unique,
    rule = 2,
    ties = "ordered"
  )

  nu_safe <- pmax(
    nu,
    1e-8
  )

  FY_hat <- function(y) {

    y <- as.numeric(
      y
    )

    standardized <- outer(
      y,
      yhat,
      "-"
    )

    standardized <- sweep(
      standardized,
      2,
      nu_safe,
      "/"
    )

    rowMeans(
      pnorm(
        standardized
      )
    )
  }

  latent_sd <- sqrt(
    sigmaY2
  )

  gaussian_limit <- qnorm(
    1 -
      latent_tail_probability /
        2
  ) *
    latent_sd

  proxy_limit_lower <- min(
    yhat -
      8 *
        nu_safe
  )

  proxy_limit_upper <- max(
    yhat +
      8 *
        nu_safe
  )

  latent_lower <- min(
    -gaussian_limit,
    proxy_limit_lower
  )

  latent_upper <- max(
    gaussian_limit,
    proxy_limit_upper
  )

  latent_grid <- seq(
    latent_lower,
    latent_upper,
    length.out = mapping_grid_size
  )

  mapping_probability <- FY_hat(
    latent_grid
  )

  mapped_grid <- Fx_inv(
    pmin(
      pmax(
        mapping_probability,
        1e-12
      ),
      1 -
        1e-12
    )
  )

  g <- stats::approxfun(
    x = latent_grid,
    y = mapped_grid,
    rule = 2,
    ties = "ordered"
  )

  corrected <- g(
    yhat
  )

  list(
    corrected = corrected,
    g = g,
    latent_range = range(
      latent_grid
    ),
    fx_bandwidth = as.numeric(
      fx_bandwidth$bw
    ),
    lme_range = range(
      lme_sample
    )
  )
}

# ------------------------------------------------------------
# 8. Measurement-equation Monte Carlo
# ------------------------------------------------------------

compute_measurement_parameters_mc <- function(
    sigmaY2,
    vhat,
    g,
    n_mc,
    seed) {

  vhat <- as.numeric(
    vhat
  )

  if (
    any(!is.finite(
      vhat
    )) ||
      any(
        vhat <
          -variance_tolerance
      ) ||
      any(
        vhat >
          sigmaY2 +
            variance_tolerance
      )
  ) {
    stop(
      "All Var(Yhat_t) values must lie in [0, sigmaY2]."
    )
  }

  vhat <- pmin(
    pmax(
      vhat,
      0
    ),
    sigmaY2
  )

  set.seed(
    seed
  )

  number_of_years <- length(
    vhat
  )

  alpha <- numeric(
    number_of_years
  )

  beta <- numeric(
    number_of_years
  )

  vdelta <- numeric(
    number_of_years
  )

  maximum_absolute_latent_draw <- 0

  for (
    t in seq_len(
      number_of_years
    )
  ) {

    current_vhat <- vhat[t]

    # This exact representation gives:
    #   Var(Yhat) = current_vhat
    #   Var(Y)    = sigmaY2
    #   Cov(Y, Yhat) = current_vhat
    #
    # No rounding or clipping is applied because Y and Yhat are
    # continuous latent Gaussian variables.
    Yhat <- sqrt(
      current_vhat
    ) *
      rnorm(
        n_mc
      )

    Y <- Yhat +
      sqrt(
        pmax(
          sigmaY2 -
            current_vhat,
          0
        )
      ) *
      rnorm(
        n_mc
      )

    maximum_absolute_latent_draw <- max(
      maximum_absolute_latent_draw,
      abs(
        Y
      ),
      abs(
        Yhat
      )
    )

    X <- g(
      Y
    )

    Xstar <- g(
      Yhat
    )

    variance_X <- stats::var(
      X
    )

    if (
      !is.finite(
        variance_X
      ) ||
        variance_X <= 0
    ) {
      stop(
        "A non-positive Monte Carlo variance was obtained at ",
        "measurement year index ",
        t,
        "."
      )
    }

    beta_t <- stats::cov(
      Xstar,
      X
    ) /
      variance_X

    alpha_t <- mean(
      Xstar
    ) -
      beta_t *
        mean(
          X
        )

    delta <- Xstar -
      alpha_t -
      beta_t *
        X

    alpha[t] <- alpha_t

    beta[t] <- beta_t

    vdelta[t] <- pmax(
      stats::var(
        delta
      ),
      1e-10
    )
  }

  list(
    alpha = alpha,
    beta = beta,
    vdelta = vdelta,
    maximum_absolute_latent_draw =
      maximum_absolute_latent_draw
  )
}

# ------------------------------------------------------------
# 9. Full-annual Kalman filter and RTS smoother
# ------------------------------------------------------------

kalman_filter_smoother_annual <- function(
    mu,
    M,
    r2,
    Xstar,
    alpha,
    beta,
    vdelta) {

  number_of_years <- length(
    mu
  )

  if (
    length(
      M
    ) !=
      number_of_years -
        1L ||
      length(
        r2
      ) !=
        number_of_years ||
      length(
        Xstar
      ) !=
        number_of_years ||
      length(
        alpha
      ) !=
        number_of_years ||
      length(
        beta
      ) !=
        number_of_years ||
      length(
        vdelta
      ) !=
        number_of_years
  ) {
    stop(
      "The state and observation vectors are not aligned to the ",
      "same annual grid."
    )
  }

  has_observation <- is.finite(
    Xstar
  )

  if (
    any(
      has_observation &
        (
          !is.finite(
            alpha
          ) |
            !is.finite(
              beta
            ) |
            !is.finite(
              vdelta
            ) |
            vdelta <= 0
        )
    )
  ) {
    stop(
      "At least one observed year has invalid measurement ",
      "parameters."
    )
  }

  X_pred <- numeric(
    number_of_years
  )

  P_pred <- numeric(
    number_of_years
  )

  X_filt <- numeric(
    number_of_years
  )

  P_filt <- numeric(
    number_of_years
  )

  kalman_gain <- rep(
    NA_real_,
    number_of_years
  )

  effective_proxy <- rep(
    NA_real_,
    number_of_years
  )

  observation_weight <- rep(
    NA_real_,
    number_of_years
  )

  innovation <- rep(
    NA_real_,
    number_of_years
  )

  X_pred[
    1
  ] <- mu[
    1
  ]

  P_pred[
    1
  ] <- r2[
    1
  ]

  for (
    t in seq_len(
      number_of_years
    )
  ) {

    if (t >= 2L) {

      transition <- M[t - 1L]

      X_pred[
        t
      ] <- mu[
        t
      ] +
        transition *
          (
            X_filt[
              t -
                1L
            ] -
              mu[
                t -
                  1L
              ]
          )

      P_pred[
        t
      ] <- transition^2 *
        P_filt[
          t -
            1L
        ] +
        r2[
          t
        ]
    }

    if (!has_observation[
      t
    ]) {

      # Prediction-only step for a year without documentary
      # information.
      X_filt[
        t
      ] <- X_pred[
        t
      ]

      P_filt[
        t
      ] <- P_pred[
        t
      ]

      next
    }

    denominator <- beta[
      t
    ]^2 *
      P_pred[
        t
      ] +
      vdelta[
        t
      ]

    if (
      !is.finite(
        denominator
      ) ||
        denominator <= 0
    ) {
      stop(
        "A non-positive Kalman update denominator was obtained ",
        "at year index ",
        t,
        "."
      )
    }

    kalman_gain[
      t
    ] <- beta[
      t
    ] *
      P_pred[
        t
      ] /
      denominator

    innovation[
      t
    ] <- Xstar[
      t
    ] -
      alpha[
        t
      ] -
      beta[
        t
      ] *
        X_pred[
          t
        ]

    X_filt[
      t
    ] <- X_pred[
      t
    ] +
      kalman_gain[
        t
      ] *
        innovation[
          t
        ]

    P_filt[
      t
    ] <- (
      1 -
        beta[
          t
        ] *
          kalman_gain[
            t
          ]
    ) *
      P_pred[
        t
      ]

    P_filt[
      t
    ] <- pmax(
      P_filt[
        t
      ],
      0
    )

    if (
      abs(
        beta[
          t
        ]
      ) >
        beta_tolerance
    ) {
      effective_proxy[
        t
      ] <- (
        Xstar[
          t
        ] -
          alpha[
            t
          ]
      ) /
        beta[
          t
        ]

      observation_weight[
        t
      ] <- beta[
        t
      ]^2 *
        P_pred[
          t
        ] /
        denominator
    }
  }

  X_smooth <- numeric(
    number_of_years
  )

  P_smooth <- numeric(
    number_of_years
  )

  smoother_gain <- numeric(
    number_of_years -
      1L
  )

  X_smooth[
    number_of_years
  ] <- X_filt[
    number_of_years
  ]

  P_smooth[
    number_of_years
  ] <- P_filt[
    number_of_years
  ]

  if (number_of_years >= 2L) {

    for (
      t in seq(
        from = number_of_years -
          1L,
        to = 1L,
        by = -1L
      )
    ) {

      if (
        !is.finite(
          P_pred[
            t +
              1L
          ]
        ) ||
          P_pred[
            t +
              1L
          ] <= 0
      ) {
        stop(
          "A non-positive predicted variance was obtained before ",
          "RTS smoothing at year index ",
          t +
            1L,
          "."
        )
      }

      smoother_gain[
        t
      ] <- P_filt[
        t
      ] *
        M[
          t
        ] /
        P_pred[
          t +
            1L
        ]

      X_smooth[
        t
      ] <- X_filt[
        t
      ] +
        smoother_gain[
          t
        ] *
          (
            X_smooth[
              t +
                1L
            ] -
              X_pred[
                t +
                  1L
              ]
          )

      P_smooth[
        t
      ] <- P_filt[
        t
      ] +
        smoother_gain[
          t
        ]^2 *
          (
            P_smooth[
              t +
                1L
            ] -
              P_pred[
                t +
                  1L
              ]
          )
    }
  }

  P_smooth <- pmax(
    P_smooth,
    0
  )

  filtered_check <- rep(
    NA_real_,
    number_of_years
  )

  check_years <- which(
    has_observation &
      is.finite(
        effective_proxy
      ) &
      is.finite(
        observation_weight
      )
  )

  filtered_check[
    check_years
  ] <- (
    1 -
      observation_weight[
        check_years
      ]
  ) *
    X_pred[
      check_years
    ] +
    observation_weight[
      check_years
    ] *
    effective_proxy[
      check_years
    ]

  maximum_filter_identity_error <- if (
    length(check_years) > 0L
  ) {
    max(
      abs(
        X_filt[check_years] -
          filtered_check[check_years]
      ),
      na.rm = TRUE
    )
  } else {
    0
  }

  if (
    length(
      check_years
    ) > 0L &&
      (
        !is.finite(
          maximum_filter_identity_error
        ) ||
          maximum_filter_identity_error >
            identity_tolerance
      )
  ) {
    stop(
      "The filtered-mean weighted-average identity failed. ",
      "Maximum absolute error: ",
      maximum_filter_identity_error
    )
  }

  list(
    has_observation = has_observation,
    X_pred = X_pred,
    P_pred = P_pred,
    X_filt = X_filt,
    P_filt = P_filt,
    X_smooth = X_smooth,
    P_smooth = P_smooth,
    kalman_gain = kalman_gain,
    smoother_gain = smoother_gain,
    innovation = innovation,
    effective_proxy = effective_proxy,
    observation_weight = observation_weight,
    filtered_check = filtered_check,
    maximum_filter_identity_error =
      maximum_filter_identity_error
  )
}

# ------------------------------------------------------------
# 7. Checkpoint helpers
# ------------------------------------------------------------

number_of_locations <- nrow(
  selected_locations
)

number_of_years <- length(
  analysis_years
)

posterior_mean_matrix <- matrix(
  NA_real_,
  nrow = number_of_locations,
  ncol = number_of_years
)

posterior_variance_matrix <- matrix(
  NA_real_,
  nrow = number_of_locations,
  ncol = number_of_years
)

completed_locations <- rep(
  FALSE,
  number_of_locations
)

diagnostic_rows <- vector(
  "list",
  number_of_locations
)

checkpoint_settings <- list(
  location_scope = location_scope,
  location_ids = selected_locations$location_id,
  analysis_years = analysis_years,
  event_years = event_years,
  lambda1 = fixed_lambda1,
  lambda2 = fixed_lambda2,
  lambda3 = fixed_lambda3,
  prior_max_iter = prior_max_iter,
  prior_tol = prior_tol,
  measurement_mc_size = measurement_mc_size,
  sigmaY2 = sigmaY2
)

checkpoint_is_compatible <- function(
    checkpoint,
    settings) {

  is.list(checkpoint) &&
    identical(
      checkpoint$settings,
      settings
    ) &&
    identical(
      dim(checkpoint$posterior_mean_matrix),
      c(
        number_of_locations,
        number_of_years
      )
    ) &&
    identical(
      dim(checkpoint$posterior_variance_matrix),
      c(
        number_of_locations,
        number_of_years
      )
    )
}

save_checkpoint <- function() {

  saveRDS(
    list(
      settings = checkpoint_settings,
      posterior_mean_matrix = posterior_mean_matrix,
      posterior_variance_matrix = posterior_variance_matrix,
      completed_locations = completed_locations,
      diagnostic_rows = diagnostic_rows
    ),
    output_files$checkpoint,
    compress = FALSE
  )
}

if (
  resume_from_checkpoint &&
    file.exists(output_files$checkpoint)
) {

  checkpoint <- readRDS(
    output_files$checkpoint
  )

  if (!checkpoint_is_compatible(
    checkpoint,
    checkpoint_settings
  )) {
    stop(
      "The existing spatial-assimilation checkpoint is not ",
      "compatible with the current settings. Delete:\n  ",
      output_files$checkpoint,
      "\nor set resume_from_checkpoint <- FALSE."
    )
  }

  posterior_mean_matrix <-
    checkpoint$posterior_mean_matrix

  posterior_variance_matrix <-
    checkpoint$posterior_variance_matrix

  completed_locations <-
    checkpoint$completed_locations

  diagnostic_rows <-
    checkpoint$diagnostic_rows

  message(
    "Resuming from checkpoint with ",
    sum(completed_locations),
    " of ",
    number_of_locations,
    " locations completed."
  )
}

# ------------------------------------------------------------
# 8. Process one spatial location
# ------------------------------------------------------------

process_location <- function(location_index) {

  location_id <- selected_locations$location_id[
    location_index
  ]

  location_long <- selected_locations$long[
    location_index
  ]

  location_lat <- selected_locations$lati[
    location_index
  ]

  lme_row <- lme_rows[
    location_index
  ]

  X <- annual_kelvin[
    lme_row,
    analysis_year_indices,
    ,
    drop = TRUE
  ]

  if (is.null(dim(X))) {
    X <- matrix(
      X,
      nrow = number_of_years,
      ncol = 13L
    )
  }

  if (!identical(
    dim(X),
    c(number_of_years, 13L)
  )) {
    X <- matrix(
      as.numeric(X),
      nrow = number_of_years,
      ncol = 13L
    )
  }

  X <- X - 273.15
  storage.mode(X) <- "double"

  if (any(!is.finite(X))) {
    stop(
      "Non-finite LME values were found at location_id ",
      location_id,
      "."
    )
  }

  prior_fit <- fit_theta_once(
    X = X,
    lambda1 = fixed_lambda1,
    lambda2 = fixed_lambda2,
    lambda3 = fixed_lambda3,
    max_iter = prior_max_iter,
    tol = prior_tol,
    verbose = prior_verbose
  )

  yhat <- as.numeric(
    reaches_mean_aligned[
      location_index,
      event_columns,
      drop = TRUE
    ]
  )

  mspe <- as.numeric(
    reaches_variance_aligned[
      location_index,
      event_columns,
      drop = TRUE
    ]
  )

  if (
    any(!is.finite(yhat)) ||
      any(!is.finite(mspe)) ||
      any(mspe < -variance_tolerance)
  ) {
    stop(
      "Invalid REACHES means or prediction variances were found ",
      "at location_id ",
      location_id,
      "."
    )
  }

  mspe <- pmax(
    mspe,
    0
  )

  nu <- sqrt(mspe)

  quantile_mapping <- build_quantile_mapping(
    lme_sample = as.numeric(X),
    yhat = yhat,
    nu = nu,
    sigmaY2 = sigmaY2,
    bandwidth_seed = 10000L + location_id,
    mapping_grid_size = mapping_grid_size,
    lme_cdf_grid_size = lme_cdf_grid_size,
    latent_tail_probability = latent_tail_probability
  )

  implied_vhat_raw <- sigmaY2 - mspe

  implied_vhat <- pmin(
    pmax(implied_vhat_raw, 0),
    sigmaY2
  )

  number_vhat_clipped <- sum(
    abs(
      implied_vhat -
        implied_vhat_raw
    ) > variance_tolerance
  )

  measurement_parameters <-
    compute_measurement_parameters_mc(
      sigmaY2 = sigmaY2,
      vhat = implied_vhat,
      g = quantile_mapping$g,
      n_mc = measurement_mc_size,
      seed = measurement_mc_seed +
        100003L * location_id
    )

  Xstar <- rep(
    NA_real_,
    number_of_years
  )

  alpha <- rep(
    NA_real_,
    number_of_years
  )

  beta <- rep(
    NA_real_,
    number_of_years
  )

  vdelta <- rep(
    NA_real_,
    number_of_years
  )

  Xstar[event_year_indices] <-
    quantile_mapping$corrected

  alpha[event_year_indices] <-
    measurement_parameters$alpha

  beta[event_year_indices] <-
    measurement_parameters$beta

  vdelta[event_year_indices] <-
    measurement_parameters$vdelta

  smoother <- kalman_filter_smoother_annual(
    mu = prior_fit$mu,
    M = prior_fit$M,
    r2 = prior_fit$r2,
    Xstar = Xstar,
    alpha = alpha,
    beta = beta,
    vdelta = vdelta
  )

  lme_annual_mean <- rowMeans(X)

  event_posterior <- smoother$X_smooth[
    event_year_indices
  ]

  event_lme <- lme_annual_mean[
    event_year_indices
  ]

  valid_weights <- smoother$observation_weight[
    event_year_indices
  ]

  valid_beta <- beta[
    event_year_indices
  ]

  mapping_limit <- max(
    abs(quantile_mapping$latent_range)
  )

  diagnostics <- data.frame(
    location_id = location_id,
    long = location_long,
    lati = location_lat,
    figure6_cluster =
      selected_locations$figure6_cluster[
        location_index
      ],
    lambda1 = fixed_lambda1,
    lambda2 = fixed_lambda2,
    lambda3 = fixed_lambda3,
    prior_converged = prior_fit$converged,
    prior_iterations = prior_fit$iters,
    prior_objective = prior_fit$obj,
    number_of_event_years = length(event_years),
    number_of_prediction_only_years =
      number_of_years -
      length(event_years),
    minimum_reaches_mspe = min(mspe),
    maximum_reaches_mspe = max(mspe),
    number_of_vhat_values_clipped =
      number_vhat_clipped,
    minimum_beta = min(valid_beta),
    median_beta = median(valid_beta),
    maximum_beta = max(valid_beta),
    number_abs_beta_below_0_01 = sum(
      abs(valid_beta) < 0.01
    ),
    minimum_observation_weight = min(
      valid_weights,
      na.rm = TRUE
    ),
    median_observation_weight = median(
      valid_weights,
      na.rm = TRUE
    ),
    maximum_observation_weight = max(
      valid_weights,
      na.rm = TRUE
    ),
    maximum_filter_identity_error =
      smoother$maximum_filter_identity_error,
    mapping_latent_lower =
      quantile_mapping$latent_range[1],
    mapping_latent_upper =
      quantile_mapping$latent_range[2],
    maximum_absolute_mc_latent_draw =
      measurement_parameters$
        maximum_absolute_latent_draw,
    mc_draw_exceeded_mapping_range =
      measurement_parameters$
        maximum_absolute_latent_draw >
      mapping_limit,
    proportion_posterior_below_both_inputs = mean(
      event_posterior <
        pmin(
          quantile_mapping$corrected,
          event_lme
        )
    ),
    proportion_posterior_above_both_inputs = mean(
      event_posterior >
        pmax(
          quantile_mapping$corrected,
          event_lme
        )
    ),
    mean_posterior_celsius = mean(
      smoother$X_smooth
    ),
    minimum_posterior_variance = min(
      smoother$P_smooth
    ),
    maximum_posterior_variance = max(
      smoother$P_smooth
    )
  )

  list(
    posterior_mean = smoother$X_smooth,
    posterior_variance = smoother$P_smooth,
    diagnostics = diagnostics
  )
}

# ------------------------------------------------------------
# 9. Run all requested locations
# ------------------------------------------------------------

start_time <- Sys.time()

for (location_index in seq_len(number_of_locations)) {

  if (completed_locations[location_index]) {
    next
  }

  message(
    "\nSpatial assimilation ",
    location_index,
    "/",
    number_of_locations,
    ": location_id = ",
    selected_locations$location_id[location_index],
    ", cluster = ",
    selected_locations$figure6_cluster[location_index]
  )

  location_result <- tryCatch(
    process_location(location_index),
    error = function(error_condition) {
      save_checkpoint()

      stop(
        "Spatial assimilation failed at location_id ",
        selected_locations$location_id[location_index],
        ": ",
        conditionMessage(error_condition),
        call. = FALSE
      )
    }
  )

  posterior_mean_matrix[location_index, ] <-
    location_result$posterior_mean

  posterior_variance_matrix[location_index, ] <-
    location_result$posterior_variance

  diagnostic_rows[[location_index]] <-
    location_result$diagnostics

  completed_locations[location_index] <- TRUE

  if (
    location_index %% checkpoint_every == 0L ||
      all(completed_locations)
  ) {
    save_checkpoint()
  }

  message(
    "Completed location_id ",
    selected_locations$location_id[location_index],
    "; elapsed minutes = ",
    round(
      as.numeric(
        difftime(
          Sys.time(),
          start_time,
          units = "mins"
        )
      ),
      2
    )
  )

  invisible(gc())
}

if (!all(completed_locations)) {
  stop(
    "The spatial assimilation ended before every location was completed."
  )
}

if (
  any(!is.finite(posterior_mean_matrix)) ||
    any(!is.finite(posterior_variance_matrix)) ||
    any(posterior_variance_matrix < 0)
) {
  stop(
    "The completed spatial posterior contains invalid values."
  )
}

# ------------------------------------------------------------
# 10. Save final spatial products
# ------------------------------------------------------------

make_grid_output <- function(
    value_matrix) {

  output <- cbind(
    selected_locations %>%
      dplyr::select(
        location_id,
        lati,
        long,
        figure6_cluster
      ),
    as.data.frame(
      value_matrix,
      check.names = FALSE
    )
  )

  names(output)[
    5:ncol(output)
  ] <- paste0(
    "x",
    analysis_years
  )

  output
}

posterior_mean_output <- make_grid_output(
  posterior_mean_matrix
)

posterior_variance_output <- make_grid_output(
  posterior_variance_matrix
)

spatial_diagnostics <- bind_rows(
  diagnostic_rows
) %>%
  arrange(
    figure6_cluster,
    location_id
  )

elapsed_minutes <- as.numeric(
  difftime(
    Sys.time(),
    start_time,
    units = "mins"
  )
)

spatial_metadata <- data.frame(
  location_scope = location_scope,
  number_of_locations = number_of_locations,
  number_of_annual_years = number_of_years,
  first_year = min(analysis_years),
  last_year = max(analysis_years),
  number_of_reaches_event_years =
    length(event_years),
  process_variance = sigmaY2,
  lambda1 = fixed_lambda1,
  lambda2 = fixed_lambda2,
  lambda3 = fixed_lambda3,
  location_specific_cross_validation = FALSE,
  prior_max_iter = prior_max_iter,
  prior_tol = prior_tol,
  measurement_mc_size = measurement_mc_size,
  mapping_grid_size = mapping_grid_size,
  lme_cdf_grid_size = lme_cdf_grid_size,
  elapsed_minutes = elapsed_minutes,
  lme_archive_file = input_files$lme_archive,
  reaches_mean_file = input_files$reaches_mean,
  reaches_variance_file = input_files$reaches_variance,
  reaches_metadata_file = input_files$reaches_metadata,
  cluster_assignment_file = if (
    location_scope == "figure6"
  ) {
    input_files$clusters
  } else {
    NA_character_
  },
  stringsAsFactors = FALSE
)

readr::write_csv(
  posterior_mean_output,
  output_files$mean
)

readr::write_csv(
  posterior_variance_output,
  output_files$variance
)

readr::write_csv(
  spatial_diagnostics,
  output_files$diagnostics
)

readr::write_csv(
  spatial_metadata,
  output_files$metadata
)

if (save_data_valid_copy) {
  readr::write_csv(
    posterior_mean_output,
    validation_files$mean
  )

  readr::write_csv(
    posterior_variance_output,
    validation_files$variance
  )
}

capture.output(
  sessionInfo(),
  file = output_files$session_info
)

if (
  remove_checkpoint_on_success &&
    file.exists(output_files$checkpoint)
) {
  file.remove(
    output_files$checkpoint
  )
}

message(
  "\nSaved spatial posterior mean: ",
  output_files$mean
)

message(
  "Saved spatial posterior variance: ",
  output_files$variance
)

message(
  "Saved spatial diagnostics: ",
  output_files$diagnostics
)

message(
  "Saved spatial metadata: ",
  output_files$metadata
)

if (save_data_valid_copy) {
  message(
    "Saved Data/Valid posterior copies for Figure S5."
  )
}
