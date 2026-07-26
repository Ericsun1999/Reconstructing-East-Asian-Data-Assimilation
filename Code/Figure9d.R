here::i_am("Code/Figure9d.R")

# ============================================================
# Generate:
#   Figure 9(d)  -- Beijing
#   Figure S3(d) -- Shanghai
#   Figure S4(d) -- Hong Kong
#
# Required inputs:
#   Data/reaches_kriging_city3_mean.csv
#   Data/reaches_kriging_city3_sd.csv
#
#   Data/par/mtB.csv, muB.csv, rtB.csv
#   Data/par/mtS.csv, muS.csv, rtS.csv
#   Data/par/mtH.csv, muH.csv, rtH.csv
#
#   Data/LME data/Figure9/d1.csv  (Hong Kong)
#   Data/LME data/Figure9/d2.csv  (Shanghai)
#   Data/LME data/Figure9/d3.csv  (Beijing)
#
# Outputs:
#   Output/Figure9/Figure9d.png
#   Output/Supplementary/FigureS3d.png
#   Output/Supplementary/FigureS4d.png
#
#   Data/Valid/tempBv5.csv
#   Data/Valid/tempSv5.csv
#   Data/Valid/tempHv5.csv
# ============================================================

library(here)
library(ggplot2)
library(mvtnorm)
library(np)
library(readr)

# ------------------------------------------------------------
# 1. Paths and analysis settings
# ------------------------------------------------------------

kriging_mean_file <- here::here(
  "Data",
  "reaches_kriging_city3_mean.csv"
)

kriging_sd_file <- here::here(
  "Data",
  "reaches_kriging_city3_sd.csv"
)

parameter_dir <- here::here(
  "Data",
  "par"
)

validation_dir <- here::here(
  "Data",
  "Valid"
)

figure9_dir <- here::here(
  "Output",
  "Figure9"
)

supplementary_dir <- here::here(
  "Output",
  "Supplementary"
)

dir.create(
  validation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure9_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  supplementary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

analysis_years <- 1368:1911
lme_columns <- 21:564

# Marginal variance of the latent REACHES process used in the
# original Figure 9(d) analysis.
sigmaY2 <- 0.887

measurement_mc_size <- 10000
measurement_mc_seed <- 1

# Preserve the original behavior of saving and plotting all
# available kriging years except the final one.
drop_final_output_year <- TRUE

candidate_lme_dirs <- c(
  here::here(
    "Data",
    "LME data",
    "Figure9"
  ),
  here::here(
    "Data",
    "LME data"
  )
)

required_lme_files <- c(
  "d1.csv",
  "d2.csv",
  "d3.csv"
)

valid_lme_dir <- vapply(
  candidate_lme_dirs,
  function(directory) {
    all(
      file.exists(
        file.path(
          directory,
          required_lme_files
        )
      )
    )
  },
  logical(1)
)

if (!any(valid_lme_dir)) {
  stop(
    "Could not find d1.csv, d2.csv, and d3.csv together in:\n",
    paste(
      paste0(
        "  - ",
        candidate_lme_dirs
      ),
      collapse = "\n"
    )
  )
}

lme_dir <- candidate_lme_dirs[
  which(valid_lme_dir)[1]
]

city_config <- list(
  HongKong = list(
    code = "H",
    lme_file = file.path(
      lme_dir,
      "d1.csv"
    ),
    long = 113.75,
    lat = 22.25,
    figure_file = file.path(
      supplementary_dir,
      "FigureS4d.png"
    ),
    valid_file = file.path(
      validation_dir,
      "tempHv5.csv"
    )
  ),
  Shanghai = list(
    code = "S",
    lme_file = file.path(
      lme_dir,
      "d2.csv"
    ),
    long = 121.25,
    lat = 31.25,
    figure_file = file.path(
      supplementary_dir,
      "FigureS3d.png"
    ),
    valid_file = file.path(
      validation_dir,
      "tempSv5.csv"
    )
  ),
  Beijing = list(
    code = "B",
    lme_file = file.path(
      lme_dir,
      "d3.csv"
    ),
    long = 116.25,
    lat = 39.75,
    figure_file = file.path(
      figure9_dir,
      "Figure9d.png"
    ),
    valid_file = file.path(
      validation_dir,
      "tempBv5.csv"
    )
  )
)

# ------------------------------------------------------------
# 2. Input helpers
# ------------------------------------------------------------

extract_year_columns <- function(data) {

  column_names <- names(data)

  year_values <- suppressWarnings(
    as.integer(
      sub(
        "^X",
        "",
        column_names
      )
    )
  )

  valid <- !is.na(year_values) &
    year_values >= 1000L &
    year_values <= 3000L

  if (!any(valid)) {
    stop(
      "No year columns were found. Expected names such as ",
      "'1368' or 'X1368'."
    )
  }

  data.frame(
    column = column_names[valid],
    year = year_values[valid],
    stringsAsFactors = FALSE
  )
}


coordinate_key <- function(
    long,
    lat) {

  sprintf(
    "%.2f_%.2f",
    long,
    lat
  )
}


read_kriging_files <- function(
    mean_file,
    sd_file) {

  if (!file.exists(mean_file)) {
    stop(
      "Kriging mean file was not found: ",
      mean_file
    )
  }

  if (!file.exists(sd_file)) {
    stop(
      "Kriging SD file was not found: ",
      sd_file
    )
  }

  mean_data <- read.csv(
    mean_file,
    check.names = FALSE
  )

  sd_data <- read.csv(
    sd_file,
    check.names = FALSE
  )

  required_coordinates <- c(
    "long",
    "lat"
  )

  if (
    !all(required_coordinates %in% names(mean_data)) ||
      !all(required_coordinates %in% names(sd_data))
  ) {
    stop(
      "Both kriging files must contain 'long' and 'lat' columns."
    )
  }

  mean_year_map <- extract_year_columns(
    mean_data
  )

  sd_year_map <- extract_year_columns(
    sd_data
  )

  common_years <- intersect(
    mean_year_map$year,
    sd_year_map$year
  )

  common_years <- sort(
    common_years
  )

  if (length(common_years) < 2L) {
    stop(
      "The kriging mean and SD files do not share at least two years."
    )
  }

  list(
    mean = mean_data,
    sd = sd_data,
    mean_year_map = mean_year_map,
    sd_year_map = sd_year_map,
    years = common_years
  )
}


extract_city_kriging <- function(
    kriging_data,
    city_name,
    long,
    lat) {

  mean_keys <- coordinate_key(
    kriging_data$mean$long,
    kriging_data$mean$lat
  )

  sd_keys <- coordinate_key(
    kriging_data$sd$long,
    kriging_data$sd$lat
  )

  target_key <- coordinate_key(
    long,
    lat
  )

  mean_row <- match(
    target_key,
    mean_keys
  )

  sd_row <- match(
    target_key,
    sd_keys
  )

  if (is.na(mean_row) || is.na(sd_row)) {
    stop(
      "The kriging location for ",
      city_name,
      " was not found at (",
      long,
      ", ",
      lat,
      ")."
    )
  }

  years <- kriging_data$years

  mean_columns <- kriging_data$mean_year_map$column[
    match(
      years,
      kriging_data$mean_year_map$year
    )
  ]

  sd_columns <- kriging_data$sd_year_map$column[
    match(
      years,
      kriging_data$sd_year_map$year
    )
  ]

  yhat <- as.numeric(
    kriging_data$mean[
      mean_row,
      mean_columns,
      drop = TRUE
    ]
  )

  nu <- as.numeric(
    kriging_data$sd[
      sd_row,
      sd_columns,
      drop = TRUE
    ]
  )

  if (
    any(!is.finite(yhat)) ||
      any(!is.finite(nu))
  ) {
    stop(
      "Non-finite kriging values were found for ",
      city_name,
      "."
    )
  }

  if (any(nu < 0)) {
    stop(
      "Negative kriging standard deviations were found for ",
      city_name,
      "."
    )
  }

  list(
    years = years,
    yhat = yhat,
    nu = nu
  )
}


read_lme_city <- function(
    input_file,
    city_name) {

  if (!file.exists(input_file)) {
    stop(
      "LME file for ",
      city_name,
      " was not found: ",
      input_file
    )
  }

  raw_data <- read.csv(
    input_file,
    row.names = 1,
    check.names = FALSE
  )

  if (max(lme_columns) > ncol(raw_data)) {
    stop(
      city_name,
      " LME file contains only ",
      ncol(raw_data),
      " data columns after removing row names; ",
      "columns 21:564 are required."
    )
  }

  selected_data <- raw_data[
    ,
    lme_columns,
    drop = FALSE
  ]

  X <- t(
    data.matrix(
      selected_data
    )
  ) -
    273.15

  storage.mode(X) <- "double"

  if (
    nrow(X) != length(analysis_years) ||
      ncol(X) != 13L
  ) {
    stop(
      city_name,
      " LME data should produce a ",
      length(analysis_years),
      " x 13 matrix, but produced ",
      nrow(X),
      " x ",
      ncol(X),
      "."
    )
  }

  if (any(!is.finite(X))) {
    stop(
      "Non-finite LME temperatures were found for ",
      city_name,
      "."
    )
  }

  rownames(X) <- as.character(
    analysis_years
  )

  list(
    matrix = X,
    sample = as.numeric(X),
    annual_mean = setNames(
      rowMeans(X),
      analysis_years
    )
  )
}


read_penalized_parameter <- function(
    filename,
    parameter_name,
    years) {

  if (!file.exists(filename)) {
    stop(
      "Parameter file was not found: ",
      filename
    )
  }

  data <- read.csv(
    filename,
    check.names = FALSE
  )

  required_columns <- c(
    "year",
    "coefficient",
    "value"
  )

  if (!all(required_columns %in% names(data))) {
    stop(
      filename,
      " must contain: ",
      paste(
        required_columns,
        collapse = ", "
      ),
      "."
    )
  }

  penalized <- data[
    data$coefficient == "Penalized ML",
    required_columns,
    drop = FALSE
  ]

  if (nrow(penalized) == 0L) {
    stop(
      "No 'Penalized ML' rows were found in ",
      filename,
      "."
    )
  }

  values <- penalized$value[
    match(
      years,
      penalized$year
    )
  ]

  values <- as.numeric(
    values
  )

  if (
    length(values) != length(years) ||
      any(!is.finite(values))
  ) {
    stop(
      "Could not align ",
      parameter_name,
      " with all requested years in ",
      filename,
      "."
    )
  }

  values
}

# ------------------------------------------------------------
# 3. Quantile-mapping helpers
# ------------------------------------------------------------

make_quantile_mapping <- function(
    lme_sample,
    yhat,
    nu,
    bandwidth_seed = 10L) {

  set.seed(
    bandwidth_seed
  )

  fx_bw <- npudistbw(
    dat = lme_sample
  )

  Fx_hat <- function(q) {

    fitted(
      npudist(
        bws = fx_bw,
        edat = data.frame(
          x = q
        )
      )
    )
  }

  FY_hat <- function(
      y,
      yhat,
      nu) {

    yhat <- as.numeric(
      yhat
    )

    nu <- pmax(
      as.numeric(nu),
      1e-8
    )

    stopifnot(
      length(yhat) ==
        length(nu)
    )

    vapply(
      y,
      function(yy) {
        mean(
          pnorm(
            (
              yy -
                yhat
            ) /
              nu
          )
        )
      },
      numeric(1)
    )
  }

  qx_min <- min(
    lme_sample
  ) -
    5 *
    sd(
      lme_sample
    )

  qx_max <- max(
    lme_sample
  ) +
    5 *
    sd(
      lme_sample
    )

  Fx_inv <- function(u) {

    u <- pmin(
      pmax(
        u,
        1e-8
      ),
      1 -
        1e-8
    )

    vapply(
      u,
      function(uu) {

        uniroot(
          function(q) {
            Fx_hat(q) -
              uu
          },
          interval = c(
            qx_min,
            qx_max
          ),
          tol = 1e-6
        )$root
      },
      numeric(1)
    )
  }

  ycorrected <- Fx_inv(
    FY_hat(
      yhat,
      yhat,
      nu
    )
  )

  bandwidth <- as.numeric(
    fx_bw$bw
  )

  if (
    length(bandwidth) != 1L ||
      !is.finite(bandwidth) ||
      bandwidth <= 0
  ) {
    stop(
      "The estimated univariate np bandwidth is invalid."
    )
  }

  g_of <- function(y) {

    u <- FY_hat(
      y,
      yhat = yhat,
      nu = nu
    )

    lower <- min(
      lme_sample
    ) -
      5 *
      bandwidth

    upper <- max(
      lme_sample
    ) +
      5 *
      bandwidth

    uniroot(
      function(q) {
        Fx_hat(q) -
          u
      },
      lower = lower,
      upper = upper,
      tol = 1e-8
    )$root
  }

  y_grid <- seq(
    -2.5,
    2.5,
    length.out = 1000
  )

  g_grid <- vapply(
    y_grid,
    g_of,
    numeric(1)
  )

  g_fast <- approxfun(
    y_grid,
    g_grid,
    rule = 2
  )

  list(
    corrected = ycorrected,
    g = g_fast,
    bandwidth = bandwidth
  )
}

# ------------------------------------------------------------
# 4. Measurement-model Monte Carlo
# ------------------------------------------------------------

round_and_clip <- function(
    x,
    lower = -2,
    upper = 1,
    digits = 2) {

  x <- round(
    x,
    digits = digits
  )

  pmin(
    pmax(
      x,
      lower
    ),
    upper
  )
}


compute_meas_params_mc <- function(
    sigmaY2,
    v,
    g,
    n_mc = 50000,
    seed = 1) {

  set.seed(seed)

  v <- as.numeric(v)

  if (
    any(!is.finite(v)) ||
      any(v < 0) ||
      any(v > sigmaY2)
  ) {
    stop(
      "All values of v must lie between 0 and sigmaY2."
    )
  }

  N <- length(v)

  alpha <- numeric(N)
  beta <- numeric(N)
  vdelta <- numeric(N)

  for (t in seq_len(N)) {

    vt <- v[t]

    Sigma <- matrix(
      c(
        sigmaY2,
        vt,
        vt,
        vt
      ),
      nrow = 2,
      byrow = TRUE
    )

    eigenvalues <- eigen(
      Sigma,
      symmetric = TRUE,
      only.values = TRUE
    )$values

    if (min(eigenvalues) <= 0) {
      Sigma <- Sigma +
        diag(
          abs(
            min(eigenvalues)
          ) +
            1e-10,
          2
        )
    }

    sample_values <- rmvnorm(
      n = n_mc,
      mean = c(
        0,
        0
      ),
      sigma = Sigma
    )

    Y <- round_and_clip(
      sample_values[, 1],
      lower = -2,
      upper = 1,
      digits = 3
    )

    Yhat <- round_and_clip(
      sample_values[, 2],
      lower = -2,
      upper = 1,
      digits = 3
    )

    X <- g(Y)
    Xstar <- g(Yhat)

    EX <- mean(X)
    EXstar <- mean(Xstar)
    VX <- var(X)
    covariance <- cov(
      Xstar,
      X
    )

    if (
      !is.finite(VX) ||
        VX <= 0
    ) {
      stop(
        "A non-positive Monte Carlo variance was obtained at index ",
        t,
        "."
      )
    }

    beta_t <- covariance /
      VX

    alpha_t <- EXstar -
      beta_t *
      EX

    delta <- Xstar -
      alpha_t -
      beta_t *
      X

    alpha[t] <- alpha_t
    beta[t] <- beta_t
    vdelta[t] <- max(
      var(delta),
      1e-10
    )
  }

  list(
    alpha = alpha,
    beta = beta,
    vdelta = vdelta
  )
}

# ------------------------------------------------------------
# 5. Kalman filter and RTS smoother
# ------------------------------------------------------------

kalman_filter_smoother <- function(
    mu,
    M,
    r2,
    Xstar,
    alpha,
    beta,
    vdelta) {

  N <- length(mu)

  stopifnot(
    length(r2) == N,
    length(Xstar) == N,
    length(alpha) == N,
    length(beta) == N,
    length(vdelta) == N,
    length(M) == N - 1L
  )

  X_pred <- numeric(N)
  P_pred <- numeric(N)
  X_filt <- numeric(N)
  P_filt <- numeric(N)

  X_pred[1] <- mu[1]
  P_pred[1] <- r2[1]

  denominator_1 <- beta[1]^2 *
    P_pred[1] +
    vdelta[1]

  K1 <- beta[1] *
    P_pred[1] /
    denominator_1

  X_filt[1] <- X_pred[1] +
    K1 *
    (
      Xstar[1] -
        alpha[1] -
        beta[1] *
        X_pred[1]
    )

  P_filt[1] <- (
    1 -
      beta[1] *
      K1
  ) *
    P_pred[1]

  for (t in 2:N) {

    transition <- M[
      t - 1
    ]

    X_pred[t] <- mu[t] +
      transition *
      (
        X_filt[t - 1] -
          mu[t - 1]
      )

    P_pred[t] <- transition^2 *
      P_filt[t - 1] +
      r2[t]

    denominator_t <- beta[t]^2 *
      P_pred[t] +
      vdelta[t]

    gain <- beta[t] *
      P_pred[t] /
      denominator_t

    X_filt[t] <- X_pred[t] +
      gain *
      (
        Xstar[t] -
          alpha[t] -
          beta[t] *
          X_pred[t]
      )

    P_filt[t] <- (
      1 -
        beta[t] *
        gain
    ) *
      P_pred[t]
  }

  X_smooth <- numeric(N)
  P_smooth <- numeric(N)
  J <- numeric(N - 1)

  X_smooth[N] <- X_filt[N]
  P_smooth[N] <- P_filt[N]

  if (N >= 2L) {
    for (t in seq(
      from = N - 1L,
      to = 1L,
      by = -1L
    )) {

      transition <- M[t]

      J[t] <- P_filt[t] *
        transition /
        P_pred[t + 1]

      X_smooth[t] <- X_filt[t] +
        J[t] *
        (
          X_smooth[t + 1] -
            X_pred[t + 1]
        )

      P_smooth[t] <- P_filt[t] +
        J[t]^2 *
        (
          P_smooth[t + 1] -
            P_pred[t + 1]
        )
    }
  }

  P_smooth <- pmax(
    P_smooth,
    0
  )

  list(
    X_pred = X_pred,
    P_pred = P_pred,
    X_filt = X_filt,
    P_filt = P_filt,
    X_smooth = X_smooth,
    P_smooth = P_smooth,
    J = J
  )
}

# ------------------------------------------------------------
# 6. Process one city
# ------------------------------------------------------------

process_city <- function(
    city_name,
    config,
    kriging_data,
    city_index) {

  message(
    "\nProcessing ",
    city_name,
    "..."
  )

  city_kriging <- extract_city_kriging(
    kriging_data = kriging_data,
    city_name = city_name,
    long = config$long,
    lat = config$lat
  )

  # Restrict to years represented in the LME prior.
  keep <- city_kriging$years %in%
    analysis_years

  years <- city_kriging$years[
    keep
  ]

  yhat <- city_kriging$yhat[
    keep
  ]

  nu <- city_kriging$nu[
    keep
  ]

  if (length(years) < 2L) {
    stop(
      "Fewer than two common years were found for ",
      city_name,
      "."
    )
  }

  lme <- read_lme_city(
    input_file = config$lme_file,
    city_name = city_name
  )

  quantile_mapping <- make_quantile_mapping(
    lme_sample = lme$sample,
    yhat = yhat,
    nu = nu,
    bandwidth_seed = 10L +
      city_index
  )

  ycorrected <- quantile_mapping$corrected

  mt_file <- file.path(
    parameter_dir,
    paste0(
      "mt",
      config$code,
      ".csv"
    )
  )

  mu_file <- file.path(
    parameter_dir,
    paste0(
      "mu",
      config$code,
      ".csv"
    )
  )

  rt_file <- file.path(
    parameter_dir,
    paste0(
      "rt",
      config$code,
      ".csv"
    )
  )

  mu_values <- read_penalized_parameter(
    filename = mu_file,
    parameter_name = "mu",
    years = years
  )

  r2_values <- read_penalized_parameter(
    filename = rt_file,
    parameter_name = "r2",
    years = years
  )

  transition_years <- years[
    -length(years)
  ]

  M_values <- read_penalized_parameter(
    filename = mt_file,
    parameter_name = "M",
    years = transition_years
  )

  # The kriging uncertainty input is a standard deviation.
  # Therefore, prediction error variance is nu^2, and
  # Var(Yhat_t) = Var(Y_t) - Var(Y_t - Yhat_t).
  prediction_error_variance <- nu^2

  v <- sigmaY2 -
    prediction_error_variance

  tolerance <- 1e-10

  if (
    any(
      v < -tolerance |
        v > sigmaY2 +
        tolerance
    )
  ) {
    warning(
      city_name,
      ": some implied Var(Yhat_t) values were outside ",
      "[0, sigmaY2] and were clipped to that interval."
    )
  }

  v <- pmin(
    pmax(
      v,
      0
    ),
    sigmaY2
  )

  measurement_parameters <- compute_meas_params_mc(
    sigmaY2 = sigmaY2,
    v = v,
    g = quantile_mapping$g,
    n_mc = measurement_mc_size,
    seed = measurement_mc_seed
  )

  smoother_output <- kalman_filter_smoother(
    mu = mu_values,
    M = M_values,
    r2 = r2_values,
    Xstar = ycorrected,
    alpha = measurement_parameters$alpha,
    beta = measurement_parameters$beta,
    vdelta = measurement_parameters$vdelta
  )

  lme_mean <- as.numeric(
    lme$annual_mean[
      as.character(years)
    ]
  )

  result <- data.frame(
    year = years,
    predicted = smoother_output$X_smooth,
    sigmasmooth2 = smoother_output$P_smooth,
    REACHES = ycorrected,
    alpha = measurement_parameters$alpha,
    beta = measurement_parameters$beta,
    vdelata = measurement_parameters$vdelta,
    LME = lme_mean,
    check.names = FALSE
  )

  if (drop_final_output_year) {
    result_for_output <- result[
      -nrow(result),
      ,
      drop = FALSE
    ]
  } else {
    result_for_output <- result
  }

  figure <- ggplot(
    data = result_for_output,
    aes(
      x = year,
      y = predicted
    )
  ) +
    geom_line(
      color = "red"
    ) +
    geom_ribbon(
      aes(
        ymin = predicted -
          sqrt(sigmasmooth2),
        ymax = predicted +
          sqrt(sigmasmooth2)
      ),
      alpha = 0.5,
      fill = "grey3"
    ) +
    xlab("year") +
    ylab("temperature") +
    theme(
      text = element_text(
        size = 11
      ),
      legend.position = "right",
      legend.key.height = grid::unit(
        1.5,
        "cm"
      ),
      plot.title = element_text(
        hjust = 0.5
      )
    )

  ggsave(
    filename = config$figure_file,
    plot = figure,
    width = 6,
    height = 3,
    units = "in",
    dpi = 300
  )

  readr::write_excel_csv(
    result_for_output,
    config$valid_file
  )

  message(
    "Saved figure: ",
    config$figure_file
  )

  message(
    "Saved validation data: ",
    config$valid_file
  )

  list(
    full_result = result,
    output_result = result_for_output,
    figure = figure,
    bandwidth = quantile_mapping$bandwidth
  )
}

# ------------------------------------------------------------
# 7. Run all three cities
# ------------------------------------------------------------

kriging_data <- read_kriging_files(
  mean_file = kriging_mean_file,
  sd_file = kriging_sd_file
)

figure9d_results <- vector(
  "list",
  length(city_config)
)

names(figure9d_results) <- names(
  city_config
)

city_names <- names(
  city_config
)

for (i in seq_along(city_names)) {

  city_name <- city_names[i]

  figure9d_results[[city_name]] <- process_city(
    city_name = city_name,
    config = city_config[[city_name]],
    kriging_data = kriging_data,
    city_index = i
  )
}

capture.output(
  sessionInfo(),
  file = file.path(
    validation_dir,
    "Figure9d_sessionInfo.txt"
  )
)

message(
  "\nCompleted Figure 9(d), Figure S3(d), and Figure S4(d)."
)
