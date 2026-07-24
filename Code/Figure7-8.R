here::i_am("Code/Figure7-8.R")

# ============================================================
# Quantile-mapping calibration for Figures 7 and 8
#
# This script produces:
#
# Figure 7:
#   Figure7a.png
#   Figure7b.png
#   Figure7c.png
#   Figure7d.png
#
# Figure 8:
#   Figure8a.png
#   Figure8b.png
#   Figure8c.png
#   Figure8d.png
#
# Figure 7(e), the spatial temperature map, is generated
# separately by Code/Figure7e.R
#
# City ordering in the REACHES files:
#   row 1 = Hong Kong
#   row 2 = Shanghai
#   row 3 = Beijing
# ============================================================

library(ggplot2)
library(np)

# ------------------------------------------------------------
# 1. Read shared REACHES inputs
# ------------------------------------------------------------

tempe_all_data <- read.csv(
  here::here("Data", "tempe_all_v3.csv"),
  header = FALSE,
  check.names = FALSE
)

# The first row contains the available REACHES years.
year3 <- suppressWarnings(
  as.integer(
    unlist(
      tempe_all_data[1, -c(1, 2), drop = FALSE ],
      use.names = FALSE
    )
  )
)

# Rows 2--4 correspond to Hong Kong, Shanghai, and Beijing.
tempe_all <- tempe_all_data[2:4, , drop = FALSE]

nu_all <- read.csv(
  here::here("Data", "tempe_all_std.csv"),
  check.names = FALSE
)


# ------------------------------------------------------------
# 2. City-specific configuration
# ------------------------------------------------------------

city_config <- list(
  Beijing = list(
    location_row = 3L,
    lme_file = here::here(
      "Data",
      "LME data",
      "d3.csv"
    ),
    figure8_panel = "b"
  ),
  Shanghai = list(
    location_row = 2L,
    lme_file = here::here(
      "Data",
      "LME data",
      "d2.csv"
    ),
    figure8_panel = "c"
  ),
  HongKong = list(
    location_row = 1L,
    lme_file = here::here(
      "Data",
      "LME data",
      "d1.csv"
    ),
    figure8_panel = "d"
  )
)

output_dir <- here::here(
  "Output",
  "Figure7-8"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 3. Quantile-mapping functions
# ------------------------------------------------------------

# Estimate the LME temperature CDF.
Fx_hat <- function(q, fx_bw) {
  fx_ob <- npudist(
    bws = fx_bw,
    edat = data.frame(x = q)
  )

  as.numeric(fitted(fx_ob))
}

# Estimate the uncertainty-aware REACHES CDF.
FY_hat <- function(y, yhat, nu) {
  yhat <- as.numeric(yhat)
  nu <- as.numeric(nu)

  if (length(yhat) != length(nu)) {
    stop("yhat and nu must have the same length.")
  }

  if (any(!is.finite(yhat))) {
    stop("yhat contains non-finite values.")
  }

  if (any(!is.finite(nu))) {
    stop("nu contains non-finite values.")
  }

  # Avoid division by zero.
  nu_safe <- pmax(nu, 1e-8)

  vapply(
    y,
    function(yy) {
      mean(
        pnorm((yy - yhat) / nu_safe)
      )
    },
    numeric(1)
  )
}

# Numerically invert the estimated LME CDF.
Fx_inv <- function(
    u,
    qx_min,
    qx_max,
    fx_bw) {

  u <- pmin(
    pmax(u, 1e-8),
    1 - 1e-8
  )

  vapply(
    u,
    function(uu) {
      objective <- function(q) {
        Fx_hat(q, fx_bw) - uu
      }

      uniroot(
        objective,
        interval = c(qx_min, qx_max),
        tol = 1e-6
      )$root
    },
    numeric(1)
  )
}

qmapping <- function(
    yhat,
    std,
    x,
    xhat,
    ymax = 1,
    ymin = -2) {

  yhat <- as.numeric(yhat)
  std <- as.numeric(std)
  x <- as.numeric(x)
  xhat <- as.numeric(xhat)

  if (any(!is.finite(x))) {
    stop("The LME calibration sample contains non-finite values.")
  }

  # Estimate F_X using npudist, with bandwidth selected
  # by least-squares cross-validation.
  fx_bw <- npudistbw(
    dat = data.frame(x = x)
  )

  qx_min <- min(x) - 5 * sd(x)
  qx_max <- max(x) + 5 * sd(x)

  y_seq <- seq(
    ymin,
    ymax,
    length.out = 150
  )

  FY_values <- FY_hat(
    y = y_seq,
    yhat = yhat,
    nu = std
  )

  FX_values <- Fx_inv(
    u = FY_values,
    qx_min = qx_min,
    qx_max = qx_max,
    fx_bw = fx_bw
  )

  ycorrected <- Fx_inv(
    u = FY_hat(
      y = xhat,
      yhat = yhat,
      nu = std
    ),
    qx_min = qx_min,
    qx_max = qx_max,
    fx_bw = fx_bw
  )

  list(
    y_seq = y_seq,
    FX_values = FX_values,
    ycorrected = ycorrected
  )
}

# ------------------------------------------------------------
# 4. Prepare quantile-mapping results for one city
# ------------------------------------------------------------

prepare_city_result <- function(
    city_name,
    config) {

  location_row <- config$location_row
  lme_file <- config$lme_file

  # ----------------------------------------------------------
  # Read city-specific LME data
  # ----------------------------------------------------------

  lme_data <- read.csv(
    lme_file,
    row.names = 1,
    check.names = FALSE
  )


  lme_temperature_data <- lme_data[
    ,
    -c(1, 2),
    drop = FALSE
  ]

  lme_temperature_data[] <- lapply(
    lme_temperature_data,
    function(x) {
      as.numeric(as.character(x))
    }
  )

  lme_temperature <- as.matrix(
    lme_temperature_data
  ) - 273.15

  # Mean across LME ensemble members for each year.
  lme_mean <- colMeans(lme_temperature)

  # Pooled LME sample used to estimate F_X.
  lme_sample <- as.numeric(
    t(lme_temperature)
  )

  # ----------------------------------------------------------
  # Extract city-specific REACHES values and uncertainties
  # ----------------------------------------------------------

  yhat <- suppressWarnings(
    as.numeric(
      unlist(
        tempe_all[
          location_row,
          -c(1, 2),
          drop = FALSE
        ],
        use.names = FALSE
      )
    )
  )

  nu <- suppressWarnings(
    as.numeric(
      unlist(
        nu_all[
          location_row,
          -c(1, 2),
          drop = FALSE
        ],
        use.names = FALSE
      )
    )
  )

  if (
    length(year3) != length(yhat) ||
      length(year3) != length(nu)
  ) {
    stop(
      "The year, REACHES prediction, and uncertainty vectors ",
      "have different lengths for ",
      city_name,
      "."
    )
  }

  # ----------------------------------------------------------
  # Quantile mapping
  # ----------------------------------------------------------

  qm <- qmapping(
    yhat = yhat,
    std = nu,
    x = lme_sample,
    xhat = yhat,
    ymax = 1,
    ymin = -2
  )

  # ----------------------------------------------------------
  # Match LME annual means to the available REACHES years
  # ----------------------------------------------------------

  lme_years <- 1350 + seq_along(lme_mean) - 1L

  lme_year_index <- match(
    year3,
    lme_years
  )

  if (anyNA(lme_year_index)) {
    stop(
      "One or more REACHES years for ",
      city_name,
      " fall outside the available LME period."
    )
  }

  list(
    city_name = city_name,
    years = year3,
    yhat = yhat,
    nu = nu,
    ycorrected = qm$ycorrected,
    transformation_index = qm$y_seq,
    transformation_temperature = qm$FX_values,
    lme_sample = lme_sample,
    lme_mean = lme_mean[lme_year_index]
  )
}

# ------------------------------------------------------------
# 5. Prepare results for all three cities
# ------------------------------------------------------------

city_results <- lapply(
  names(city_config),
  function(city_name) {
    message(
      "Preparing quantile-mapping result for ",
      city_name,
      "..."
    )

    prepare_city_result(
      city_name = city_name,
      config = city_config[[city_name]]
    )
  }
)

names(city_results) <- names(city_config)

beijing_result <- city_results[["Beijing"]]

# ------------------------------------------------------------
# 6. Figure 7(a): Beijing REACHES index time series
# ------------------------------------------------------------

df_figure7a <- data.frame(
  year = beijing_result$years,
  reach = beijing_result$yhat
)

p_figure7a <- ggplot(
  data = df_figure7a,
  aes(x = year, y = reach)
) +
  geom_line(
    colour = "darkorchid1"
  ) +
  labs(
    x = "year",
    y = "level"
  ) +
  theme(
    text = element_text(size = 19),
    legend.position = "right",
    legend.key.height = grid::unit(1.5, "cm"),
    plot.title = element_text(hjust = 0.5)
  )

# ------------------------------------------------------------
# 7. Figure 7(b): Beijing calibration function
# ------------------------------------------------------------

df_figure7b <- data.frame(
  index = beijing_result$transformation_index,
  temperature = beijing_result$transformation_temperature
)

p_figure7b <- ggplot(
  data = df_figure7b,
  aes(x = index, y = temperature)
) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE
  ) +
  labs(
    x = "index",
    y = "temperature"
  ) +
  theme(
    text = element_text(size = 19),
    legend.position = "right",
    legend.key.height = grid::unit(1.5, "cm"),
    plot.title = element_text(hjust = 0.5)
  )

# ------------------------------------------------------------
# 8. Figure 7(c): Beijing REACHES index distribution
# ------------------------------------------------------------

p_figure7c <- ggplot(
  data = df_figure7a,
  aes(x = reach)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = seq(
      -2.5,
      1.5,
      by = 0.2
    ),
    fill = "darkorchid1",
    colour = "white"
  ) +
  geom_density(
    colour = "black",
    linewidth = 0.4
  ) +
  labs(
    x = "temperature",
    y = "density"
  ) +
  theme(
    text = element_text(size = 18),
    legend.position = "right",
    legend.key.height = grid::unit(1.5, "cm"),
    plot.title = element_text(hjust = 0.5)
  )

# ------------------------------------------------------------
# 9. Figure 7(d): Beijing LME temperature distribution
# ------------------------------------------------------------

df_figure7d <- data.frame(
  temper = beijing_result$lme_sample
)

p_figure7d <- ggplot(
  data = df_figure7d,
  aes(x = temper)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = seq(
      9,
      14,
      by = 0.1
    ),
    colour = "white",
    fill = "deepskyblue"
  ) +
  geom_density(
    colour = "black",
    linewidth = 0.4
  ) +
  scale_x_continuous(
    limits = c(9, 14)
  ) +
  labs(
    x = "temperature",
    y = "density"
  ) +
  theme(
    text = element_text(size = 18),
    legend.position = "right",
    legend.key.height = grid::unit(1.5, "cm"),
    plot.title = element_text(hjust = 0.5)
  )

# ------------------------------------------------------------
# 10. Figure 8(a): Beijing reconstruction with uncertainty
# ------------------------------------------------------------

df_figure8a <- data.frame(
  year = beijing_result$years,
  REACHES = beijing_result$ycorrected,
  std = beijing_result$nu
)

p_figure8a <- ggplot(
  data = df_figure8a,
  aes(x = year, y = REACHES)
) +
  geom_ribbon(
    aes(
      ymin = REACHES - std,
      ymax = REACHES + std
    ),
    alpha = 0.5,
    fill = "grey3"
  ) +
  geom_line(
    colour = "red"
  ) +
  labs(
    x = "year",
    y = "temperature"
  ) +
  theme(
    text = element_text(size = 11),
    legend.position = "right",
    legend.key.height = grid::unit(1.5, "cm"),
    plot.title = element_text(hjust = 0.5)
  )

# ------------------------------------------------------------
# 11. Function for Figure 8(b)--(d)
# ------------------------------------------------------------

make_figure8_city_plot <- function(city_result) {
  corrected_data <- data.frame(
    year = city_result$years,
    temperature = city_result$ycorrected
  )

  lme_data <- data.frame(
    year = city_result$years,
    temperature = city_result$lme_mean
  )

  ggplot() +
    geom_line(
      data = corrected_data,
      aes(
        x = year,
        y = temperature
      ),
      colour = "palevioletred1",
      alpha = 0.75
    ) +
    geom_smooth(
      data = corrected_data,
      aes(
        x = year,
        y = temperature
      ),
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 1.1,
      colour = "firebrick",
      linetype = "dashed",
      span = 0.75
    ) +
    geom_line(
      data = lme_data,
      aes(
        x = year,
        y = temperature
      ),
      colour = "skyblue1",
      alpha = 0.75
    ) +
    geom_smooth(
      data = lme_data,
      aes(
        x = year,
        y = temperature
      ),
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 1.1,
      colour = "deepskyblue",
      linetype = "dashed"
    ) +
    labs(
      x = "year",
      y = "temperature"
    ) +
    theme(
      text = element_text(size = 12),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
}

p_figure8b <- make_figure8_city_plot(
  city_results[["Beijing"]]
)

p_figure8c <- make_figure8_city_plot(
  city_results[["Shanghai"]]
)

p_figure8d <- make_figure8_city_plot(
  city_results[["HongKong"]]
)

# ------------------------------------------------------------
# 12. Save all Figure 7--8 panels
# ------------------------------------------------------------

save_panel <- function(
    plot_object,
    filename,
    width,
    height) {

  output_file <- file.path(
    output_dir,
    filename
  )

  ggsave(
    filename = output_file,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )

  message("Saved: ", output_file)

  invisible(output_file)
}

output_files <- c(
  save_panel(
    p_figure7a,
    "Figure7a.png",
    width = 6,
    height = 4
  ),
  save_panel(
    p_figure7b,
    "Figure7b.png",
    width = 6,
    height = 4
  ),
  save_panel(
    p_figure7c,
    "Figure7c.png",
    width = 6,
    height = 3.5
  ),
  save_panel(
    p_figure7d,
    "Figure7d.png",
    width = 6,
    height = 3.5
  ),
  save_panel(
    p_figure8a,
    "Figure8a.png",
    width = 6,
    height = 3
  ),
  save_panel(
    p_figure8b,
    "Figure8b.png",
    width = 6,
    height = 3
  ),
  save_panel(
    p_figure8c,
    "Figure8c.png",
    width = 6,
    height = 3
  ),
  save_panel(
    p_figure8d,
    "Figure8d.png",
    width = 6,
    height = 3
  )
)
