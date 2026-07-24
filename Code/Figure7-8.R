here::i_am("Code/Figure7-8.R")

# ============================================================
# Quantile mapping and city-specific outputs for Figures 7--8
#
# Cities:
#   1. Hong Kong
#   2. Shanghai
#   3. Beijing
#
# All figures are generated automatically and saved under:
#   Output/Figure7-8/
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

# The first row contains the years.
year3 <- as.integer(
  tempe_all_data[1, -c(1, 2)]
)

# Rows 2--4 correspond to Hong Kong, Shanghai, and Beijing.
tempe_all <- tempe_all_data[
  2:4,
  ,
  drop = FALSE
]

nu_all <- read.csv(
  here::here("Data", "tempe_all_std.csv"),
  check.names = FALSE
)

if (nrow(tempe_all) < 3L) {
  stop(
    "tempe_all_v3.csv must contain rows for ",
    "Hong Kong, Shanghai, and Beijing."
  )
}

if (nrow(nu_all) < 3L) {
  stop(
    "tempe_all_std.csv must contain rows for ",
    "Hong Kong, Shanghai, and Beijing."
  )
}

# ------------------------------------------------------------
# 2. City-specific configuration
# ------------------------------------------------------------

city_config <- list(
  HongKong = list(
    location_row = 1L,
    lme_file = here::here(
      "Data",
      "LME data",
      "d1.csv"
    ),
    lme_temperature_limits = c(20, 24),
    figure8_panel = "b"
  ),
  Shanghai = list(
    location_row = 2L,
    lme_file = here::here(
      "Data",
      "LME data",
      "d2.csv"
    ),
    lme_temperature_limits = c(13, 18),
    figure8_panel = "c"
  ),
  Beijing = list(
    location_row = 3L,
    lme_file = here::here(
      "Data",
      "LME data",
      "d3.csv"
    ),
    lme_temperature_limits = c(9, 14),
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
    stop(
      "yhat and nu must have the same length."
    )
  }

  if (any(!is.finite(yhat))) {
    stop("yhat contains non-finite values.")
  }

  if (any(!is.finite(nu))) {
    stop("nu contains non-finite values.")
  }

  # Prevent division by zero.
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

  # Select the LME CDF bandwidth by least-squares
  # cross-validation.
  fx_bw <- npudistbw(
    dat = c(x)
  )

  qx_min <- min(x, na.rm = TRUE) -
    5 * sd(x, na.rm = TRUE)

  qx_max <- max(x, na.rm = TRUE) +
    5 * sd(x, na.rm = TRUE)

  y_seq <- seq(
    ymin,
    ymax,
    length.out = 150
  )

  FY_values <- FY_hat(
    y_seq,
    yhat,
    std
  )

  FX_values <- Fx_inv(
    FY_values,
    qx_min,
    qx_max,
    fx_bw
  )

  ycorrected <- Fx_inv(
    FY_hat(xhat, yhat, std),
    qx_min,
    qx_max,
    fx_bw
  )

  list(
    y_seq = y_seq,
    FX_values = FX_values,
    ycorrected = ycorrected
  )
}

# ------------------------------------------------------------
# 4. Save one ggplot object
# ------------------------------------------------------------

save_plot <- function(
    plot_object,
    filename,
    width,
    height) {

  ggsave(
    filename = filename,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )

  message("Saved: ", filename)

  invisible(filename)
}

# ------------------------------------------------------------
# 5. Generate Figures 7--8 for one city
# ------------------------------------------------------------

make_city_figures <- function(
    city_name,
    config,
    output_dir = output_dir) {

  location_row <- config$location_row
  lme_file <- config$lme_file
  lme_limits <- config$lme_temperature_limits
  figure8_panel <- config$figure8_panel

  if (!file.exists(lme_file)) {
    stop(
      "LME input file for ",
      city_name,
      " was not found: ",
      lme_file
    )
  }

  # ----------------------------------------------------------
  # Read city-specific LME data
  # ----------------------------------------------------------

  lme_data <- read.csv(
    lme_file,
    row.names = 1,
    check.names = FALSE
  )

  if (ncol(lme_data) < 3L) {
    stop(
      "The LME input for ",
      city_name,
      " must contain two metadata columns followed by ",
      "yearly temperature columns."
    )
  }

  # First two columns contain metadata.
  lme_temperature <- as.matrix(
    lme_data[
      ,
      -c(1, 2),
      drop = FALSE
    ]
  ) - 273.15

  if (!is.numeric(lme_temperature)) {
    stop(
      "The LME temperature values for ",
      city_name,
      " are not numeric."
    )
  }

  # Mean across the 13 LME simulations for each year.
  lme_mean <- colMeans(
    lme_temperature,
    na.rm = TRUE
  )

  # Pooled historical LME sample used to estimate F_X.
  lme_sample <- as.numeric(
    t(lme_temperature)
  )

  # ----------------------------------------------------------
  # Extract city-specific REACHES prediction and uncertainty
  # ----------------------------------------------------------

  tempe_use <- tempe_all[
    location_row,
    ,
    drop = FALSE
  ]

  yhat <- as.numeric(
    tempe_use[
      ,
      -c(1, 2),
      drop = TRUE
    ]
  )

  nu <- as.numeric(
    nu_all[
      location_row,
      -c(1, 2),
      drop = TRUE
    ]
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

  valid_years <- is.finite(year3)

  city_years <- year3[valid_years]
  yhat <- yhat[valid_years]
  nu <- nu[valid_years]

  if (
    any(!is.finite(yhat)) ||
      any(!is.finite(nu))
  ) {
    stop(
      "Non-finite REACHES predictions or uncertainty values ",
      "were found for ",
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

  ycorrected <- qm$ycorrected

  # ----------------------------------------------------------
  # Prepare LME annual means matching the REACHES years
  # ----------------------------------------------------------

  lme_years <- 1350:(
    1349 + length(lme_mean)
  )

  lme_year_index <- match(
    city_years,
    lme_years
  )

  if (anyNA(lme_year_index)) {
    stop(
      "One or more REACHES years for ",
      city_name,
      " fall outside the available LME period."
    )
  }

  lme_mean_matching <- data.frame(
    temper = lme_mean[lme_year_index],
    year = city_years
  )

  # ----------------------------------------------------------
  # Figure 7(a): REACHES index time series
  # ----------------------------------------------------------

  df_reaches <- data.frame(
    reach = yhat,
    year = city_years
  )

  p_figure7a <- ggplot(
    data = df_reaches,
    aes(x = year, y = reach)
  ) +
    geom_line(
      colour = "darkorchid1"
    ) +
    labs(
      x = "Year",
      y = "Level"
    ) +
    theme(
      text = element_text(size = 19),
      legend.position = "right",
      legend.key.height = grid::unit(1.5, "cm"),
      plot.title = element_text(hjust = 0.5)
    )

  # ----------------------------------------------------------
  # Figure 7(b): transformation function
  # ----------------------------------------------------------

  df_transformation <- data.frame(
    index = qm$y_seq,
    temperature = qm$FX_values
  )

  p_figure7b <- ggplot(
    data = df_transformation,
    aes(x = index, y = temperature)
  ) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE
    ) +
    labs(
      x = "Index",
      y = "Temperature"
    ) +
    theme(
      text = element_text(size = 19),
      legend.position = "right",
      legend.key.height = grid::unit(1.5, "cm"),
      plot.title = element_text(hjust = 0.5)
    )

  # ----------------------------------------------------------
  # Figure 7(c): REACHES index distribution
  # ----------------------------------------------------------

  p_figure7c <- ggplot(
    df_reaches,
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
      x = "Temperature",
      y = "Density"
    ) +
    theme(
      text = element_text(size = 18),
      legend.position = "right",
      legend.key.height = grid::unit(1.5, "cm"),
      plot.title = element_text(hjust = 0.5)
    )

  # ----------------------------------------------------------
  # Figure 7(d): LME temperature distribution
  # ----------------------------------------------------------

  df_lme_distribution <- data.frame(
    temper = lme_sample
  )

  lme_breaks <- seq(
    lme_limits[1],
    lme_limits[2],
    by = 0.1
  )

  p_figure7d <- ggplot(
    data = df_lme_distribution,
    aes(x = temper)
  ) +
    geom_histogram(
      aes(y = after_stat(density)),
      breaks = lme_breaks,
      colour = "white",
      fill = "deepskyblue"
    ) +
    geom_density(
      colour = "black",
      linewidth = 0.4
    ) +
    scale_x_continuous(
      limits = lme_limits
    ) +
    labs(
      x = "Temperature",
      y = "Density"
    ) +
    theme(
      text = element_text(size = 18),
      legend.position = "right",
      legend.key.height = grid::unit(1.5, "cm"),
      plot.title = element_text(hjust = 0.5)
    )

  # ----------------------------------------------------------
  # Figure 8 city panel: corrected REACHES and LME means
  # ----------------------------------------------------------

  df_corrected <- data.frame(
    YEAR = city_years,
    temperature = ycorrected
  )

  p_figure8_city <- ggplot() +
    geom_line(
      data = df_corrected,
      aes(
        x = YEAR,
        y = temperature
      ),
      colour = "palevioletred1",
      alpha = 0.75
    ) +
    geom_smooth(
      data = df_corrected,
      aes(
        x = YEAR,
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
      data = lme_mean_matching,
      aes(
        x = year,
        y = temper
      ),
      colour = "skyblue1",
      alpha = 0.75
    ) +
    geom_smooth(
      data = lme_mean_matching,
      aes(
        x = year,
        y = temper
      ),
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 1.1,
      colour = "deepskyblue",
      linetype = "dashed"
    ) +
    labs(
      x = "Year",
      y = "Temperature"
    ) +
    theme(
      text = element_text(size = 12),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank()
    )

  # ----------------------------------------------------------
  # Figure 8(a): corrected REACHES with uncertainty
  # ----------------------------------------------------------

  df_uncertainty <- data.frame(
    year = city_years,
    REACHES = ycorrected,
    std = nu
  )

  p_figure8a <- ggplot(
    data = df_uncertainty,
    aes(x = year, y = REACHES)
  ) +
    geom_line(
      colour = "red"
    ) +
    geom_ribbon(
      aes(
        ymin = REACHES - std,
        ymax = REACHES + std
      ),
      alpha = 0.5,
      fill = "grey3"
    ) +
    labs(
      x = "Year",
      y = "Temperature"
    ) +
    theme(
      text = element_text(size = 11),
      legend.position = "right",
      legend.key.height = grid::unit(1.5, "cm"),
      plot.title = element_text(hjust = 0.5)
    )

  # ----------------------------------------------------------
  # Save city-specific outputs
  # ----------------------------------------------------------

  city_files <- c(
    Figure7a = file.path(
      output_dir,
      paste0("Figure7a_", city_name, ".png")
    ),
    Figure7b = file.path(
      output_dir,
      paste0("Figure7b_", city_name, ".png")
    ),
    Figure7c = file.path(
      output_dir,
      paste0("Figure7c_", city_name, ".png")
    ),
    Figure7d = file.path(
      output_dir,
      paste0("Figure7d_", city_name, ".png")
    ),
    Figure8a = file.path(
      output_dir,
      paste0("Figure8a_", city_name, ".png")
    ),
    Figure8City = file.path(
      output_dir,
      paste0(
        "Figure8",
        figure8_panel,
        "_",
        city_name,
        ".png"
      )
    )
  )

  save_plot(
    p_figure7a,
    city_files[["Figure7a"]],
    width = 6,
    height = 4
  )

  save_plot(
    p_figure7b,
    city_files[["Figure7b"]],
    width = 6,
    height = 4
  )

  save_plot(
    p_figure7c,
    city_files[["Figure7c"]],
    width = 6,
    height = 3.5
  )

  save_plot(
    p_figure7d,
    city_files[["Figure7d"]],
    width = 6,
    height = 3.5
  )

  save_plot(
    p_figure8a,
    city_files[["Figure8a"]],
    width = 6,
    height = 3
  )

  save_plot(
    p_figure8_city,
    city_files[["Figure8City"]],
    width = 6,
    height = 3
  )

  invisible(city_files)
}

# ------------------------------------------------------------
# 6. Run all three cities
# ------------------------------------------------------------

all_output_files <- lapply(
  names(city_config),
  function(city_name) {
    make_city_figures(
      city_name = city_name,
      config = city_config[[city_name]]
    )
  }
)

names(all_output_files) <- names(city_config)
