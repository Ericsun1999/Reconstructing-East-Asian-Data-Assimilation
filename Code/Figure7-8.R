here::i_am("Code/Figure7-8.R")

# ============================================================
# Quantile-mapping calibration for Figures 7 and 8
#
# Outputs:
#   Figure 7(a)--(d): Beijing
#   Figure 8(a): Beijing
#   Figure 8(b): Beijing
#   Figure 8(c): Shanghai
#   Figure 8(d): Hong Kong
#
# Figure 7(e) is generated separately by Code/Figure7e.R.
#
# City ordering:
#   Data1 / row 1 = Hong Kong
#   Data2 / row 2 = Shanghai
#   Data3 / row 3 = Beijing
# ============================================================

library(ggplot2)
library(np)

# ------------------------------------------------------------
# 1. Read kriged REACHES data and uncertainties
# ------------------------------------------------------------

tempe_all_data <- read.csv(
  here::here("Data", "tempe_all_v3.csv"),
  header = FALSE
)

year3 <- as.integer(
  tempe_all_data[1, -c(1, 2)]
)

tempe_all <- tempe_all_data[
  2:4,
  ,
  drop = FALSE
]

nu1 <- read.csv(
  here::here("Data", "tempe_all_std.csv")
)

# ------------------------------------------------------------
# 2. Read city-specific LME data
# ------------------------------------------------------------

Data1 <- read.csv(
  here::here("Data", "LME data", "d1.csv"),
  row.names = 1
)

Data2 <- read.csv(
  here::here("Data", "LME data", "d2.csv"),
  row.names = 1
)

Data3 <- read.csv(
  here::here("Data", "LME data", "d3.csv"),
  row.names = 1
)

# ------------------------------------------------------------
# 3. Quantile-mapping functions
# ------------------------------------------------------------

# Estimate F_X.
Fx_hat <- function(q, fx_bw) {
  fx_ob <- npudist(
    bws = fx_bw,
    edat = data.frame(x = q)
  )

  fitted(fx_ob)
}

# Estimate F_Y.
FY_hat <- function(y, yhat, nu) {
  yhat <- as.numeric(yhat)
  nu <- as.numeric(nu)

  stopifnot(length(yhat) == length(nu))

  # Avoid division by zero.
  nu_safe <- pmax(nu, 1e-8)

  sapply(
    y,
    function(yy) {
      mean(
        pnorm((yy - yhat) / nu_safe)
      )
    }
  )
}

# Numerically invert F_X.
Fx_inv <- function(
    u,
    qx_min,
    qx_max,
    fx_bw) {

  u <- pmin(
    pmax(u, 1e-8),
    1 - 1e-8
  )

  sapply(
    u,
    function(uu) {
      f <- function(q) {
        Fx_hat(q, fx_bw) - uu
      }

      uniroot(
        f,
        interval = c(qx_min, qx_max),
        tol = 1e-6
      )$root
    }
  )
}

qmapping <- function(
    yhat = yhat,
    std = nu,
    x = haave2,
    ymax = 1,
    ymin = -2,
    xhat = tempe_use[, -c(1, 2)]) {

  # Estimate F_X using npudist with bandwidth chosen by LS-CV.
  fx_bw <- npudistbw(
    dat = c(x)
  )

  # Inverse of F_X using uniroot.
  qx_min <- min(x) - 5 * sd(x)
  qx_max <- max(x) + 5 * sd(x)

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
    FY_hat(
      as.numeric(xhat),
      yhat,
      std
    ),
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
# 4. Prepare one city's quantile-mapping results
#
# This preserves the calculations used in the original script.
# ------------------------------------------------------------

prepare_city <- function(
    city_name,
    DData,
    loc) {

  tempe_use <- tempe_all[
    loc,
    ,
    drop = FALSE
  ]

  nu <- nu1[
    loc,
    -c(1, 2)
  ]

  haa <- (
    DData[, -c(1, 2), drop = FALSE]
  ) - 273.15

  haave <- colMeans(haa)
  haave2 <- as.numeric(t(haa))

  yhat <- as.numeric(
    as.matrix(
      tempe_use[, -c(1, 2), drop = FALSE]
    )
  )

  nu <- as.numeric(
    as.matrix(nu)
  )

  qm <- qmapping(
    yhat = yhat,
    std = nu,
    x = haave2,
    ymax = 1,
    ymin = -2,
    xhat = tempe_use[, -c(1, 2), drop = FALSE]
  )

  list(
    city_name = city_name,
    loc = loc,
    tempe_use = tempe_use,
    nu = nu,
    haa = haa,
    haave = haave,
    haave2 = haave2,
    yhat = yhat,
    qm = qm,
    ycorrected = qm$ycorrected
  )
}

# ------------------------------------------------------------
# 5. Automatically prepare all three cities
# ------------------------------------------------------------

city_config <- list(
  Beijing = list(
    DData = Data3,
    loc = 3L
  ),
  Shanghai = list(
    DData = Data2,
    loc = 2L
  ),
  HongKong = list(
    DData = Data1,
    loc = 1L
  )
)

city_results <- lapply(
  names(city_config),
  function(city_name) {
    config <- city_config[[city_name]]

    prepare_city(
      city_name = city_name,
      DData = config$DData,
      loc = config$loc
    )
  }
)

names(city_results) <- names(city_config)

beijing <- city_results[["Beijing"]]
shanghai <- city_results[["Shanghai"]]
hongkong <- city_results[["HongKong"]]

# ------------------------------------------------------------
# 6. Figure 7(a)--(d): Beijing only
# ------------------------------------------------------------

df111 <- data.frame(
  reach = as.numeric(
    beijing$tempe_use[, -c(1, 2)]
  ),
  year = c(year3)
)

# Figure 7(a)
p_figure7a <- ggplot(
  data = df111,
  aes(x = year, y = reach)
) +
  geom_line(
    color = "darkorchid1"
  ) +
  theme(
    text = element_text(size = 19),
    legend.position = "right",
    legend.key.height = grid::unit(1.5, "cm"),
    plot.title = element_text(hjust = 0.5)
  ) +
  xlab("year") +
  ylab("level")

# Figure 7(b)
df_transformation <- data.frame(
  z = beijing$qm$y_seq,
  y1 = beijing$qm$FX_values
)

p_figure7b <- ggplot(
  data = df_transformation,
  aes(x = z, y = y1)
) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE
  ) +
  xlab("index") +
  ylab("temperature") +
  theme(
    text = element_text(size = 19),
    legend.position = "right",
    legend.key.height = grid::unit(1.5, "cm"),
    plot.title = element_text(hjust = 0.5)
  )

# Figure 7(c)
p_figure7c <- ggplot(
  df111,
  aes(x = reach)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = seq(
      -2.5,
      1.5,
      by = 0.2
    ),
    binwidth = 0.1,
    fill = "darkorchid1",
    color = "white"
  ) +
  geom_density(
    color = "black",
    size = 0.4
  ) +
  theme(
    text = element_text(size = 18),
    legend.position = "right",
    legend.key.height = grid::unit(1.5, "cm"),
    plot.title = element_text(hjust = 0.5)
  ) +
  xlab("temperature") +
  ylab("density")

# Figure 7(d)
dsnorm.fit2 <- data.frame(
  temper = c(beijing$haave2)
)

p_figure7d <- ggplot(
  data = dsnorm.fit2,
  aes(x = temper)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = seq(
      9,
      14,
      by = 0.1
    ),
    binwidth = 0.1,
    colour = "white",
    fill = "deepskyblue"
  ) +
  geom_density(
    color = "black",
    size = 0.4
  ) +
  theme(
    text = element_text(size = 18),
    legend.position = "right",
    legend.key.height = grid::unit(1.5, "cm"),
    plot.title = element_text(hjust = 0.5)
  ) +
  xlim(9, 14) +
  xlab("temperature")

# ------------------------------------------------------------
# 7. Figure 8(a): Beijing only
# ------------------------------------------------------------

df8 <- data.frame(
  year = c(year3),
  REACHES = beijing$ycorrected,
  std = beijing$nu
)

p_figure8a <- ggplot(
  data = df8,
  aes(x = year, y = REACHES)
) +
  # Preserve the original layer order.
  geom_line(
    color = "red"
  ) +
  geom_ribbon(
    aes(
      ymin = REACHES - std,
      ymax = REACHES + std
    ),
    alpha = 0.5,
    fill = "grey3"
  ) +
  theme(
    text = element_text(size = 11),
    legend.position = "right",
    legend.key.height = grid::unit(1.5, "cm"),
    plot.title = element_text(hjust = 0.5)
  ) +
  xlab("year") +
  ylab("temperature")

# ------------------------------------------------------------
# 8. Figure 8(b)--(d): Beijing, Shanghai, and Hong Kong
#
# This preserves the original year indexing and the first
# 524-row plotting range used to create the manuscript panels.
# ------------------------------------------------------------

make_figure8_city_plot <- function(city_result) {

  df7 <- data.frame(
    YEAR = c(year3),
    y = city_result$ycorrected
  )

  df7 <- cbind(
    df7,
    t(as.vector("REACH"))
  )

  colnames(df7) <- c(
    "YEAR",
    "temperature",
    "type"
  )

  dsnorm.fit1 <- data.frame(
    temper = c(city_result$haave)[
      year3 - 1350
    ],
    year = c(year3)
  )

  ggplot() +
    geom_line(
      data = df7[1:524, ],
      aes(
        x = YEAR,
        y = temperature
      ),
      color = "palevioletred1",
      alpha = 0.75
    ) +
    geom_smooth(
      data = df7[1:524, ],
      aes(
        x = YEAR,
        y = temperature
      ),
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      size = 1.1,
      color = "firebrick",
      linetype = "dashed",
      span = 0.75
    ) +
    geom_line(
      data = dsnorm.fit1[1:524, ],
      aes(
        x = year,
        y = temper
      ),
      color = "skyblue1",
      alpha = 0.75
    ) +
    geom_smooth(
      data = dsnorm.fit1[1:524, ],
      aes(
        x = year,
        y = temper
      ),
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      size = 1.1,
      color = "deepskyblue",
      linetype = "dashed"
    ) +
    xlab("year") +
    ylab("temperature") +
    theme(
      text = element_text(size = 12),
      legend.position = "bottom",
      legend.key.height = grid::unit(-0.5, "cm"),
      plot.title = element_text(hjust = 0.5)
    )
}

p_figure8b <- make_figure8_city_plot(beijing)
p_figure8c <- make_figure8_city_plot(shanghai)
p_figure8d <- make_figure8_city_plot(hongkong)

# ------------------------------------------------------------
# 9. Save all panels
# ------------------------------------------------------------

output_dir <- here::here(
  "Output",
  "Figure7-8"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

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

message(
  "All Figure 7--8 panels were saved to: ",
  output_dir
)
