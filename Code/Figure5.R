here::i_am("Code/Figure5.R")

# ============================================================
# REACHES observations and kriging results for Figure 5
#
# Required input:
#   Output/Intermediate/calibration_parameters.rds
#
# Outputs:
#   Output/Figure5/Figure5(a1).png
#   Output/Figure5/Figure5(a2).png
#   Output/Figure5/Figure5(a3).png
#   Output/Figure5/Figure5(b1).png
#   Output/Figure5/Figure5(b2).png
#   Output/Figure5/Figure5(b3).png
# ============================================================

library(sp)
library(RColorBrewer)
library(ggplot2)
library(mvtnorm)
library(maps)

calibration_file <- here::here(
  "Output",
  "Intermediate",
  "calibration_parameters.rds"
)

if (!file.exists(calibration_file)) {
  stop(
    "Required calibration file was not found: ",
    calibration_file,
    "\nRun Code/prepare_calibration.R first."
  )
}

calibration_results <- readRDS(calibration_file)

temp2 <- calibration_results$temp2
var.fit2 <- calibration_results$var_fit2
vario.fit2 <- calibration_results$vario_fit2

# ------------------------------------------------------------
# 1. Select years and define the 0.5 x 0.5 grid
# ------------------------------------------------------------

year_start <- min(temp2$year, na.rm = TRUE)
year_end <- max(temp2$year, na.rm = TRUE)
year_all <- year_start:year_end

year_positions <- c(98L, 484L)

if (max(year_positions) > length(year_all)) {
  stop(
    "The available year range is too short for the original ",
    "Figure 5 year positions 98 and 484."
  )
}

year2 <- year_all[year_positions]

loc <- expand.grid(
  long = seq(98.25, 124.25, by = 0.5),
  lat = seq(18.25, 42.25, by = 0.5)
)

coordinates(loc) <- ~ long + lat
proj4string(loc) <- CRS("+proj=longlat +datum=WGS84")

n_long <- length(unique(loc@coords[, 1]))
n_lat <- length(unique(loc@coords[, 2]))

arr.pred <- arr.std <- array(
  NA_real_,
  dim = c(
    n_long,
    n_lat,
    length(year2)
  )
)

y2 <- var.fit2$y2

figure5_output_dir <- here::here(
  "Output",
  "Figure5"
)

dir.create(
  figure5_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

save_figure5_panel <- function(
    plot_object,
    filename,
    width = 5,
    height = 4,
    dpi = 300) {

  output_file <- file.path(
    figure5_output_dir,
    filename
  )

  ggsave(
    filename = output_file,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )

  invisible(output_file)
}

# ------------------------------------------------------------
# 2. Figure 5(a1) and Figure 5(b1):
#    original REACHES observations
# ------------------------------------------------------------

panel_rows <- c("a", "b")

for (i in seq_along(year2)) {
  temp.dat <- as.data.frame(
    y2[y2$year == year2[i], ]
  )

  p_observed <- ggplot(
    temp.dat,
    aes(long, lat)
  ) +
    geom_point(
      aes(colour = level),
      cex = 2
    ) +
    ggtitle(
      paste0("Year ", year2[i])
    ) +
    coord_map(
      xlim = c(98, 124.5),
      ylim = c(18, 42.5)
    ) +
    scale_colour_gradientn(
      colours = rev(
        brewer.pal(
          n = 9,
          name = "RdBu"
        )
      ),
      limits = c(-2, 2),
      na.value = "transparent",
      guide = "colourbar"
    ) +
    borders(
      database = "world",
      xlim = c(76, 132),
      ylim = c(19, 52),
      fill = NA,
      colour = "grey50"
    ) +
    theme(
      text = element_text(size = 15),
      legend.title = element_blank(),
      legend.position = c(1.15, 0.7),
      legend.key.height = grid::unit(
        2.4,
        "cm"
      )
    )

  save_figure5_panel(
    plot_object = p_observed,
    filename = paste0(
      "Figure5(",
      panel_rows[i],
      "1).png"
    )
  )
}

# ------------------------------------------------------------
# 3. Kriging helper functions
# ------------------------------------------------------------

# These definitions preserve the original calculation used in
# Figure 5. The square-root transformation is intentionally kept
# because it is part of the original analysis code.
sigmay <- sqrt(vario.fit2$psill2)
sigmae <- sqrt(vario.fit2$psill1)
alpha <- vario.fit2$range / 100

cuts <- c(-Inf, -1.5, -0.5, 0.5, Inf)
vals <- c(-2, -1, 0, 1)

Ez_discrete <- function(sigma, cuts, vals) {
  stopifnot(length(vals) == length(cuts) - 1)

  a <- head(cuts, -1) / sigma
  b <- tail(cuts, -1) / sigma
  probs <- pnorm(b) - pnorm(a)

  sum(vals * probs)
}

EZstar_h <- function(sigma, cuts, vals) {
  stopifnot(length(vals) == length(cuts) - 1)

  a <- head(cuts, -1) / sigma
  b <- tail(cuts, -1) / sigma
  part <- sigma * (dnorm(a) - dnorm(b))

  sum(vals * part)
}

cov_Z_pair <- function(
    rho,
    sigma,
    cuts,
    vals,
    Eh = NULL) {

  if (is.null(Eh)) {
    Eh <- Ez_discrete(
      sigma,
      cuts,
      vals
    )
  }

  K <- length(vals)
  a_std <- head(cuts, -1) / sigma
  b_std <- tail(cuts, -1) / sigma

  Sigma2 <- matrix(
    c(1, rho, rho, 1),
    2,
    2
  )

  Eh2 <- 0

  for (k in seq_len(K)) {
    for (l in seq_len(K)) {
      lower <- c(
        a_std[k],
        a_std[l]
      )
      upper <- c(
        b_std[k],
        b_std[l]
      )

      pij <- as.numeric(
        pmvnorm(
          lower = lower,
          upper = upper,
          mean = c(0, 0),
          sigma = Sigma2
        )
      )

      Eh2 <- Eh2 +
        vals[k] * vals[l] * pij
    }
  }

  Eh2 - Eh^2
}

cZY_vector <- function(
    s_coords,
    s0,
    sigma_Y2,
    sigma_E2,
    alpha,
    cuts = c(-Inf, -1.5, -0.5, 0.5, Inf),
    vals = c(-2, -1, 0, 1)) {

  sigma2 <- sigma_Y2 + sigma_E2
  sigma <- sqrt(sigma2)

  Ezstar_h <- EZstar_h(
    sigma,
    cuts,
    vals
  )

  dists <- sqrt(
    rowSums(
      (
        s_coords -
          matrix(
            s0,
            nrow(s_coords),
            ncol(s_coords),
            byrow = TRUE
          )
      )^2
    )
  )

  covYY <- sigma_Y2 *
    exp(-dists / alpha)

  (covYY / sigma2) * Ezstar_h
}

SigmaZ_matrix <- function(
    s_coords,
    sigma_Y2,
    sigma_E2,
    alpha,
    cuts = c(-Inf, -1.5, -0.5, 0.5, Inf),
    vals = c(-2, -1, 0, 1)) {

  n <- nrow(s_coords)
  sigma2 <- sigma_Y2 + sigma_E2
  sigma <- sqrt(sigma2)

  Eh <- Ez_discrete(
    sigma,
    cuts,
    vals
  )

  dmat <- as.matrix(
    dist(
      s_coords,
      method = "euclidean",
      upper = TRUE,
      diag = TRUE
    )
  )

  rho <- (
    sigma_Y2 *
      exp(-dmat / alpha)
  ) / sigma2

  diag(rho) <- 1

  Sig <- matrix(
    NA_real_,
    n,
    n
  )

  for (i in seq_len(n)) {
    for (j in i:n) {
      cij <- cov_Z_pair(
        rho[i, j],
        sigma,
        cuts,
        vals,
        Eh = Eh
      )

      Sig[i, j] <- cij

      if (j != i) {
        Sig[j, i] <- cij
      }
    }
  }

  Sig
}

y_pred <- function(
    s_coords,
    s0,
    z_obs = NULL,
    sigma_Y2,
    sigma_E2,
    alpha,
    cuts = c(-Inf, -1.5, -0.5, 0.5, Inf),
    vals = c(-2, -1, 0, 1),
    tol = 1e-2,
    return = c("mean", "var")) {

  return <- unique(return)

  if (!is.null(z_obs)) {
    stopifnot(
      length(z_obs) ==
        nrow(s_coords)
    )
  }

  c_zy <- cZY_vector(
    s_coords,
    s0,
    sigma_Y2,
    sigma_E2,
    alpha,
    cuts,
    vals
  )

  SigmaZ <- SigmaZ_matrix(
    s_coords,
    sigma_Y2,
    sigma_E2,
    alpha,
    cuts,
    vals
  )

  SigmaZ <- (
    SigmaZ + t(SigmaZ)
  ) / 2

  out <- list()

  if ("mean" %in% return) {
    sigma2 <- sigma_Y2 + sigma_E2
    sigma <- sqrt(sigma2)

    Ez <- Ez_discrete(
      sigma,
      cuts,
      vals
    )

    rhs <- z_obs - Ez
    w_mean <- solve(
      SigmaZ,
      rhs,
      tol = tol
    )

    out$mean <- drop(
      crossprod(
        c_zy,
        w_mean
      )
    )
  }

  if ("var" %in% return) {
    w_var <- solve(
      SigmaZ,
      c_zy,
      tol = tol
    )

    var_pred <- sigma_Y2 -
      drop(
        crossprod(
          c_zy,
          w_var
        )
      )

    out$var <- max(
      var_pred,
      0
    )
  }

  if (length(out) == 1) {
    return(out[[1]])
  }

  out
}

# ------------------------------------------------------------
# 4. Figure 5(a2), 5(a3), 5(b2), and 5(b3)
# ------------------------------------------------------------

for (i in seq_along(year2)) {
  temp14 <- y2[
    y2$year == year2[i],
  ]

  temp15 <- numeric(
    nrow(loc@coords)
  )
  temp16 <- numeric(
    nrow(loc@coords)
  )

  locations <- temp14@coords

  for (j in seq_len(nrow(loc@coords))) {
    s0 <- as.numeric(
      loc@coords[j, ]
    )

    res <- y_pred(
      s_coords = locations,
      s0 = s0,
      z_obs = as.matrix(temp14$level),
      sigma_Y2 = sigmay,
      sigma_E2 = sigmae,
      alpha = alpha,
      return = c("mean", "var")
    )

    temp15[j] <- res$mean
    temp16[j] <- res$var
  }

  arr.pred[, , i] <- matrix(
    temp15,
    nrow = n_long,
    ncol = n_lat
  )

  arr.std[, , i] <- matrix(
    temp16,
    nrow = n_long,
    ncol = n_lat
  )

  temp.dat1 <- cbind(
    as.data.frame(loc@coords)[, 1:2],
    mu = c(arr.pred[, , i]),
    std = c(arr.std[, , i])
  )

  p_mean <- ggplot(
    temp.dat1,
    aes(long, lat)
  ) +
    geom_point(
      aes(colour = mu),
      cex = 9
    ) +
    ggtitle(
      paste0("Year ", year2[i])
    ) +
    coord_map(
      xlim = c(98, 124.5),
      ylim = c(18, 42.5)
    ) +
    scale_colour_gradientn(
      colours = rev(
        brewer.pal(
          n = 9,
          name = "RdBu"
        )
      ),
      limits = c(-2, 2),
      na.value = "transparent",
      guide = "colourbar"
    ) +
    borders(
      database = "world",
      xlim = c(76, 132),
      ylim = c(19, 52),
      fill = NA,
      colour = "grey50"
    ) +
    theme(
      text = element_text(size = 15),
      legend.title = element_blank(),
      legend.position = c(1.15, 0.7),
      legend.key.height = grid::unit(
        2.4,
        "cm"
      )
    )

  p_std <- ggplot(
    temp.dat1,
    aes(long, lat)
  ) +
    geom_point(
      aes(colour = std),
      cex = 9
    ) +
    ggtitle(
      paste0("Year ", year2[i])
    ) +
    coord_map(
      xlim = c(98, 124.5),
      ylim = c(18, 42.5)
    ) +
    scale_colour_gradient(
      low = "skyblue3",
      high = "white",
      limits = c(0, 1),
      na.value = "transparent",
      guide = "colourbar"
    ) +
    borders(
      database = "world",
      xlim = c(76, 132),
      ylim = c(19, 52),
      fill = NA,
      colour = "grey50"
    ) +
    theme(
      text = element_text(size = 15),
      legend.title = element_blank(),
      legend.position = c(1.15, 0.5),
      legend.key.height = grid::unit(
        1.745,
        "cm"
      )
    )

  save_figure5_panel(
    plot_object = p_mean,
    filename = paste0(
      "Figure5(",
      panel_rows[i],
      "2).png"
    )
  )

  save_figure5_panel(
    plot_object = p_std,
    filename = paste0(
      "Figure5(",
      panel_rows[i],
      "3).png"
    )
  )
}

message(
  "All Figure 5 panels were saved to: ",
  figure5_output_dir
)
