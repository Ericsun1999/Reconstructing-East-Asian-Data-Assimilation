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

predict_grid_one_year <- function(
    current_year,
    observation_data,
    prediction_coordinates,
    sigma_Y2,
    sigma_E2,
    alpha,
    cuts = c(-Inf, -1.5, -0.5, 0.5, Inf),
    vals = c(-2, -1, 0, 1),
    tol = 1e-2) {

  temp14 <- observation_data[
    observation_data$year == current_year,
  ]

  if (nrow(temp14) == 0L) {
    stop(
      "No REACHES observations were found for year ",
      current_year,
      "."
    )
  }

  locations <- temp14@coords
  z_obs <- as.numeric(temp14$level)

  if (length(z_obs) != nrow(locations)) {
    stop(
      "The observation values and coordinates have different ",
      "lengths for year ",
      current_year,
      "."
    )
  }

  # SigmaZ depends on the observation locations and model
  # parameters, not on the prediction location. Construct it
  # only once for this year.
  SigmaZ <- SigmaZ_matrix(
    s_coords = locations,
    sigma_Y2 = sigma_Y2,
    sigma_E2 = sigma_E2,
    alpha = alpha,
    cuts = cuts,
    vals = vals
  )

  SigmaZ <- (
    SigmaZ + t(SigmaZ)
  ) / 2

  sigma2 <- sigma_Y2 + sigma_E2
  sigma <- sqrt(sigma2)

  Ez <- Ez_discrete(
    sigma = sigma,
    cuts = cuts,
    vals = vals
  )

  rhs_mean <- z_obs - Ez

  # Each column is the cross-covariance vector for one
  # prediction grid point.
  CZY <- vapply(
    seq_len(nrow(prediction_coordinates)),
    function(j) {

      s0 <- as.numeric(
        prediction_coordinates[j, ]
      )

      cZY_vector(
        s_coords = locations,
        s0 = s0,
        sigma_Y2 = sigma_Y2,
        sigma_E2 = sigma_E2,
        alpha = alpha,
        cuts = cuts,
        vals = vals
      )
    },
    numeric(nrow(locations))
  )

  if (!is.matrix(CZY)) {
    CZY <- matrix(
      CZY,
      nrow = nrow(locations)
    )
  }

  # Solve the mean system and all prediction-variance systems
  # together. This performs one matrix factorization per year.
  all_rhs <- cbind(
    rhs_mean,
    CZY
  )

  all_solutions <- tryCatch(
    solve(
      SigmaZ,
      all_rhs,
      tol = tol
    ),
    error = function(e) {
      stop(
        "Failed to solve the kriging systems for year ",
        current_year,
        ". Reciprocal condition number = ",
        signif(rcond(SigmaZ), 6),
        ". Original error: ",
        conditionMessage(e)
      )
    }
  )

  w_mean <- all_solutions[, 1]

  W_variance <- all_solutions[
    ,
    -1,
    drop = FALSE
  ]

  prediction_mean <- drop(
    crossprod(
      CZY,
      w_mean
    )
  )

  prediction_variance <- sigma_Y2 -
    colSums(
      CZY * W_variance
    )

  prediction_variance <- pmax(
    prediction_variance,
    0
  )

  list(
    mean = prediction_mean,
    variance = prediction_variance
  )
}

# ------------------------------------------------------------
# 4. Figure 5(a2), 5(a3), 5(b2), and 5(b3)
# ------------------------------------------------------------

for (i in seq_along(year2)) {

  message(
    "Kriging year ",
    year2[i],
    " (",
    i,
    "/",
    length(year2),
    ")"
  )

  temp14 <- y2[
    y2$year == year2[i],
  ]

  yearly_result <- predict_grid_one_year(
    current_year = year2[i],
    observation_data = y2,
    prediction_coordinates = loc@coords,
    sigma_Y2 = sigmay,
    sigma_E2 = sigmae,
    alpha = alpha,
    cuts = cuts,
    vals = vals,
    tol = 1e-2
  )

  temp15 <- yearly_result$mean
  temp16 <- yearly_result$variance

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
