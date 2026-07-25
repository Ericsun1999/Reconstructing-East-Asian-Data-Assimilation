here::i_am("Code/prepare_calibration.R")

# ============================================================
# Shared calibration analysis used by Figures 4, 5, and 7(e)
#
# Output:
#   Output/Intermediate/calibration_parameters.rds
# ============================================================

library(readxl)
library(dplyr)
library(sp)
library(spacetime)
library(zoo)
library(gstat)
library(mgcv)
library(ggplot2)
library(mvtnorm)

temperature <- read_excel(
  here::here("Data", "temperature index value.v1.xlsx"),
  col_types = c(
    "skip", "skip", "numeric", "numeric", "skip", "skip", "skip",
    "skip", "skip", "numeric", "numeric", "skip", "skip"
  )
)
colnames(temperature) <- c("level", "year", "long", "lat")

temp2 <- temperature %>%
  group_by(year, long, lat)

par_est_initial <- function(temp2) {
  y2 <- temp2[temp2$year > 100 & temp2$year <= 2000, ]

  coord <- unique(
    round(
      as.data.frame(cbind(long = y2$long, lat = y2$lat)),
      4
    )
  )
  coordinates(coord) <- ~ long + lat
  proj4string(coord) <- CRS("+proj=longlat +datum=WGS84")

  year <- sort(unique(y2$year))
  YM <- as.yearmon(year)

  temp23 <- as.data.frame(round(as.matrix(y2), 4))
  ind <- matrix(NA, nrow(temp23), 2)
  ind[, 2] <- match(temp23$year, year)
  ind[, 1] <- match(
    temp23$long * 10^7 + temp23$lat,
    coord@coords[, 1] * 10^7 + coord@coords[, 2]
  )

  stsdf2 <- STSDF(
    sp = coord,
    time = YM,
    data = y2,
    index = ind
  )

  vario.stsdf2 <- variogramST(
    level ~ 1,
    data = stsdf2,
    tlags = 0,
    width = 10,
    na.omit = TRUE
  )
  vario.stsdf3 <- subset(vario.stsdf2, dist > 1)

  coordinates(y2) <- ~ long + lat
  proj4string(y2) <- CRS("+proj=longlat +datum=WGS84")

  vario2 <- variogram(level ~ 1, data = y2, width = 10)

  ind <- match(
    as.numeric(rownames(vario.stsdf3)) - 1,
    rownames(vario2)
  )
  ind2 <- rep(FALSE, nrow(vario2))
  ind2[ind] <- TRUE
  vario2 <- vario2[ind2, ]

  if (sum(vario.stsdf3$dist == 0) > 0) {
    vario2[, 1:3] <- vario.stsdf3[ind[-1], 1:3]
  } else {
    vario2[, 1:3] <- vario.stsdf3[ind, 1:3]
  }

  vario.fit2 <- fit.variogram(
    vario2,
    vgm(model = "Exp", nugget = NA),
    fit.kappa = TRUE,
    fit.method = 2
  )

  list(
    psill = vario.fit2$psill,
    range = vario.fit2$range[2],
    y2 = y2
  )
}

findgam <- function(y, z, y_target) {
  df <- data.frame(vary = y, varz = z)
  gam_fit <- gam(varz ~ s(vary), data = df)

  x_values <- seq(
    min(df$vary),
    max(df$vary),
    length.out = 1000
  )
  predicted_y <- predict(
    gam_fit,
    newdata = data.frame(vary = x_values)
  )

  closest_index <- which.min(abs(predicted_y - y_target))
  x_closest <- x_values[closest_index]

  message(
    "The x-value corresponding to y = ",
    y_target,
    " is: ",
    x_closest
  )

  x_closest
}

par_estimation_calibrate <- function(
    psill1,
    psill,
    range,
    n1 = 500,
    n2 = 200,
    n3 = 40,
    plot = FALSE) {

  # Var(Y + epsilon) versus Var(Z)
  vary <- numeric(n1)
  varz <- numeric(n1)

  for (i in seq_len(n1)) {
    k <- 2 * i / n1
    variance_y <- k

    y <- rnorm(10000, 0, sqrt(variance_y))
    z <- round(y, digits = 0)
    z[z > 1] <- 1
    z[z < -2] <- -2

    vary[i] <- var(y)
    varz[i] <- var(z)
  }

  x_closest <- findgam(
    vary,
    varz,
    psill[1] + psill[2]
  )
  df <- data.frame(vary = vary, varz = varz)

  p1 <- ggplot(df, aes(x = vary, y = varz)) +
    geom_point(
      aes(x = x_closest, y = psill[1] + psill[2]),
      colour = "red",
      cex = 3
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 25)
    ) +
    xlab(expression(sigma[Y]^2 + sigma[epsilon]^2)) +
    ylab(
      expression(
        f[1] ~ "(" ~ sigma[Y]^2 + sigma[epsilon]^2 ~ ")"
      )
    ) +
    ylim(0, 1.42) +
    xlim(0, 1.42) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "blue"
    ) +
    geom_abline(
      intercept = psill[1] + psill[2],
      slope = 0,
      linetype = "longdash",
      color = "red"
    ) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(x),
      se = FALSE,
      size = 0.7,
      color = "firebrick",
      linetype = "solid"
    )

  if (isTRUE(plot)) print(p1)

  # Var(epsilon) versus rounded Var(epsilon)
  sigmaepsilon <- numeric(n2)
  sigmaepsilon1 <- numeric(n2)

  for (i in seq_len(n2)) {
    k <- i * x_closest / n2
    variance_y <- x_closest - k

    y <- rnorm(10000, sd = sqrt(variance_y))
    e <- rmvnorm(
      10000,
      sigma = matrix(
        c(
          x_closest - variance_y, 0,
          0, x_closest - variance_y
        ),
        ncol = 2
      )
    )

    z1 <- round(y + e[, 1], digits = 0)
    z1[z1 > 1] <- 1
    z1[z1 < -2] <- -2

    z2 <- round(y + e[, 2], digits = 0)
    z2[z2 > 1] <- 1
    z2[z2 < -2] <- -2

    sigmaepsilon[i] <- var(e[, 1] - e[, 2]) / 2
    sigmaepsilon1[i] <- var(z1 - z2) / 2
  }

  x_closest1 <- findgam(
    sigmaepsilon,
    sigmaepsilon1,
    psill1
  )
  df1 <- data.frame(
    sigmaepsilon = sigmaepsilon,
    sigmaepsilon1 = sigmaepsilon1
  )

  p2 <- ggplot(
    df1,
    aes(x = sigmaepsilon, y = sigmaepsilon1)
  ) +
    geom_point(
      aes(x = x_closest1, y = psill1),
      colour = "red",
      cex = 3
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 25)
    ) +
    xlab(expression(sigma[epsilon]^2)) +
    ylab(expression(f[2] ~ "(" ~ sigma[epsilon]^2 ~ ")")) +
    ylim(0, 0.9) +
    xlim(0, 0.9) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "blue"
    ) +
    geom_abline(
      intercept = psill1,
      slope = 0,
      linetype = "longdash",
      color = "red"
    ) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(x),
      se = FALSE,
      size = 0.7,
      color = "firebrick",
      linetype = "solid"
    )

  if (isTRUE(plot)) print(p2)

  # alpha versus calibrated alpha
  alpha2 <- numeric(n3 + 1)
  alphastar <- numeric(n3 + 1)
  cz <- array(0, dim = c(241, 151))

  variance_y <- x_closest - x_closest1

  for (aaa in 1:151) {
    aa <- 95 + 5 * aaa

    for (hh in 0:240) {
      a <- array(
        c(0, 0, 0, hh * 5),
        dim = c(2, 2)
      )
      sigma22 <- matrix(0, nrow = 2, ncol = 2)

      for (k in 1:2) {
        for (j in 1:2) {
          sigma22[k, j] <-
            variance_y *
            exp(
              -sqrt(
                t(a[, k] - a[, j]) %*%
                  (a[, k] - a[, j])
              ) / aa
            )

          if (k == j) {
            sigma22[k, j] <- sigma22[k, j] + x_closest1
          }
        }
      }

      y <- rmvnorm(1000, sigma = sigma22)
      z1 <- round(y, digits = 0)
      z1[z1 > 1] <- 1
      z1[z1 < -2] <- -2

      cz[hh + 1, aaa] <- cov(z1[, 1], z1[, 2])
    }
  }

  kk <- 0

  for (i in 0:n3) {
    alpha1 <- 100 + i * 10
    sum2 <- 100000000
    kk <- kk + 1

    for (aaa in kk:151) {
      sum1 <- 0
      aa <- 95 + 5 * aaa

      for (hh in 0:240) {
        sum1 <- sum1 +
          (
            cz[hh + 1, aaa] -
              variance_y * exp(-(hh * 5) / alpha1)
          )^2
      }

      if (sum2 >= sum1) {
        sum2 <- sum1
        np <- aa
      }
    }

    kk <- (np - 95) / 5
    if (kk > 150) kk <- 150

    alphastar[i + 1] <- np
    alpha2[i + 1] <- alpha1
  }

  x_closest2 <- findgam(alpha2, alphastar, range)
  df2 <- data.frame(
    alphastar = alphastar,
    alpha2 = alpha2
  )

  p3 <- ggplot(
    df2,
    aes(x = alpha2, y = alphastar)
  ) +
    geom_point(
      aes(x = x_closest2, y = range),
      colour = "red",
      cex = 3
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 27)
    ) +
    xlab(expression(alpha)) +
    ylab(expression(f[3] ~ "(" ~ alpha ~ ")")) +
    ylim(90, 810) +
    xlim(90, 600) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "blue"
    ) +
    geom_abline(
      intercept = range,
      slope = 0,
      linetype = "longdash",
      color = "red"
    ) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(x),
      se = FALSE,
      size = 0.7,
      color = "firebrick",
      linetype = "solid"
    )

  if (isTRUE(plot)) print(p3)

  list(
    psill1 = x_closest1,
    psill2 = x_closest - x_closest1,
    range = x_closest2,
    plot1 = p1,
    plot2 = p2,
    plot3 = p3,
    calibration_data1 = df,
    calibration_data2 = df1,
    calibration_data3 = df2
  )
}

var.fit2 <- par_est_initial(temp2)

set.seed(10)

vario.fit2 <- par_estimation_calibrate(
  psill1 = var.fit2$psill[1],
  psill = var.fit2$psill,
  range = var.fit2$range,
  n1 = 500,
  n2 = 200,
  n3 = 40,
  plot = FALSE
)

intermediate_dir <- here::here("Output", "Intermediate")
dir.create(
  intermediate_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

calibration_file <- file.path(
  intermediate_dir,
  "calibration_parameters.rds"
)

saveRDS(
  list(
    temp2 = temp2,
    var_fit2 = var.fit2,
    vario_fit2 = vario.fit2
  ),
  calibration_file
)

message(
  "Calibration results saved to: ",
  calibration_file
)
