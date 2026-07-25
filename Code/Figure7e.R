here::i_am("Code/Figure7e.R")

# ============================================================
# Spatial quantile-mapping calibration for Figure 7(e)
#
# Required inputs:
#   Output/Intermediate/calibration_parameters.rds
#   Data/LME data/b0.csv
#
# Output:
#   Output/Figure7-8/Figure7e.png
# ============================================================

library(sp)
library(sf)
library(rnaturalearth)
library(dplyr)
library(stringr)
library(np)
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
# 1. Define the coarse grid and available REACHES years
# ------------------------------------------------------------

year_start <- min(temp2$year, na.rm = TRUE)
year_end <- max(temp2$year, na.rm = TRUE)
year_all <- year_start:year_end

ncase1 <- vapply(
  year_all,
  function(current_year) {
    sum(temp2$year == current_year)
  },
  numeric(1)
)

year2 <- year_all[ncase1 >= 1]

loc <- expand.grid(
  long = seq(97.5, 125, by = 2.5),
  lat = seq(
    18,
    42.63158,
    by = 1.89473692308
  )
)

coordinates(loc) <- ~ long + lat
proj4string(loc) <- CRS("+proj=longlat +datum=WGS84")

n_long <- length(unique(loc@coords[, 1]))
n_lat <- length(unique(loc@coords[, 2]))

if (n_long != 12L || n_lat != 14L) {
  stop(
    "Expected a 12 x 14 grid, but obtained ",
    n_long,
    " x ",
    n_lat,
    "."
  )
}

arr.pred <- arr.std <- array(
  NA_real_,
  dim = c(
    n_long,
    n_lat,
    length(year2)
  )
)

y2 <- var.fit2$y2

# ------------------------------------------------------------
# 2. Kriging helper functions
# ------------------------------------------------------------

# These definitions preserve the original calculation used in
# Figures 5 and 7(e).
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
    tol = 1e-8,
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
# 3. Krige all available years on the coarse grid
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

  if (i == 1L) {
    sk.all <- as.data.frame(loc@coords)
    sk.all1 <- as.data.frame(loc@coords)
  }

  sk.all <- cbind(
    sk.all,
    c(arr.pred[, , i])
  )

  sk.all1 <- cbind(
    sk.all1,
    c(arr.std[, , i])
  )
}

names(sk.all)[-c(1, 2)] <- year2
names(sk.all1)[-c(1, 2)] <- year2

tempe_all <- sk.all
nu_reaches <- sk.all1

# ------------------------------------------------------------
# 4. Keep locations in or near China, Taiwan, Hong Kong,
#    and Macao
# ------------------------------------------------------------

world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

targets <- world %>%
  filter(
    str_detect(
      admin,
      regex(
        "China|Taiwan|Hong Kong|Macao|Macau",
        ignore_case = TRUE
      )
    )
  )

china_tw_hk_mo <- st_union(targets)

china_tw_hk_mo <- st_as_sf(
  data.frame(
    name = "China+Taiwan+HK+MO"
  ),
  geometry = china_tw_hk_mo
)

st_crs(china_tw_hk_mo) <- 4326

coast_km <- 80

china_tw_hk_mo_valid <- st_make_valid(
  china_tw_hk_mo
)

poly_m <- st_transform(
  china_tw_hk_mo_valid,
  3857
)

poly_m_buf <- st_buffer(
  poly_m,
  dist = coast_km * 1000
)

china_tw_hk_mo_buf <- st_transform(
  poly_m_buf,
  4326
)

pts_sf <- st_as_sf(
  tempe_all,
  coords = c("long", "lat"),
  crs = 4326
)

keep_location <- st_within(
  pts_sf,
  china_tw_hk_mo_buf,
  sparse = FALSE
)[, 1]

tempe_all1 <- tempe_all[
  keep_location,
  ,
  drop = FALSE
]

nu_reaches1 <- nu_reaches[
  keep_location,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 5. Load LME data
# ------------------------------------------------------------

lme_file <- here::here(
  "Data",
  "LME data",
  "b0.csv"
)

if (!file.exists(lme_file)) {
  stop(
    "The LME file required for Figure 7(e) was not found: ",
    lme_file,
    "\nPlace b0.csv under Data/LME data/ or update lme_file."
  )
}

lme_df <- read.csv(
  lme_file,
  row.names = 1,
  check.names = FALSE
)

# ------------------------------------------------------------
# 6. Spatial quantile mapping at zero
# ------------------------------------------------------------

n_loc <- 121L
n_rep <- 13L

if (nrow(lme_df) < n_loc * n_rep) {
  stop(
    "Expected at least ",
    n_loc * n_rep,
    " LME rows, but found ",
    nrow(lme_df),
    "."
  )
}

get_block <- function(df, r) {
  first_row <- (r - 1L) * n_loc + 1L
  last_row <- r * n_loc

  df[
    first_row:last_row,
    ,
    drop = FALSE
  ]
}

get_xs_all_years <- function(lme_df, s) {
  xs_list <- lapply(
    seq_len(n_rep),
    function(r) {
      blk <- get_block(lme_df, r)
      as.numeric(
        blk[s, -(1:2)]
      )
    }
  )

  xs_kelvin <- unlist(
    xs_list,
    use.names = FALSE
  )

  xs_kelvin - 273.15
}

build_Fx_inv_local <- function(x_s) {
  x_s <- as.numeric(x_s)
  x_s <- x_s[is.finite(x_s)]

  bw <- npudistbw(dat = x_s)

  Fx_hat <- function(q) {
    fit <- npudist(
      bws = bw,
      edat = data.frame(x = q)
    )

    as.numeric(fitted(fit))
  }

  Fx_inv <- function(u) {
    u <- pmin(
      pmax(u, 1e-8),
      1 - 1e-8
    )

    xmin <- min(x_s) - 5 * sd(x_s)
    xmax <- max(x_s) + 5 * sd(x_s)

    sapply(
      u,
      function(uu) {
        uniroot(
          function(q) {
            Fx_hat(q) - uu
          },
          interval = c(xmin, xmax),
          tol = 1e-6
        )$root
      }
    )
  }

  Fx_inv
}

FY_hat_1loc <- function(y, yhat, nu) {
  yhat <- as.numeric(yhat)
  nu <- pmax(
    as.numeric(nu),
    1e-8
  )

  stopifnot(
    length(yhat) ==
      length(nu)
  )

  mean(
    pnorm(
      (y - yhat) / nu
    )
  )
}

lme_base <- get_block(
  lme_df,
  1
) %>%
  transmute(
    lat = lati,
    lon = long
  )

reaches_base <- tempe_all1 %>%
  transmute(
    lon = long,
    lat = lat
  )

coord_key_lme <- lme_base %>%
  mutate(
    key = paste0(
      round(lon, 3),
      "_",
      round(lat, 3)
    )
  )

coord_key_rea <- reaches_base %>%
  mutate(
    key = paste0(
      round(lon, 3),
      "_",
      round(lat, 3)
    )
  )

idx_map <- match(
  coord_key_lme$key,
  coord_key_rea$key
)

if (anyNA(idx_map)) {
  stop(
    "At least one LME grid location could not be matched ",
    "to the kriged REACHES grid."
  )
}

common_year_columns <- intersect(
  names(tempe_all1)[-c(1, 2)],
  names(nu_reaches1)[-c(1, 2)]
)

if (length(common_year_columns) == 0L) {
  stop(
    "No common year columns were found between the ",
    "kriged means and variances."
  )
}

# npudistbw uses a numerical multistart search, so fix the
# random sequence for repeatable bandwidth selection.
set.seed(10)

g0_121 <- sapply(
  seq_len(n_loc),
  function(s) {
    x_s <- get_xs_all_years(
      lme_df,
      s
    )

    Fx_inv_s <- build_Fx_inv_local(
      x_s
    )

    i_rea <- idx_map[s]

    yhat_s <- as.numeric(
      tempe_all1[
        i_rea,
        common_year_columns,
        drop = TRUE
      ]
    )

    nu_s <- as.numeric(
      nu_reaches1[
        i_rea,
        common_year_columns,
        drop = TRUE
      ]
    )

    u0 <- FY_hat_1loc(
      y = 0,
      yhat = yhat_s,
      nu = nu_s
    )

    Fx_inv_s(u0)
  }
)

plot_df <- lme_base %>%
  mutate(
    temp0_c = g0_121
  )

# ------------------------------------------------------------
# 7. Save Figure 7(e)
# ------------------------------------------------------------

p_figure7e <- ggplot(
  plot_df,
  aes(lon, lat)
) +
  geom_point(
    aes(colour = temp0_c),
    cex = 8.99,
    shape = 15
  ) +
  coord_map(
    xlim = c(98, 124.5),
    ylim = c(18, 42.5)
  ) +
  scale_color_gradientn(
    colors = c(
      "blue",
      "cyan",
      "green",
      "yellow",
      "red"
    ),
    limits = c(-10, 25),
    na.value = "transparent",
    guide = "colourbar"
  ) +
  borders(
    database = "world",
    xlim = c(76, 126),
    ylim = c(18, 45),
    fill = NA,
    colour = "grey50"
  ) +
  theme(
    text = element_text(size = 15),
    legend.title = element_blank(),
    legend.position = "right"
  )

figure7e_output_dir <- here::here(
  "Output",
  "Figure7-8"
)

dir.create(
  figure7e_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

figure7e_output_file <- file.path(
  figure7e_output_dir,
  "Figure7e.png"
)

ggsave(
  filename = figure7e_output_file,
  plot = p_figure7e,
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)

message(
  "Figure 7(e) saved to: ",
  figure7e_output_file
)
