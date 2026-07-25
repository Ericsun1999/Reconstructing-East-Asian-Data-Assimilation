here::i_am("Code/Figure7e.R")

# ============================================================
# Quantile mapping and spatial smoothing for Figure 7(e)
#
# Inputs:
#   Data/reaches_kriging_grid12x14_mean.csv
#     Kriging means at the 12 x 14 coarse grid.
#
#   Data/reaches_kriging_grid12x14_variance.csv
#     Kriging prediction VARIANCES at the same grid and years.
#     They are converted to standard deviations before evaluating F_Y.
#
#   Data/LME data/b0.csv
#     LME temperatures for 121 locations and 13 ensemble members.
#
# Outputs:
#   Output/Figure7-8/Figure7e.png
#   Output/Figure7-8/Figure7e_qm_locations.csv
#   Output/Figure7-8/Figure7e_smoothed_grid.csv
# ============================================================

library(here)
library(dplyr)
library(sf)
library(rnaturalearth)
library(stringr)
library(np)
library(mgcv)
library(ggplot2)
library(maps)

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

mean_file <- here::here(
  "Data",
  "reaches_kriging_grid12x14_mean.csv"
)

variance_file <- here::here(
  "Data",
  "reaches_kriging_grid12x14_variance.csv"
)

lme_file <- here::here(
  "Data",
  "LME data",
  "b0.csv"
)

required_files <- c(mean_file, variance_file, lme_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "The following required files were not found:\n",
    paste0("  - ", missing_files, collapse = "\n")
  )
}

# ------------------------------------------------------------
# 2. Read and validate the kriging files
# ------------------------------------------------------------

drop_accidental_index_column <- function(df) {
  if (ncol(df) == 0L) {
    return(df)
  }

  first_name <- names(df)[1]
  looks_like_index_name <- first_name %in% c("", "X", "...1", "row.names")
  first_values <- suppressWarnings(as.integer(df[[1]]))

  looks_like_row_numbers <-
    length(first_values) == nrow(df) &&
    all(!is.na(first_values)) &&
    identical(first_values, seq_len(nrow(df)))

  if (looks_like_index_name && looks_like_row_numbers) {
    df <- df[, -1, drop = FALSE]
  }

  df
}

standardize_coordinate_names <- function(df) {
  if ("lon" %in% names(df) && !"long" %in% names(df)) {
    names(df)[names(df) == "lon"] <- "long"
  }
  if ("longitude" %in% names(df) && !"long" %in% names(df)) {
    names(df)[names(df) == "longitude"] <- "long"
  }
  if ("lati" %in% names(df) && !"lat" %in% names(df)) {
    names(df)[names(df) == "lati"] <- "lat"
  }
  if ("latitude" %in% names(df) && !"lat" %in% names(df)) {
    names(df)[names(df) == "latitude"] <- "lat"
  }

  if (!all(c("long", "lat") %in% names(df))) {
    stop(
      "Each kriging file must contain coordinate columns named ",
      "'long' and 'lat' (or recognizable equivalents)."
    )
  }

  df
}

extract_year_columns <- function(df) {
  year_values <- suppressWarnings(
    as.integer(sub("^X", "", names(df)))
  )
  valid <- !is.na(year_values)

  if (!any(valid)) {
    stop("No yearly columns were detected in a kriging file.")
  }

  output <- names(df)[valid]
  names(output) <- as.character(year_values[valid])
  output
}

read_kriging_file <- function(path) {
  df <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  df <- drop_accidental_index_column(df)
  df <- standardize_coordinate_names(df)
  df$long <- as.numeric(df$long)
  df$lat <- as.numeric(df$lat)
  df
}

tempe_all <- read_kriging_file(mean_file)
variance_reaches <- read_kriging_file(variance_file)

if (nrow(tempe_all) != nrow(variance_reaches)) {
  stop(
    "The mean and prediction-variance files have different row counts: ",
    nrow(tempe_all), " versus ", nrow(variance_reaches), "."
  )
}

coordinate_difference <- max(
  abs(tempe_all$long - variance_reaches$long),
  abs(tempe_all$lat - variance_reaches$lat),
  na.rm = TRUE
)

if (!is.finite(coordinate_difference) || coordinate_difference > 1e-8) {
  stop(
    "The coordinates in the mean and prediction-variance files ",
    "are not aligned row by row."
  )
}

mean_year_map <- extract_year_columns(tempe_all)
variance_year_map <- extract_year_columns(variance_reaches)
common_years <- intersect(names(mean_year_map), names(variance_year_map))
common_years <- as.character(sort(as.integer(common_years)))

if (length(common_years) == 0L) {
  stop("The mean and prediction-variance files have no common years.")
}

mean_year_columns <- unname(mean_year_map[common_years])
variance_year_columns <- unname(variance_year_map[common_years])

variance_values <- as.matrix(
  variance_reaches[, variance_year_columns, drop = FALSE]
)
storage.mode(variance_values) <- "double"

if (any(variance_values < -1e-10, na.rm = TRUE)) {
  stop(
    "Negative values were found in reaches_kriging_grid12x14_variance.csv. ",
    "This file must contain nonnegative prediction variances, apart from negligible numerical error."
  )
}

# ------------------------------------------------------------
# 3. Keep locations in or near China, Taiwan, Hong Kong,
#    and Macao
# ------------------------------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")

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

if (nrow(targets) == 0L) {
  stop("The requested East Asian regions were not found.")
}

china_tw_hk_mo <- st_union(targets)
china_tw_hk_mo <- st_as_sf(
  data.frame(name = "China+Taiwan+HK+MO"),
  geometry = china_tw_hk_mo
)
st_crs(china_tw_hk_mo) <- 4326
china_tw_hk_mo_valid <- st_make_valid(china_tw_hk_mo)

coast_buffer_km <- 80
region_buffer <- china_tw_hk_mo_valid %>%
  st_transform(3857) %>%
  st_buffer(dist = coast_buffer_km * 1000) %>%
  st_transform(4326)

mean_points_sf <- st_as_sf(
  tempe_all,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

keep_location <- lengths(
  st_intersects(mean_points_sf, region_buffer)
) > 0

tempe_all1 <- tempe_all[keep_location, , drop = FALSE]
variance_reaches1 <- variance_reaches[keep_location, , drop = FALSE]

if (nrow(tempe_all1) == 0L) {
  stop("No kriging locations remained after the regional mask.")
}

# ------------------------------------------------------------
# 4. Load and validate LME data
# ------------------------------------------------------------

lme_df <- read.csv(
  lme_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if ("lon" %in% names(lme_df) && !"long" %in% names(lme_df)) {
  names(lme_df)[names(lme_df) == "lon"] <- "long"
}
if ("lat" %in% names(lme_df) && !"lati" %in% names(lme_df)) {
  names(lme_df)[names(lme_df) == "lat"] <- "lati"
}

if (!all(c("long", "lati") %in% names(lme_df))) {
  stop("b0.csv must contain coordinate columns named 'long' and 'lati'.")
}

n_loc <- 121L
n_rep <- 13L
expected_lme_rows <- n_loc * n_rep

if (nrow(lme_df) < expected_lme_rows) {
  stop(
    "Expected at least ", expected_lme_rows,
    " rows in b0.csv, but found ", nrow(lme_df), "."
  )
}

get_block <- function(df, r) {
  first_row <- (r - 1L) * n_loc + 1L
  last_row <- r * n_loc
  df[first_row:last_row, , drop = FALSE]
}

get_xs_all_years <- function(lme_df, s) {
  xs_list <- lapply(
    seq_len(n_rep),
    function(r) {
      block <- get_block(lme_df, r)
      as.numeric(block[s, -c(1, 2), drop = TRUE])
    }
  )

  xs_kelvin <- unlist(xs_list, use.names = FALSE)
  xs_kelvin - 273.15
}

# ------------------------------------------------------------
# 5. Quantile-mapping functions
# ------------------------------------------------------------

build_Fx_inv_local <- function(x_s) {
  x_s <- as.numeric(x_s)
  x_s <- x_s[is.finite(x_s)]

  if (length(x_s) < 2L) {
    stop("Too few finite LME values for a local quantile-mapping fit.")
  }

  x_sd <- sd(x_s)
  if (!is.finite(x_sd) || x_sd <= 0) {
    stop("The local LME sample has zero or invalid variation.")
  }

  bw <- npudistbw(dat = x_s)

  Fx_hat <- function(q) {
    fit <- npudist(
      bws = bw,
      edat = data.frame(x = q)
    )
    as.numeric(fitted(fit))
  }

  Fx_inv <- function(u) {
    u <- pmin(pmax(u, 1e-8), 1 - 1e-8)
    xmin <- min(x_s) - 5 * x_sd
    xmax <- max(x_s) + 5 * x_sd

    sapply(
      u,
      function(uu) {
        uniroot(
          function(q) Fx_hat(q) - uu,
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
  nu <- as.numeric(nu)

  valid <- is.finite(yhat) & is.finite(nu)
  yhat <- yhat[valid]
  nu <- nu[valid]

  if (length(yhat) == 0L) {
    stop(
      "No valid paired kriging means and prediction standard deviations ",
      "were available at a location."
    )
  }

  # nu contains prediction standard deviations obtained by
  # taking square roots of the stored prediction variances.
  nu <- pmax(nu, 1e-8)
  mean(pnorm((y - yhat) / nu))
}

# ------------------------------------------------------------
# 6. Match the 121 LME locations to the kriged REACHES grid
# ------------------------------------------------------------

lme_base <- get_block(lme_df, 1) %>%
  transmute(
    lon = as.numeric(long),
    lat = as.numeric(lati)
  )

reaches_base <- tempe_all1 %>%
  transmute(
    lon = as.numeric(long),
    lat = as.numeric(lat)
  )

make_coordinate_key <- function(lon, lat, digits = 3) {
  paste0(round(lon, digits), "_", round(lat, digits))
}

coord_key_lme <- make_coordinate_key(lme_base$lon, lme_base$lat)
coord_key_reaches <- make_coordinate_key(reaches_base$lon, reaches_base$lat)

if (anyDuplicated(coord_key_reaches)) {
  stop("The filtered REACHES grid contains duplicated coordinate keys.")
}

idx_map <- match(coord_key_lme, coord_key_reaches)

if (anyNA(idx_map)) {
  missing_locations <- lme_base[is.na(idx_map), , drop = FALSE]
  stop(
    "Some LME locations could not be matched to the filtered REACHES grid. ",
    "First unmatched location: ",
    paste(missing_locations[1, ], collapse = ", ")
  )
}

# ------------------------------------------------------------
# 7. Compute g_s(0) for all 121 locations
# ------------------------------------------------------------

# npudistbw uses multistart numerical optimization.
# The fixed seed makes the updated figure reproducible.
set.seed(10)

g0_121 <- vapply(
  seq_len(n_loc),
  function(s) {
    x_s <- get_xs_all_years(lme_df, s)
    Fx_inv_s <- build_Fx_inv_local(x_s)
    reaches_row <- idx_map[s]

    yhat_s <- as.numeric(
      tempe_all1[
        reaches_row,
        mean_year_columns,
        drop = TRUE
      ]
    )

    # The variance file stores kriging prediction variances.
    # FY_hat_1loc() requires prediction standard deviations.
    prediction_variance_s <- as.numeric(
      variance_reaches1[
        reaches_row,
        variance_year_columns,
        drop = TRUE
      ]
    )

    if (any(
      prediction_variance_s < -1e-10,
      na.rm = TRUE
    )) {
      stop(
        "Negative prediction variances were found for LME location ",
        s,
        "."
      )
    }

    nu_s <- sqrt(
      pmax(
        prediction_variance_s,
        0
      )
    )

    u0 <- FY_hat_1loc(
      y = 0,
      yhat = yhat_s,
      nu = nu_s
    )

    Fx_inv_s(u0)
  },
  numeric(1)
)

plot_df <- lme_base %>%
  mutate(temp0_c = g0_121)

# ------------------------------------------------------------
# 8. Smooth the 121 mapped values to a 0.25-degree grid
# ------------------------------------------------------------

smooth_fit <- gam(
  temp0_c ~ s(lon, lat, bs = "tp", fx = TRUE),
  data = plot_df,
  method = "REML"
)

lon_seq <- seq(98, 124.5, by = 0.25)
lat_seq <- seq(18, 42.5, by = 0.25)

smooth_grid <- expand.grid(
  lon = lon_seq,
  lat = lat_seq
)

smooth_grid$temp0_c <- as.numeric(
  predict(smooth_fit, newdata = smooth_grid)
)

# Keep only land locations in the target regions, using the
# original 1-km coastal buffer from the smoothing code.
smooth_buffer_km <- 1
smooth_region_buffer <- china_tw_hk_mo_valid %>%
  st_transform(3857) %>%
  st_buffer(dist = smooth_buffer_km * 1000) %>%
  st_transform(4326)

smooth_points_sf <- st_as_sf(
  smooth_grid,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

keep_smooth_location <- lengths(
  st_intersects(smooth_points_sf, smooth_region_buffer)
) > 0

smooth_grid <- smooth_grid[
  keep_smooth_location,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 9. Save Figure 7(e) and its underlying values
# ------------------------------------------------------------

p_figure7e <- ggplot(
  smooth_grid,
  aes(lon, lat)
) +
  geom_point(
    aes(colour = temp0_c),
    cex = 4.5,
    shape = 15
  ) +
  coord_map(
    xlim = c(98, 124.5),
    ylim = c(18, 42.5)
  ) +
  scale_color_gradientn(
    colors = c("blue", "cyan", "green", "yellow", "red"),
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

figure7e_output_dir <- here::here("Output", "Figure7-8")
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

write.csv(
  plot_df,
  file.path(
    figure7e_output_dir,
    "Figure7e_qm_locations.csv"
  ),
  row.names = FALSE
)

write.csv(
  smooth_grid,
  file.path(
    figure7e_output_dir,
    "Figure7e_smoothed_grid.csv"
  ),
  row.names = FALSE
)

message(
  "Figure 7(e) saved to: ",
  figure7e_output_file
)
