here::i_am("Code/Figure6.R")

# ============================================================
# REACHES-only analysis
#
# Main outputs:
#   1. p_reach_map              : map of the five REACHES clusters
#   2. plot_reaches_clusters()  : five REACHES functional boxplots
#
# No posterior or LME code is included.
# ============================================================

library(readxl)
library(dplyr)
library(sf)
library(rnaturalearth)
library(fda)
library(mclust)
library(ggplot2)
library(maps)

# ------------------------------------------------------------
# 1. Read REACHES kriged data
# ------------------------------------------------------------
tempe_all_v4 <- read.csv(here::here("Data", "tempe_figure6.csv"))

# Keep grid cells within 150 km of China or Taiwan
pts <- tempe_all_v4 %>%
  mutate(row_id = row_number()) %>%
  st_as_sf(
    coords = c("long", "lati"),
    crs = 4326,
    remove = FALSE
  )

world <- ne_countries(scale = "medium", returnclass = "sf")
easia <- world %>%
  filter(admin %in% c("China", "Taiwan"))

buffer_km <- 150

easia_m <- st_transform(easia, 3857)
pts_m <- st_transform(pts, 3857)
easia_buffer_m <- st_buffer(easia_m, dist = buffer_km * 1000)

inside_easia_km <- lengths(st_intersects(pts_m, easia_buffer_m)) > 0
tempe_all_v4 <- tempe_all_v4[inside_easia_km, ]

# ------------------------------------------------------------
# 2. Read original REACHES observations
# ------------------------------------------------------------
temperature <- read_excel(
  here::here("Data", "temperature index value.v1.xlsx"),
  col_types = c(
    "skip", "skip", "numeric", "numeric", "skip", "skip", "skip",
    "skip", "skip", "numeric", "numeric", "skip", "skip"
  )
)

colnames(temperature) <- c("level", "year", "long", "lat")

# ------------------------------------------------------------
# 3. Restrict analysis to 1368--1911
# ------------------------------------------------------------
use_years <- 1368:1911
target_cols <- paste0("x", use_years)

reach_cols <- grep("^x[0-9]+$", names(tempe_all_v4), value = TRUE)
reach_cols_use <- intersect(target_cols, reach_cols)
time_grid <- as.numeric(sub("^x", "", reach_cols_use))


# ------------------------------------------------------------
# 4. Keep only grid cells associated with at least one original
#    REACHES observation during 1368--1911
# ------------------------------------------------------------
obs_sf <- temperature %>%
  filter(
    year >= 1368,
    year <= 1911,
    is.finite(level),
    is.finite(long),
    is.finite(lat)
  ) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326, remove = FALSE)

grid_sf <- tempe_all_v4 %>%
  mutate(grid_id = row_number()) %>%
  st_as_sf(coords = c("long", "lati"), crs = 4326, remove = FALSE)

obs_m <- st_transform(obs_sf, 3857)
grid_m <- st_transform(grid_sf, 3857)

nearest_grid_id <- st_nearest_feature(obs_m, grid_m)
used_grid_id <- sort(unique(nearest_grid_id))

reaches_grid <- tempe_all_v4 %>%
  mutate(grid_id = row_number()) %>%
  filter(grid_id %in% used_grid_id) %>%
  dplyr::select(long, lati, all_of(reach_cols_use))

# ------------------------------------------------------------
# 5. REACHES FDA + FPCA + Mclust
# ------------------------------------------------------------
Y_R <- reaches_grid %>%
  dplyr::select(all_of(reach_cols_use)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

basis_R <- create.bspline.basis(
  rangeval = range(time_grid),
  nbasis = 20
)

fd_R <- Data2fd(
  argvals = time_grid,
  y = t(Y_R),
  basisobj = basis_R
)

pca_R <- pca.fd(fd_R, nharm = 5)
scores_R <- pca_R$scores

# Keep the same model-selection setup as the original code.
# In the original analysis, BIC selects five clusters.
mc_R <- Mclust(scores_R, G = 1:8)
cl_R <- mc_R$classification


# Relabel the five clusters to match the geographic ordering
# in the original analysis:
# old 2, 4, 5, 1, 3 -> new 1, 2, 3, 4, 5
relabel_R <- c(
  "2" = 1,
  "4" = 2,
  "5" = 3,
  "1" = 4,
  "3" = 5
)

cl_R_geo <- unname(relabel_R[as.character(cl_R)])


# Evaluate the fitted REACHES functions on the yearly grid.
# Rows = grid cells; columns = years.
eval_R <- t(eval.fd(time_grid, fd_R))

# Convenient list of grid-cell indices for clusters 1--5
reaches_cluster_indices <- lapply(1:5, function(k) which(cl_R_geo == k))
names(reaches_cluster_indices) <- paste0("cluster_", 1:5)

# ------------------------------------------------------------
# 6. Five REACHES-only functional boxplots
# ------------------------------------------------------------
year_ticks <- seq(1400, 1900, by = 100)
ylim_R <- c(-1.7, 0.3)
ytick_R <- c(-1.5, -1, -0.5, 0)

draw_reaches_cluster <- function(k) {
  if (!k %in% 1:5) {
    stop("k must be one of 1, 2, 3, 4, or 5.")
  }

  idx <- reaches_cluster_indices[[k]]

  if (length(idx) == 0L) {
    plot.new()
    title(main = paste0("REACHES Cluster ", k, " (no grid cells)"))
    return(invisible(NULL))
  }

  fbplot(
    t(eval_R[idx, , drop = FALSE]),
    x = time_grid,
    xlim = range(time_grid),
    ylim = ylim_R,
    axes = FALSE,
    main = paste0("REACHES Cluster ", k),
    xlab = "",
    ylab = ""
  )

  axis(
    side = 1,
    at = year_ticks,
    labels = year_ticks,
    cex.axis = 0.8
  )

  axis(
    side = 2,
    at = ytick_R,
    labels = ytick_R,
    las = 1,
    cex.axis = 0.8
  )

  box()
  invisible(NULL)
}

# Draw all five clusters in one 5 x 1 figure.
plot_reaches_cluster <- function(k) {
  stopifnot(k %in% 1:5)

  idx <- which(cl_R_geo == k)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(
    mfrow = c(1, 1),
    mar = c(3.5, 4, 2, 1),
    mgp = c(2, 0.7, 0),
    tcl = -0.2
  )

  draw_fb(
    eval_mat = eval_R,
    idx = idx,
    ylim_use = ylim_R,
    ytick_use = ytick_R
  )

  title(main = paste("REACHES Cluster", k))
}

plot_reaches_clusters <- function() {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(
    mfrow = c(5, 1),
    mar = c(2.5, 4, 1.5, 1),
    mgp = c(1.8, 0.5, 0),
    tcl = -0.2
  )

  for (k in 1:5) {
    idx <- which(cl_R_geo == k)
    draw_fb(
      eval_mat = eval_R,
      idx = idx,
      ylim_use = ylim_R,
      ytick_use = ytick_R
    )
    title(main = paste("REACHES Cluster", k))
  }
}

# Save all five clusters to one image.
save_reaches_clusters <- function(
    filename = "reaches_only_5_clusters.jpg",
    width = 6,
    height = 12,
    res = 300) {

  jpeg(
    filename = filename,
    width = width,
    height = height,
    res = res,
    units = "in"
  )

  on.exit(dev.off(), add = TRUE)
  plot_reaches_clusters()
  invisible(filename)
}

# Save the five clusters as five separate image files.
save_reaches_clusters_separately <- function(
    output_dir = "reaches_only_clusters",
    width = 6,
    height = 4,
    res = 300) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  for (k in 1:5) {
    filename <- file.path(
      output_dir,
      paste0("reaches_cluster_", k, ".jpg")
    )

    jpeg(
      filename = filename,
      width = width,
      height = height,
      res = res,
      units = "in"
    )

    draw_reaches_cluster(k)
    dev.off()
  }

  invisible(output_dir)
}

# ------------------------------------------------------------
# 7. REACHES-only five-cluster map: p_reach_map
# ------------------------------------------------------------
get_grid_size <- function(df) {
  dx <- min(diff(sort(unique(df$long))), na.rm = TRUE)
  dy <- min(diff(sort(unique(df$lat))), na.rm = TRUE)
  list(dx = dx, dy = dy)
}

df_map_reach <- data.frame(
  long = reaches_grid$long,
  lat = reaches_grid$lati,
  cluster = factor(cl_R_geo, levels = 1:5)
)

gs_reach <- get_grid_size(df_map_reach)

p_reach_map <- ggplot(
  df_map_reach,
  aes(x = long, y = lat, fill = cluster)
) +
  geom_tile(
    width = gs_reach$dx * 0.98,
    height = gs_reach$dy * 0.98
  ) +
  borders(
    "world",
    xlim = c(76, 132),
    ylim = c(18, 52),
    fill = NA,
    colour = "grey40"
  ) +
  coord_map(
    xlim = c(98, 130.5),
    ylim = c(18, 42.5)
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Cluster"
  ) +
  theme_minimal()

# ------------------------------------------------------------
# 8. Commands to display or save the requested outputs
# ------------------------------------------------------------

# Display the map:
print(p_reach_map)

plot_reaches_cluster <- function(k) {
  draw_reaches_cluster(k)
}

# Display the five REACHES functional boxplots.
plot_reaches_cluster(1)
plot_reaches_cluster(2)
plot_reaches_cluster(3)
plot_reaches_cluster(4)
plot_reaches_cluster(5)
