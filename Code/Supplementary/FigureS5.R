here::i_am("Code/Supplementary/FigureS5.R")

# ============================================================
# Generate 12 separate panels for Supplementary Figure S5:
#
#   REACHES clusters:
#     FigureS5_reaches_cluster_1.png
#     FigureS5_reaches_cluster_2.png
#     FigureS5_reaches_cluster_3.png
#     FigureS5_reaches_cluster_4.png
#
#   LME clusters:
#     FigureS5_lme_cluster_1.png
#     FigureS5_lme_cluster_2.png
#     FigureS5_lme_cluster_3.png
#     FigureS5_lme_cluster_4.png
#
#   Posterior clusters:
#     FigureS5_posterior_cluster_1.png
#     FigureS5_posterior_cluster_2.png
#     FigureS5_posterior_cluster_3.png
#     FigureS5_posterior_cluster_4.png
#
# Required inputs:
#   1. Figure 6 cluster assignments
#   2. REACHES kriged data on the LME 53x49 grid
#   3. LME ensemble-mean annual data on the same grid
#   4. Assimilated posterior annual mean on the same grid
#
# This script assumes the current analysis uses 4 clusters.
# ============================================================

library(here)
library(dplyr)
library(fda)

# ------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) {
    return(NA_character_)
  }
  hit[1]
}

read_grid_time_series <- function(file_path) {
  if (is.na(file_path) || !file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  dat <- read.csv(file_path, check.names = FALSE)

  # latitude column name
  lat_col <- if ("lati" %in% names(dat)) {
    "lati"
  } else if ("lat" %in% names(dat)) {
    "lat"
  } else {
    stop("Could not find a latitude column named 'lati' or 'lat' in: ", file_path)
  }

  if (!"long" %in% names(dat)) {
    stop("Could not find longitude column 'long' in: ", file_path)
  }

  year_cols <- grep("^(x)?[0-9]{4}$", names(dat), value = TRUE)

  if (length(year_cols) == 0L) {
    stop("No yearly columns found in: ", file_path)
  }

  year_num <- as.numeric(sub("^x", "", year_cols))
  keep <- year_num >= 1368 & year_num <= 1911

  year_cols <- year_cols[keep]
  year_num <- year_num[keep]

  if (length(year_cols) == 0L) {
    stop("No year columns between 1368 and 1911 were found in: ", file_path)
  }

  out <- dat %>%
    dplyr::select(long, !!lat_col, all_of(year_cols))

  names(out)[2] <- "lati"

  list(
    data = out,
    year_cols = year_cols,
    years = year_num
  )
}

prepare_cluster_assignment <- function(file_path) {
  if (is.na(file_path) || !file.exists(file_path)) {
    stop("Cluster-assignment file not found: ", file_path)
  }

  dat <- read.csv(file_path, check.names = FALSE)

  lat_col <- if ("lati" %in% names(dat)) {
    "lati"
  } else if ("lat" %in% names(dat)) {
    "lat"
  } else {
    stop("Could not find 'lati' or 'lat' in cluster-assignment file.")
  }

  cluster_col <- NULL
  candidate_cluster_cols <- c(
    "cluster",
    "figure6_cluster",
    "cluster_id",
    "Cluster"
  )

  for (nm in candidate_cluster_cols) {
    if (nm %in% names(dat)) {
      cluster_col <- nm
      break
    }
  }

  if (is.null(cluster_col)) {
    stop(
      "Could not find a cluster column. Tried: ",
      paste(candidate_cluster_cols, collapse = ", ")
    )
  }

  out <- dat %>%
    dplyr::select(long, !!lat_col, !!cluster_col)

  names(out) <- c("long", "lati", "cluster")
  out$cluster <- as.integer(out$cluster)

  out
}

add_grid_key <- function(dat) {
  dat %>%
    mutate(
      grid_key = paste0(
        sprintf("%.6f", as.numeric(long)),
        "_",
        sprintf("%.6f", as.numeric(lati))
      )
    )
}

row_center_matrix <- function(mat) {
  mat <- as.matrix(mat)
  mat - rowMeans(mat, na.rm = TRUE)
}

get_ylim_with_padding <- function(mat, pad_ratio = 0.05) {
  rng <- range(mat, finite = TRUE, na.rm = TRUE)
  span <- rng[2] - rng[1]
  if (!is.finite(span) || span <= 0) {
    span <- 1
  }
  c(rng[1] - pad_ratio * span, rng[2] + pad_ratio * span)
}

draw_cluster_fbplot <- function(
    mat,
    years,
    cluster_id,
    title_text,
    ylim_use = NULL,
    ylab_text = ""
) {
  idx <- which(cluster_vector == cluster_id)

  if (length(idx) == 0L) {
    plot.new()
    title(main = paste0(title_text, " (no locations)"))
    return(invisible(NULL))
  }

  mat_use <- mat[idx, , drop = FALSE]

  if (is.null(ylim_use)) {
    ylim_use <- get_ylim_with_padding(mat_use)
  }

  fbplot(
    t(mat_use),
    x = years,
    xlim = range(years),
    ylim = ylim_use,
    axes = FALSE,
    main = title_text,
    xlab = "",
    ylab = ""
  )

  axis(
    side = 1,
    at = seq(1400, 1900, by = 100),
    labels = seq(1400, 1900, by = 100),
    cex.axis = 0.8
  )

  axis(
    side = 2,
    las = 1,
    cex.axis = 0.8
  )

  mtext("Year", side = 1, line = 2.2, cex = 0.9)
  mtext(ylab_text, side = 2, line = 2.5, cex = 0.9)

  box()
  invisible(NULL)
}

save_cluster_panel <- function(
    mat,
    years,
    cluster_id,
    title_text,
    output_file,
    ylim_use = NULL,
    ylab_text = "",
    width = 6,
    height = 4,
    res = 300
) {
  png(
    filename = output_file,
    width = width,
    height = height,
    units = "in",
    res = res
  )

  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  par(
    mfrow = c(1, 1),
    mar = c(3.5, 4, 2, 1),
    mgp = c(2, 0.7, 0),
    tcl = -0.2
  )

  draw_cluster_fbplot(
    mat = mat,
    years = years,
    cluster_id = cluster_id,
    title_text = title_text,
    ylim_use = ylim_use,
    ylab_text = ylab_text
  )

  invisible(output_file)
}

# ------------------------------------------------------------
# 2. Locate input files
# ------------------------------------------------------------

cluster_file <- first_existing_file(c(
  here("Output", "Figure6", "Figure6_cluster_assignments.csv"),
  here("Output", "Intermediate", "FigureS5", "FigureS5_cluster_assignments.csv")
))

reaches_file <- first_existing_file(c(
  here("Data", "reaches_kriging_grid53x49_mean.csv"),
  here("Data", "reaches_kriging_lme_grid_mean.csv"),
  here("Output", "Intermediate", "REACHES", "reaches_kriging_lme_grid_mean.csv")
))

lme_file <- first_existing_file(c(
  here("Data", "LME data", "lme_ensemble_mean_1368_1911.csv"),
  here("Output", "Intermediate", "LME", "lme_ensemble_mean_1368_1911.csv")
))

posterior_file <- first_existing_file(c(
  here("Data", "Valid", "assimilated_posterior_lme_grid_mean.csv"),
  here("Output", "Intermediate", "Assimilation", "assimilated_posterior_lme_grid_mean.csv")
))

if (is.na(cluster_file)) {
  stop("Could not find Figure 6 cluster assignments csv.")
}
if (is.na(reaches_file)) {
  stop("Could not find REACHES kriged grid csv.")
}
if (is.na(lme_file)) {
  stop("Could not find LME ensemble-mean csv.")
}
if (is.na(posterior_file)) {
  stop("Could not find assimilated posterior mean csv.")
}

# ------------------------------------------------------------
# 3. Read inputs
# ------------------------------------------------------------

cluster_df <- prepare_cluster_assignment(cluster_file)

reaches_obj <- read_grid_time_series(reaches_file)
lme_obj <- read_grid_time_series(lme_file)
posterior_obj <- read_grid_time_series(posterior_file)

common_years <- Reduce(
  intersect,
  list(reaches_obj$years, lme_obj$years, posterior_obj$years)
)

common_years <- sort(common_years)

if (length(common_years) == 0L) {
  stop("No common analysis years across REACHES, LME, and posterior.")
}

reaches_year_cols <- reaches_obj$year_cols[reaches_obj$years %in% common_years]
lme_year_cols <- lme_obj$year_cols[lme_obj$years %in% common_years]
posterior_year_cols <- posterior_obj$year_cols[posterior_obj$years %in% common_years]

cluster_df <- add_grid_key(cluster_df)
reaches_df <- add_grid_key(reaches_obj$data)
lme_df <- add_grid_key(lme_obj$data)
posterior_df <- add_grid_key(posterior_obj$data)

merged_df <- cluster_df %>%
  dplyr::select(grid_key, long, lati, cluster) %>%
  inner_join(
    reaches_df %>% dplyr::select(grid_key, all_of(reaches_year_cols)),
    by = "grid_key"
  ) %>%
  inner_join(
    lme_df %>% dplyr::select(grid_key, all_of(lme_year_cols)),
    by = "grid_key"
  ) %>%
  inner_join(
    posterior_df %>% dplyr::select(grid_key, all_of(posterior_year_cols)),
    by = "grid_key"
  )

if (nrow(merged_df) == 0L) {
  stop("No common grid cells remain after merging all inputs.")
}

cluster_vector <- merged_df$cluster
n_clusters <- length(sort(unique(cluster_vector)))

if (n_clusters != 4L) {
  warning(
    "The current merged data contain ", n_clusters,
    " clusters, not 4. The script will still run using the observed clusters."
  )
}

cluster_ids <- sort(unique(cluster_vector))

# Rebuild separate matrices from merged data
reaches_mat <- as.matrix(merged_df[, reaches_year_cols, drop = FALSE])
lme_mat <- as.matrix(merged_df[, lme_year_cols, drop = FALSE])
posterior_mat <- as.matrix(merged_df[, posterior_year_cols, drop = FALSE])

storage.mode(reaches_mat) <- "double"
storage.mode(lme_mat) <- "double"
storage.mode(posterior_mat) <- "double"

# Keep REACHES as-is.
# Center LME and posterior by grid cell, matching the earlier logic.
lme_mat_centered <- row_center_matrix(lme_mat)
posterior_mat_centered <- row_center_matrix(posterior_mat)

# Global y-limits within each family
ylim_reaches <- get_ylim_with_padding(reaches_mat)
ylim_lme <- get_ylim_with_padding(lme_mat_centered)
ylim_posterior <- get_ylim_with_padding(posterior_mat_centered)

# ------------------------------------------------------------
# 4. Save 12 separate panels
# ------------------------------------------------------------

output_dir <- here("Output", "Supplementary", "FigureS5_panels")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

saved_files <- character(0)

for (k in cluster_ids) {
  out_file <- file.path(
    output_dir,
    paste0("FigureS5_reaches_cluster_", k, ".png")
  )

  save_cluster_panel(
    mat = reaches_mat,
    years = common_years,
    cluster_id = k,
    title_text = paste0("REACHES cluster ", k),
    output_file = out_file,
    ylim_use = ylim_reaches,
    ylab_text = "Index"
  )

  saved_files <- c(saved_files, out_file)
}

for (k in cluster_ids) {
  out_file <- file.path(
    output_dir,
    paste0("FigureS5_lme_cluster_", k, ".png")
  )

  save_cluster_panel(
    mat = lme_mat_centered,
    years = common_years,
    cluster_id = k,
    title_text = paste0("LME cluster ", k),
    output_file = out_file,
    ylim_use = ylim_lme,
    ylab_text = "Centered temperature"
  )

  saved_files <- c(saved_files, out_file)
}

for (k in cluster_ids) {
  out_file <- file.path(
    output_dir,
    paste0("FigureS5_posterior_cluster_", k, ".png")
  )

  save_cluster_panel(
    mat = posterior_mat_centered,
    years = common_years,
    cluster_id = k,
    title_text = paste0("Posterior cluster ", k),
    output_file = out_file,
    ylim_use = ylim_posterior,
    ylab_text = "Centered temperature"
  )

  saved_files <- c(saved_files, out_file)
}

message("Saved the following Figure S5 panels:")
message(paste(saved_files, collapse = "\n"))
