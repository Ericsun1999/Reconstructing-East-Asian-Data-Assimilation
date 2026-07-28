here::i_am("Code/Figure6.R")

# ============================================================
# REACHES-only functional clustering for Figure 6
#
# The REACHES temperature-index field must be evaluated at the
# native LME grid locations. This is not the same as using every
# point of the complete 53 x 49 REACHES prediction grid.
#
# Supported input modes:
#   1. "generated":
#      Output/Intermediate/REACHES/
#        reaches_kriging_lme_grid_mean.csv
#
#   2. "precomputed":
#      Data/REACHES/precomputed/
#        reaches_kriging_lme_grid_mean.csv
#
#   3. "auto" (default):
#      Use the generated file when available; otherwise use the
#      precomputed repository copy.
#
# Outputs:
#   Output/Figure6/Figure6_cluster_map.jpg
#   Output/Figure6/Figure6_cluster_1.jpg
#   ...
#   Output/Figure6/Figure6_cluster_K.jpg
#   Output/Figure6/Figure6_cluster_assignments.csv
#   Output/Figure6/Figure6_mclust_bic.csv
#   Output/Figure6/Figure6_cluster_selection.csv
#
# K is selected by BIC by default. Set cluster_selection_mode
# to "fixed" only when a prespecified number of clusters is
# scientifically required.
# ============================================================

library(dplyr)
library(sf)
library(rnaturalearth)
library(fda)
library(mclust)
library(ggplot2)
library(maps)
library(readr)
library(readxl)

# ------------------------------------------------------------
# 1. Select the REACHES-at-LME-grid input
# ------------------------------------------------------------

input_mode <- "auto"

cluster_selection_mode <- "bic"
fixed_number_of_clusters <- 5L
candidate_numbers_of_clusters <- 1:8

allowed_cluster_selection_modes <- c(
  "bic",
  "fixed"
)

if (!cluster_selection_mode %in%
    allowed_cluster_selection_modes) {
  stop(
    "cluster_selection_mode must be one of: ",
    paste(
      allowed_cluster_selection_modes,
      collapse = ", "
    )
  )
}

allowed_input_modes <- c(
  "auto",
  "generated",
  "precomputed"
)

if (!input_mode %in% allowed_input_modes) {
  stop(
    "input_mode must be one of: ",
    paste(
      allowed_input_modes,
      collapse = ", "
    )
  )
}

input_files <- c(
  generated = here::here(
    "Output",
    "Intermediate",
    "REACHES",
    "reaches_kriging_lme_grid_mean.csv"
  ),
  precomputed = here::here(
    "Data",
    "reaches_kriging_lme_grid_mean.csv"
  )
)

select_input_file <- function(
    input_mode,
    input_files) {

  generated_file <- unname(
    input_files["generated"]
  )

  precomputed_file <- unname(
    input_files["precomputed"]
  )

  generated_available <- file.exists(
    generated_file
  )

  precomputed_available <- file.exists(
    precomputed_file
  )

  if (input_mode == "generated") {
    if (!generated_available) {
      stop(
        "The generated Figure 6 input was not found:\n  ",
        generated_file,
        "\nRun Code/Get_tempe_all_data.R first."
      )
    }

    return(
      generated_file
    )
  }

  if (input_mode == "precomputed") {
    if (!precomputed_available) {
      stop(
        "The precomputed Figure 6 input was not found:\n  ",
        precomputed_file
      )
    }

    return(
      precomputed_file
    )
  }

  if (generated_available) {
    return(
      generated_file
    )
  }

  if (precomputed_available) {
    return(
      precomputed_file
    )
  }

  stop(
    "Neither Figure 6 input file was found.\n",
    "Generated location:\n  ",
    generated_file,
    "\nPrecomputed location:\n  ",
    precomputed_file,
    "\nRun Code/Get_tempe_all_data.R or add the ",
    "precomputed file under Data/REACHES/precomputed/."
  )
}

figure6_input_file <- select_input_file(
  input_mode = input_mode,
  input_files = input_files
)

message(
  "Figure 6 input selected: ",
  figure6_input_file
)

# ------------------------------------------------------------
# 2. Read and standardize the kriged REACHES field
# ------------------------------------------------------------

tempe_all_v4 <- readr::read_csv(
  figure6_input_file,
  show_col_types = FALSE
)

if (
  !"lati" %in%
    names(tempe_all_v4) &&
    "lat" %in%
      names(tempe_all_v4)
) {
  tempe_all_v4 <- tempe_all_v4 %>%
    rename(
      lati = lat
    )
}

required_coordinate_columns <- c(
  "long",
  "lati"
)

missing_coordinate_columns <- setdiff(
  required_coordinate_columns,
  names(
    tempe_all_v4
  )
)

if (length(missing_coordinate_columns) > 0L) {
  stop(
    "The Figure 6 input is missing coordinate columns: ",
    paste(
      missing_coordinate_columns,
      collapse = ", "
    )
  )
}

# Accept either x1368 or 1368 as an input column name and
# standardize all annual columns to xYYYY.
raw_year_columns <- grep(
  "^x?[0-9]{4}$",
  names(
    tempe_all_v4
  ),
  value = TRUE
)

if (length(raw_year_columns) == 0L) {
  stop(
    "No annual REACHES columns were found in the Figure 6 input."
  )
}

raw_year_values <- as.integer(
  sub(
    "^x",
    "",
    raw_year_columns
  )
)

names(
  tempe_all_v4
)[
  match(
    raw_year_columns,
    names(
      tempe_all_v4
    )
  )
] <- paste0(
  "x",
  raw_year_values
)

tempe_all_v4 <- tempe_all_v4 %>%
  mutate(
    long = as.numeric(long),
    lati = as.numeric(lati)
  )

if (
  any(!is.finite(tempe_all_v4$long)) ||
    any(!is.finite(tempe_all_v4$lati))
) {
  stop(
    "The Figure 6 input contains invalid coordinates."
  )
}

if (anyDuplicated(
  tempe_all_v4[
    ,
    c(
      "long",
      "lati"
    )
  ]
)) {
  stop(
    "The Figure 6 input contains duplicated LME-grid ",
    "coordinates."
  )
}

# ------------------------------------------------------------
# 3. Keep LME grid cells within 150 km of China or Taiwan
# ------------------------------------------------------------

sf::sf_use_s2(
  TRUE
)

pts <- tempe_all_v4 %>%
  mutate(
    row_id = row_number()
  ) %>%
  st_as_sf(
    coords = c(
      "long",
      "lati"
    ),
    crs = 4326,
    remove = FALSE
  )

world <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
)

easia <- world %>%
  filter(
    admin %in%
      c(
        "China",
        "Taiwan"
      )
  ) %>%
  st_make_valid() %>%
  st_union()

buffer_distance_m <- 150000

inside_easia_150km <- lengths(
  st_is_within_distance(
    pts,
    easia,
    dist = buffer_distance_m
  )
) > 0L

tempe_all_v4 <- tempe_all_v4[
  inside_easia_150km,
  ,
  drop = FALSE
]

if (nrow(tempe_all_v4) == 0L) {
  stop(
    "No LME grid locations remained after the 150-km ",
    "geographic restriction."
  )
}

# ------------------------------------------------------------
# 4. Read original REACHES observations
# ------------------------------------------------------------

calibration_file <- here::here(
  "Output",
  "Intermediate",
  "calibration_parameters.rds"
)

raw_reaches_file <- here::here(
  "Data",
  "temperature index value.v1.xlsx"
)

if (file.exists(calibration_file)) {

  calibration_results <- readRDS(
    calibration_file
  )

  if (!"temp2" %in% names(calibration_results)) {
    stop(
      "The calibration archive does not contain temp2."
    )
  }

  temperature <- calibration_results$temp2 %>%
    as.data.frame() %>%
    dplyr::select(
      level,
      year,
      long,
      lat
    )

} else {

  if (!file.exists(raw_reaches_file)) {
    stop(
      "Neither the calibration archive nor the raw REACHES ",
      "workbook was found."
    )
  }

  temperature <- readxl::read_excel(
    raw_reaches_file,
    col_types = c(
      "skip",
      "skip",
      "numeric",
      "numeric",
      "skip",
      "skip",
      "skip",
      "skip",
      "skip",
      "numeric",
      "numeric",
      "skip",
      "skip"
    )
  )

  colnames(
    temperature
  ) <- c(
    "level",
    "year",
    "long",
    "lat"
  )
}

temperature <- temperature %>%
  transmute(
    level = as.numeric(level),
    year = as.integer(year),
    long = as.numeric(long),
    lat = as.numeric(lat)
  )

# ------------------------------------------------------------
# 5. Restrict the functional analysis to 1368--1911
# ------------------------------------------------------------

use_years <- 1368:1911

target_columns <- paste0(
  "x",
  use_years
)

available_reaches_columns <- grep(
  "^x[0-9]{4}$",
  names(
    tempe_all_v4
  ),
  value = TRUE
)

reaches_columns_used <- intersect(
  target_columns,
  available_reaches_columns
)

time_grid <- as.integer(
  sub(
    "^x",
    "",
    reaches_columns_used
  )
)

if (length(time_grid) < 2L) {
  stop(
    "At least two REACHES event years are required for FDA."
  )
}

if (is.unsorted(time_grid)) {
  stop(
    "The REACHES year columns are not chronologically ordered."
  )
}

# ------------------------------------------------------------
# 6. Keep LME-grid cells associated with an original REACHES
#    observation during 1368--1911
# ------------------------------------------------------------

observation_sf <- temperature %>%
  filter(
    year >= 1368L,
    year <= 1911L,
    is.finite(level),
    is.finite(long),
    is.finite(lat)
  ) %>%
  st_as_sf(
    coords = c(
      "long",
      "lat"
    ),
    crs = 4326,
    remove = FALSE
  )

grid_sf <- tempe_all_v4 %>%
  mutate(
    grid_id = row_number()
  ) %>%
  st_as_sf(
    coords = c(
      "long",
      "lati"
    ),
    crs = 4326,
    remove = FALSE
  )

nearest_grid_id <- st_nearest_feature(
  observation_sf,
  grid_sf
)

used_grid_id <- sort(
  unique(
    nearest_grid_id
  )
)

reaches_grid <- tempe_all_v4 %>%
  mutate(
    grid_id = row_number()
  ) %>%
  filter(
    grid_id %in%
      used_grid_id
  ) %>%
  dplyr::select(
    any_of(
      "location_id"
    ),
    long,
    lati,
    all_of(
      reaches_columns_used
    )
  )

minimum_required_locations <- max(
  5L,
  if (
    cluster_selection_mode == "fixed"
  ) {
    as.integer(
      fixed_number_of_clusters
    )
  } else {
    1L
  }
)

if (
  nrow(reaches_grid) <
    minimum_required_locations
) {
  stop(
    "Too few LME-grid locations remained for functional ",
    "clustering."
  )
}

# ------------------------------------------------------------
# 7. REACHES FDA, FPCA, and Gaussian-mixture clustering
# ------------------------------------------------------------

Y_R <- reaches_grid %>%
  dplyr::select(
    all_of(
      reaches_columns_used
    )
  ) %>%
  mutate(
    across(
      everything(),
      as.numeric
    )
  ) %>%
  as.matrix()

if (anyNA(Y_R)) {
  stop(
    "The Figure 6 functional data contain missing values."
  )
}

basis_R <- fda::create.bspline.basis(
  rangeval = range(
    time_grid
  ),
  nbasis = 20
)

fd_R <- fda::Data2fd(
  argvals = time_grid,
  y = t(
    Y_R
  ),
  basisobj = basis_R
)

pca_R <- fda::pca.fd(
  fd_R,
  nharm = 5
)

scores_R <- pca_R$scores

# First fit the complete candidate set so the BIC-selected
# solution and the full BIC table are always retained.
mc_R_bic <- mclust::Mclust(
  scores_R,
  G = candidate_numbers_of_clusters
)

bic_selected_number_of_clusters <- as.integer(
  mc_R_bic$G
)

if (cluster_selection_mode == "bic") {

  mc_R <- mc_R_bic

} else {

  if (
    length(fixed_number_of_clusters) != 1L ||
      !is.finite(fixed_number_of_clusters) ||
      fixed_number_of_clusters < 1L
  ) {
    stop(
      "fixed_number_of_clusters must be one positive integer."
    )
  }

  mc_R <- mclust::Mclust(
    scores_R,
    G = as.integer(
      fixed_number_of_clusters
    )
  )
}

number_of_clusters <- as.integer(
  mc_R$G
)

cl_R <- as.integer(
  mc_R$classification
)

# Evaluate the fitted functions before relabeling. Mixture
# component labels are arbitrary, so assign reproducible Figure 6
# labels by ordering components from the lowest to the highest
# overall fitted REACHES temperature index.
eval_R <- t(
  fda::eval.fd(
    time_grid,
    fd_R
  )
)

overall_curve_mean <- rowMeans(
  eval_R
)

component_mean_index <- tapply(
  overall_curve_mean,
  cl_R,
  mean
)

component_order <- as.integer(
  names(
    sort(
      component_mean_index,
      decreasing = FALSE
    )
  )
)

relabel_R <- setNames(
  seq_len(
    number_of_clusters
  ),
  component_order
)

cl_R_geo <- unname(
  relabel_R[
    as.character(
      cl_R
    )
  ]
)

if (
  anyNA(cl_R_geo) ||
    length(
      unique(
        cl_R_geo
      )
    ) != number_of_clusters
) {
  stop(
    "Deterministic cluster relabeling failed."
  )
}

reaches_cluster_indices <- lapply(
  seq_len(
    number_of_clusters
  ),
  function(k) {
    which(
      cl_R_geo == k
    )
  }
)

names(
  reaches_cluster_indices
) <- paste0(
  "cluster_",
  seq_len(
    number_of_clusters
  )
)

message(
  "BIC selected ",
  bic_selected_number_of_clusters,
  " clusters. Figure 6 is using ",
  number_of_clusters,
  " clusters under cluster_selection_mode = \"",
  cluster_selection_mode,
  "\"."
)

# ------------------------------------------------------------
# 8. Draw the functional boxplots
# ------------------------------------------------------------

year_ticks <- seq(
  1400,
  1900,
  by = 100
)

ylim_R <- c(
  -1.7,
  0.3
)

ytick_R <- c(
  -1.5,
  -1,
  -0.5,
  0
)

draw_reaches_cluster <- function(k) {

  if (!k %in%
      seq_len(
        number_of_clusters
      )) {
    stop(
      "k must be an integer from 1 to ",
      number_of_clusters,
      "."
    )
  }

  cluster_indices <- reaches_cluster_indices[[k]]

  if (length(cluster_indices) == 0L) {
    plot.new()

    title(
      main = paste0(
        "REACHES Cluster ",
        k,
        " (no grid cells)"
      )
    )

    return(
      invisible(
        NULL
      )
    )
  }

  fda::fbplot(
    t(
      eval_R[
        cluster_indices,
        ,
        drop = FALSE
      ]
    ),
    x = time_grid,
    xlim = range(
      time_grid
    ),
    ylim = ylim_R,
    axes = FALSE,
    main = paste0(
      "REACHES Cluster ",
      k
    ),
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

  invisible(
    NULL
  )
}

save_reaches_clusters_separately <- function(
    output_dir,
    width = 6,
    height = 4,
    res = 300) {

  dir.create(
    output_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  for (
    k in seq_len(
      number_of_clusters
    )
  ) {

    filename <- file.path(
      output_dir,
      paste0(
        "Figure6_cluster_",
        k,
        ".jpg"
      )
    )

    jpeg(
      filename = filename,
      width = width,
      height = height,
      res = res,
      units = "in"
    )

    old_parameters <- par(
      no.readonly = TRUE
    )

    tryCatch(
      {
        par(
          mfrow = c(
            1,
            1
          ),
          mar = c(
            3.5,
            4,
            2,
            1
          ),
          mgp = c(
            2,
            0.7,
            0
          ),
          tcl = -0.2
        )

        draw_reaches_cluster(
          k
        )
      },
      finally = {
        par(
          old_parameters
        )

        dev.off()
      }
    )
  }

  invisible(
    output_dir
  )
}

# ------------------------------------------------------------
# 9. Create the cluster map
# ------------------------------------------------------------

get_grid_size <- function(
    data,
    coordinate) {

  coordinate_values <- sort(
    unique(
      data[[coordinate]]
    )
  )

  coordinate_differences <- diff(
    coordinate_values
  )

  coordinate_differences <- coordinate_differences[
    coordinate_differences > 0
  ]

  if (length(coordinate_differences) == 0L) {
    stop(
      "Unable to determine the grid spacing for ",
      coordinate,
      "."
    )
  }

  min(
    coordinate_differences
  )
}

df_map_reach <- data.frame(
  long = reaches_grid$long,
  lat = reaches_grid$lati,
  cluster = factor(
    cl_R_geo,
    levels = seq_len(
      number_of_clusters
    )
  )
)

grid_width <- get_grid_size(
  df_map_reach,
  "long"
)

grid_height <- get_grid_size(
  df_map_reach,
  "lat"
)

p_reach_map <- ggplot(
  df_map_reach,
  aes(
    x = long,
    y = lat,
    fill = cluster
  )
) +
  geom_tile(
    width = grid_width * 0.98,
    height = grid_height * 0.98
  ) +
  borders(
    "world",
    xlim = c(
      76,
      132
    ),
    ylim = c(
      18,
      52
    ),
    fill = NA,
    colour = "grey40"
  ) +
  coord_map(
    xlim = c(
      98,
      130.5
    ),
    ylim = c(
      18,
      42.5
    )
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Cluster"
  ) +
  theme_minimal()

# ------------------------------------------------------------
# 10. Save Figure 6 panels and assignments
# ------------------------------------------------------------

figure6_output_dir <- here::here(
  "Output",
  "Figure6"
)

dir.create(
  figure6_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

ggsave(
  filename = file.path(
    figure6_output_dir,
    "Figure6_cluster_map.jpg"
  ),
  plot = p_reach_map,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)

save_reaches_clusters_separately(
  output_dir = figure6_output_dir
)

cluster_assignment_output <- reaches_grid %>%
  dplyr::select(
    any_of(
      "location_id"
    ),
    long,
    lati
  ) %>%
  mutate(
    original_mclust_component =
      as.integer(
        cl_R
      ),
    figure6_cluster =
      as.integer(
        cl_R_geo
      )
  )

readr::write_csv(
  cluster_assignment_output,
  file.path(
    figure6_output_dir,
    "Figure6_cluster_assignments.csv"
  )
)

bic_matrix <- as.matrix(
  mc_R_bic$BIC
)

bic_output <- as.data.frame(
  as.table(
    bic_matrix
  ),
  stringsAsFactors = FALSE
)

names(
  bic_output
) <- c(
  "number_of_clusters",
  "model_name",
  "BIC"
)

bic_output <- bic_output %>%
  mutate(
    number_of_clusters =
      as.integer(
        as.character(
          number_of_clusters
        )
      ),
    BIC = as.numeric(
      BIC
    )
  ) %>%
  filter(
    is.finite(
      BIC
    )
  ) %>%
  arrange(
    desc(
      BIC
    )
  )

readr::write_csv(
  bic_output,
  file.path(
    figure6_output_dir,
    "Figure6_mclust_bic.csv"
  )
)

cluster_selection_output <- data.frame(
  cluster_selection_mode =
    cluster_selection_mode,
  bic_selected_number_of_clusters =
    bic_selected_number_of_clusters,
  figure_number_of_clusters =
    number_of_clusters,
  selected_model_name =
    mc_R$modelName,
  relabeling_rule =
    "ascending overall fitted REACHES temperature index"
)

readr::write_csv(
  cluster_selection_output,
  file.path(
    figure6_output_dir,
    "Figure6_cluster_selection.csv"
  )
)

message(
  "All Figure 6 panels, BIC values, and cluster assignments ",
  "were saved to: ",
  figure6_output_dir
)
