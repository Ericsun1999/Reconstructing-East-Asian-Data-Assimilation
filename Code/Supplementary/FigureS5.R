here::i_am("Code/Supplementary/FigureS5.R")

# ============================================================
# Supplementary Figure S5
#
# Functional boxplots within the five REACHES-defined clusters.
#
# Rows:
#   REACHES clusters 1--5
#
# Columns:
#   1. Kriged REACHES anomalies
#   2. Cellwise-centered LME ensemble mean
#   3. Cellwise-centered assimilated posterior mean
#
# Clustering is performed ONLY on the kriged REACHES field.
# LME and posterior trajectories are summarized within those
# same REACHES-defined regions and are not reclustered.
# ============================================================

library(here)
library(readxl)
library(readr)
library(dplyr)
library(sf)
library(fda)
library(mclust)

# ------------------------------------------------------------
# 1. Input and output locations
# ------------------------------------------------------------

reaches_grid_candidates <- c(
  here::here("Data", "reaches_kriging_grid53x49_mean.csv"),
  here::here("Data", "tempe_all_v4.csv")
)

posterior_grid_candidates <- c(
  here::here("Data", "assimilated_posterior_grid53x49_mean.csv"),
  here::here("Data", "kalman_mean_v4.csv")
)

lme_mean_candidates <- c(
  here::here("Output", "Intermediate", "FigureS5",
             "lme_ensemble_mean_grid53x49.csv"),
  here::here("Data", "lme_ensemble_mean_grid53x49.csv"),
  here::here("Data", "LME data", "FigureS5",
             "lme_ensemble_mean_grid53x49.csv")
)

lme_member_directories <- c(
  here::here("Data", "LME data", "FigureS5"),
  here::here("Data", "LME data", "Figure6"),
  here::here("Data", "LME data")
)

reaches_observation_file <- here::here(
  "Data",
  "temperature index value.v1.xlsx"
)

output_dir <- here::here(
  "Output",
  "Supplementary"
)

intermediate_dir <- here::here(
  "Output",
  "Intermediate",
  "FigureS5"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  intermediate_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

figure_png <- file.path(
  output_dir,
  "FigureS5.png"
)

figure_pdf <- file.path(
  output_dir,
  "FigureS5.pdf"
)

cluster_file <- file.path(
  intermediate_dir,
  "reaches_cluster_assignments.csv"
)

diagnostic_file <- file.path(
  intermediate_dir,
  "FigureS5_diagnostics.csv"
)

analysis_years <- 1368:1911

# Relabeling used for the geographic cluster order in Figure 6.
cluster_relabel <- c(
  "2" = 1,
  "4" = 2,
  "5" = 3,
  "1" = 4,
  "3" = 5
)

# ------------------------------------------------------------
# 2. General helpers
# ------------------------------------------------------------

first_existing_file <- function(
    candidates,
    description) {

  existing <- candidates[
    file.exists(candidates)
  ]

  if (length(existing) == 0L) {
    stop(
      description,
      " was not found. Searched:\n",
      paste(
        paste0("  - ", candidates),
        collapse = "\n"
      )
    )
  }

  existing[1]
}


standardize_grid_data <- function(
    data,
    object_name) {

  names(data) <- sub(
    "^X(?=[0-9]{4}$)",
    "x",
    names(data),
    perl = TRUE
  )

  if (
    ncol(data) > 0L &&
      names(data)[1] %in% c(
        "X",
        "...1",
        "row.names"
      )
  ) {
    first_column <- suppressWarnings(
      as.integer(data[[1]])
    )

    if (
      length(first_column) == nrow(data) &&
        all(
          first_column ==
            seq_len(nrow(data)),
          na.rm = TRUE
        )
    ) {
      data <- data[
        ,
        -1,
        drop = FALSE
      ]
    }
  }

  longitude_name <- intersect(
    c(
      "long",
      "lon",
      "longitude"
    ),
    names(data)
  )

  latitude_name <- intersect(
    c(
      "lat",
      "lati",
      "latitude"
    ),
    names(data)
  )

  if (
    length(longitude_name) == 0L ||
      length(latitude_name) == 0L
  ) {
    stop(
      object_name,
      " must contain longitude and latitude columns."
    )
  }

  names(data)[
    names(data) ==
      longitude_name[1]
  ] <- "long"

  names(data)[
    names(data) ==
      latitude_name[1]
  ] <- "lat"

  data <- data %>%
    mutate(
      long = as.numeric(long),
      lat = as.numeric(lat)
    ) %>%
    filter(
      is.finite(long),
      is.finite(lat)
    ) %>%
    distinct(
      long,
      lat,
      .keep_all = TRUE
    )

  year_columns <- grep(
    "^x[0-9]{4}$",
    names(data),
    value = TRUE
  )

  if (length(year_columns) == 0L) {
    stop(
      object_name,
      " does not contain columns named xYYYY."
    )
  }

  data
}


extract_year_map <- function(
    data) {

  columns <- grep(
    "^x[0-9]{4}$",
    names(data),
    value = TRUE
  )

  years <- as.integer(
    sub(
      "^x",
      "",
      columns
    )
  )

  data.frame(
    year = years,
    column = columns
  ) %>%
    arrange(year)
}


coordinate_key <- function(
    long,
    lat) {

  paste(
    sprintf(
      "%.6f",
      long
    ),
    sprintf(
      "%.6f",
      lat
    ),
    sep = "_"
  )
}


center_rows <- function(
    matrix_data) {

  sweep(
    matrix_data,
    1,
    rowMeans(
      matrix_data,
      na.rm = TRUE
    ),
    FUN = "-"
  )
}


find_lme_member_file <- function(
    directory,
    member_id) {

  candidates <- c(
    file.path(
      directory,
      paste0(
        "c",
        member_id,
        ".csv.gz"
      )
    ),
    file.path(
      directory,
      paste0(
        "c",
        member_id,
        ".csv"
      )
    )
  )

  existing <- candidates[
    file.exists(candidates)
  ]

  if (length(existing) == 0L) {
    return(NA_character_)
  }

  existing[1]
}


read_lme_ensemble_mean <- function() {

  existing_mean <- lme_mean_candidates[
    file.exists(
      lme_mean_candidates
    )
  ]

  if (length(existing_mean) > 0L) {

    message(
      "Using archived LME ensemble mean: ",
      existing_mean[1]
    )

    return(
      standardize_grid_data(
        read.csv(
          existing_mean[1],
          check.names = FALSE
        ),
        "LME ensemble mean"
      )
    )
  }

  valid_directory <- vapply(
    lme_member_directories,
    function(directory) {

      files <- vapply(
        1:13,
        function(member_id) {
          find_lme_member_file(
            directory,
            member_id
          )
        },
        character(1)
      )

      all(
        !is.na(files)
      )
    },
    logical(1)
  )

  if (!any(valid_directory)) {
    stop(
      paste0(
        "Neither an archived LME ensemble-mean file nor ",
        "c1--c13 files were found.\n",
        "Directories searched for c1--c13:\n",
        paste(
          paste0(
            "  - ",
            lme_member_directories
          ),
          collapse = "\n"
        )
      )
    )
  }

  lme_directory <- lme_member_directories[
    which(valid_directory)[1]
  ]

  message(
    "Constructing the LME ensemble mean from: ",
    lme_directory
  )

  member_data <- lapply(
    1:13,
    function(member_id) {

      input_file <- find_lme_member_file(
        lme_directory,
        member_id
      )

      standardize_grid_data(
        read.csv(
          input_file,
          check.names = FALSE
        ),
        paste0(
          "LME member c",
          member_id
        )
      )
    }
  )

  reference_keys <- coordinate_key(
    member_data[[1]]$long,
    member_data[[1]]$lat
  )

  reference_year_map <- extract_year_map(
    member_data[[1]]
  )

  for (member_id in 2:13) {

    member_keys <- coordinate_key(
      member_data[[member_id]]$long,
      member_data[[member_id]]$lat
    )

    if (!setequal(
      reference_keys,
      member_keys
    )) {
      stop(
        "The spatial grid in LME member c",
        member_id,
        " does not match c1."
      )
    }

    member_data[[member_id]] <- member_data[[member_id]][
      match(
        reference_keys,
        member_keys
      ),
      ,
      drop = FALSE
    ]

    member_years <- extract_year_map(
      member_data[[member_id]]
    )$year

    if (!identical(
      reference_year_map$year,
      member_years
    )) {
      stop(
        "The year columns in LME member c",
        member_id,
        " do not match c1."
      )
    }
  }

  year_columns <- reference_year_map$column

  lme_sum <- matrix(
    0,
    nrow = nrow(
      member_data[[1]]
    ),
    ncol = length(
      year_columns
    )
  )

  for (member_id in 1:13) {

    current_matrix <- data.matrix(
      member_data[[member_id]][
        ,
        year_columns,
        drop = FALSE
      ]
    )

    if (anyNA(current_matrix)) {
      stop(
        "Missing LME values were found in member c",
        member_id,
        "."
      )
    }

    lme_sum <- lme_sum +
      current_matrix
  }

  lme_mean <- data.frame(
    long = member_data[[1]]$long,
    lat = member_data[[1]]$lat,
    lme_sum / 13,
    check.names = FALSE
  )

  names(lme_mean)[
    3:ncol(lme_mean)
  ] <- year_columns

  readr::write_csv(
    lme_mean,
    file.path(
      intermediate_dir,
      "lme_ensemble_mean_grid53x49.csv"
    )
  )

  lme_mean
}

# ------------------------------------------------------------
# 3. Read the three gridded products
# ------------------------------------------------------------

reaches_grid_file <- first_existing_file(
  reaches_grid_candidates,
  "The kriged REACHES grid"
)

posterior_grid_file <- first_existing_file(
  posterior_grid_candidates,
  "The assimilated posterior grid"
)

if (!file.exists(
  reaches_observation_file
)) {
  stop(
    "The original REACHES file was not found: ",
    reaches_observation_file
  )
}

message(
  "Using REACHES grid: ",
  reaches_grid_file
)

message(
  "Using posterior grid: ",
  posterior_grid_file
)

reaches_grid_full <- standardize_grid_data(
  read.csv(
    reaches_grid_file,
    check.names = FALSE
  ),
  "Kriged REACHES grid"
)

posterior_grid_full <- standardize_grid_data(
  read.csv(
    posterior_grid_file,
    check.names = FALSE
  ),
  "Assimilated posterior grid"
)

lme_grid_full <- read_lme_ensemble_mean()

# ------------------------------------------------------------
# 4. Determine REACHES event years shared by all products
# ------------------------------------------------------------

reaches_year_map <- extract_year_map(
  reaches_grid_full
)

posterior_year_map <- extract_year_map(
  posterior_grid_full
)

lme_year_map <- extract_year_map(
  lme_grid_full
)

event_years <- sort(
  Reduce(
    intersect,
    list(
      reaches_year_map$year,
      posterior_year_map$year,
      lme_year_map$year,
      analysis_years
    )
  )
)

if (length(event_years) == 0L) {
  stop(
    "The three gridded products have no common years ",
    "between 1368 and 1911."
  )
}

if (length(event_years) != 524L) {
  warning(
    "The supplementary manuscript describes 524 REACHES ",
    "event years, but ",
    length(event_years),
    " common years were found."
  )
}

year_columns <- paste0(
  "x",
  event_years
)

# ------------------------------------------------------------
# 5. Identify REACHES-supported grid cells
# ------------------------------------------------------------

temperature <- readxl::read_excel(
  reaches_observation_file,
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

names(temperature) <- c(
  "level",
  "year",
  "long",
  "lat"
)

temperature <- temperature %>%
  transmute(
    level = as.numeric(level),
    year = as.integer(year),
    long = as.numeric(long),
    lat = as.numeric(lat)
  ) %>%
  filter(
    year %in% event_years,
    is.finite(level),
    is.finite(long),
    is.finite(lat)
  ) %>%
  distinct(
    year,
    long,
    lat,
    .keep_all = TRUE
  )

observation_sf <- sf::st_as_sf(
  temperature,
  coords = c(
    "long",
    "lat"
  ),
  crs = 4326,
  remove = FALSE
)

reaches_grid_sf <- reaches_grid_full %>%
  mutate(
    grid_row = row_number()
  ) %>%
  sf::st_as_sf(
    coords = c(
      "long",
      "lat"
    ),
    crs = 4326,
    remove = FALSE
  )

observation_projected <- sf::st_transform(
  observation_sf,
  3857
)

grid_projected <- sf::st_transform(
  reaches_grid_sf,
  3857
)

nearest_grid_rows <- sf::st_nearest_feature(
  observation_projected,
  grid_projected
)

supported_grid_rows <- sort(
  unique(
    nearest_grid_rows
  )
)

reaches_supported <- reaches_grid_full[
  supported_grid_rows,
  ,
  drop = FALSE
]

reaches_supported$cell_id <- coordinate_key(
  reaches_supported$long,
  reaches_supported$lat
)

# ------------------------------------------------------------
# 6. Align posterior and LME to the supported REACHES cells
# ------------------------------------------------------------

posterior_grid_full$cell_id <- coordinate_key(
  posterior_grid_full$long,
  posterior_grid_full$lat
)

lme_grid_full$cell_id <- coordinate_key(
  lme_grid_full$long,
  lme_grid_full$lat
)

posterior_match <- match(
  reaches_supported$cell_id,
  posterior_grid_full$cell_id
)

lme_match <- match(
  reaches_supported$cell_id,
  lme_grid_full$cell_id
)

if (anyNA(posterior_match)) {
  stop(
    sum(
      is.na(
        posterior_match
      )
    ),
    " REACHES-supported cells are missing from the posterior grid."
  )
}

if (anyNA(lme_match)) {
  stop(
    sum(
      is.na(
        lme_match
      )
    ),
    " REACHES-supported cells are missing from the LME grid."
  )
}

posterior_supported <- posterior_grid_full[
  posterior_match,
  ,
  drop = FALSE
]

lme_supported <- lme_grid_full[
  lme_match,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    reaches_supported$cell_id,
    coordinate_key(
      posterior_supported$long,
      posterior_supported$lat
    )
  ),
  identical(
    reaches_supported$cell_id,
    coordinate_key(
      lme_supported$long,
      lme_supported$lat
    )
  )
)

# ------------------------------------------------------------
# 7. Construct the trajectory matrices
# ------------------------------------------------------------

Y_reaches <- data.matrix(
  reaches_supported[
    ,
    year_columns,
    drop = FALSE
  ]
)

Y_lme <- data.matrix(
  lme_supported[
    ,
    year_columns,
    drop = FALSE
  ]
)

Y_posterior <- data.matrix(
  posterior_supported[
    ,
    year_columns,
    drop = FALSE
  ]
)

if (
  anyNA(Y_reaches) ||
    anyNA(Y_lme) ||
    anyNA(Y_posterior)
) {
  stop(
    "At least one gridded trajectory contains missing values ",
    "over the common event years."
  )
}

# REACHES is already on a zero-mean anomaly-index scale.
# LME and posterior are centered separately within each cell.
Y_lme_centered <- center_rows(
  Y_lme
)

Y_posterior_centered <- center_rows(
  Y_posterior
)

# ------------------------------------------------------------
# 8. Define the five clusters from REACHES only
# ------------------------------------------------------------

reaches_basis <- fda::create.bspline.basis(
  rangeval = range(
    event_years
  ),
  nbasis = 20
)

reaches_fd <- fda::Data2fd(
  argvals = event_years,
  y = t(
    Y_reaches
  ),
  basisobj = reaches_basis
)

reaches_fpca <- fda::pca.fd(
  reaches_fd,
  nharm = 5
)

reaches_mclust <- mclust::Mclust(
  reaches_fpca$scores,
  G = 1:8,
  verbose = FALSE
)

if (
  length(reaches_mclust$G) != 1L ||
    reaches_mclust$G != 5L
) {
  stop(
    "BIC selected ",
    paste(
      reaches_mclust$G,
      collapse = ", "
    ),
    " REACHES clusters rather than five."
  )
}

original_cluster <- reaches_mclust$classification

reaches_cluster <- unname(
  cluster_relabel[
    as.character(
      original_cluster
    )
  ]
)

if (anyNA(reaches_cluster)) {
  stop(
    "The saved geographic relabeling does not cover all ",
    "Mclust cluster labels."
  )
}

cluster_assignments <- data.frame(
  cell_id = reaches_supported$cell_id,
  long = reaches_supported$long,
  lat = reaches_supported$lat,
  original_cluster = original_cluster,
  cluster = reaches_cluster
) %>%
  arrange(
    cluster,
    lat,
    long
  )

readr::write_csv(
  cluster_assignments,
  cluster_file
)

# ------------------------------------------------------------
# 9. Smooth trajectories for the functional boxplots
# ------------------------------------------------------------

lme_basis <- fda::create.bspline.basis(
  rangeval = range(
    event_years
  ),
  nbasis = 30
)

posterior_basis <- fda::create.bspline.basis(
  rangeval = range(
    event_years
  ),
  nbasis = 30
)

lme_fd <- fda::Data2fd(
  argvals = event_years,
  y = t(
    Y_lme_centered
  ),
  basisobj = lme_basis
)

posterior_fd <- fda::Data2fd(
  argvals = event_years,
  y = t(
    Y_posterior_centered
  ),
  basisobj = posterior_basis
)

reaches_evaluated <- t(
  fda::eval.fd(
    event_years,
    reaches_fd
  )
)

lme_evaluated <- t(
  fda::eval.fd(
    event_years,
    lme_fd
  )
)

posterior_evaluated <- t(
  fda::eval.fd(
    event_years,
    posterior_fd
  )
)

# ------------------------------------------------------------
# 10. Draw Figure S5
# ------------------------------------------------------------

year_ticks <- seq(
  1400,
  1900,
  by = 100
)

column_settings <- list(
  reaches = list(
    matrix = reaches_evaluated,
    limits = c(
      -1.7,
      0.3
    ),
    ticks = c(
      -1.5,
      -1,
      -0.5,
      0
    ),
    title = "REACHES kriged\nanomaly"
  ),
  lme = list(
    matrix = lme_evaluated,
    limits = c(
      -0.5,
      0.5
    ),
    ticks = c(
      -0.5,
      0,
      0.5
    ),
    title = "LME ensemble mean\n(centered)"
  ),
  posterior = list(
    matrix = posterior_evaluated,
    limits = c(
      -1,
      1
    ),
    ticks = c(
      -1,
      0,
      1
    ),
    title = "Assimilated posterior\n(centered)"
  )
)


draw_functional_boxplot <- function(
    evaluated_matrix,
    selected_rows,
    y_limits,
    y_ticks,
    show_x_labels = TRUE) {

  if (length(selected_rows) < 2L) {
    stop(
      "A functional boxplot requires at least two trajectories."
    )
  }

  fda::fbplot(
    t(
      evaluated_matrix[
        selected_rows,
        ,
        drop = FALSE
      ]
    ),
    x = event_years,
    xlim = range(
      event_years
    ),
    ylim = y_limits,
    axes = FALSE,
    main = "",
    xlab = "",
    ylab = ""
  )

  axis(
    side = 1,
    at = year_ticks,
    labels = if (
      show_x_labels
    ) {
      year_ticks
    } else {
      FALSE
    },
    cex.axis = 0.72
  )

  axis(
    side = 2,
    at = y_ticks,
    labels = y_ticks,
    las = 1,
    cex.axis = 0.72
  )

  box()
}


draw_complete_figure <- function() {

  old_parameters <- par(
    no.readonly = TRUE
  )

  on.exit(
    par(
      old_parameters
    ),
    add = TRUE
  )

  par(
    mfrow = c(
      5,
      3
    ),
    mar = c(
      2.2,
      2.8,
      1.2,
      0.6
    ),
    oma = c(
      0.2,
      3.4,
      2.7,
      0.2
    ),
    mgp = c(
      1.4,
      0.42,
      0
    ),
    tcl = -0.2
  )

  settings_order <- c(
    "reaches",
    "lme",
    "posterior"
  )

  for (cluster_id in 1:5) {

    selected_rows <- which(
      reaches_cluster ==
        cluster_id
    )

    for (column_index in seq_along(
      settings_order
    )) {

      setting <- column_settings[
        [
          settings_order[
            column_index
          ]
        ]
      ]

      draw_functional_boxplot(
        evaluated_matrix = setting$matrix,
        selected_rows = selected_rows,
        y_limits = setting$limits,
        y_ticks = setting$ticks,
        show_x_labels = TRUE
      )

      if (cluster_id == 1L) {
        mtext(
          setting$title,
          side = 3,
          line = 1.25,
          cex = 0.9,
          font = 2
        )
      }

      if (column_index == 1L) {
        mtext(
          paste(
            "Cluster",
            cluster_id
          ),
          side = 2,
          line = 3.9,
          cex = 0.9,
          font = 2
        )
      }
    }
  }
}


png(
  filename = figure_png,
  width = 10,
  height = 12.5,
  units = "in",
  res = 300
)

draw_complete_figure()

dev.off()


pdf(
  file = figure_pdf,
  width = 10,
  height = 12.5,
  onefile = TRUE
)

draw_complete_figure()

dev.off()

# ------------------------------------------------------------
# 11. Save diagnostics
# ------------------------------------------------------------

diagnostics <- data.frame(
  quantity = c(
    "number_of_event_years",
    "number_of_reaches_supported_cells",
    "selected_number_of_clusters",
    "mclust_model",
    paste0(
      "cluster_",
      1:5,
      "_cells"
    )
  ),
  value = c(
    length(
      event_years
    ),
    nrow(
      reaches_supported
    ),
    reaches_mclust$G,
    reaches_mclust$modelName,
    as.numeric(
      table(
        factor(
          reaches_cluster,
          levels = 1:5
        )
      )
    )
  )
)

readr::write_csv(
  diagnostics,
  diagnostic_file
)

message(
  "Saved Figure S5 PNG: ",
  figure_png
)

message(
  "Saved Figure S5 PDF: ",
  figure_pdf
)

message(
  "Saved cluster assignments: ",
  cluster_file
)
