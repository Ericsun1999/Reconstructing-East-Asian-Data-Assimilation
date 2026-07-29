here::i_am("Code/Supplementary/FigureS5.R")

# ============================================================
# Supplementary Figure S5
#
# Functional boxplots within the REACHES-defined clusters used
# in Figure 6.
#
# Rows:
#   Figure 6 clusters (currently four)
#
# Columns:
#   1. Kriged REACHES index trajectories
#   2. Cellwise-centered LME ensemble mean
#   3. Cellwise-centered assimilated posterior mean
#
# IMPORTANT:
#   Figure S5 does not repeat FPCA or Mclust. It reuses the exact
#   location set and cluster labels saved by Figure6.R.
#
# Required generated inputs:
#   Output/Figure6/Figure6_cluster_assignments.csv
#   Output/Intermediate/REACHES/
#     reaches_kriging_lme_grid_mean.csv
#   Output/Intermediate/LME/
#     lme_ensemble_mean_1368_1911.csv
#   Output/Intermediate/Assimilation/
#     assimilated_posterior_lme_grid_mean.csv
#
# Equivalent precomputed inputs may be placed under:
#   Data/REACHES/precomputed/
#   Data/LME data/precomputed/
#   Data/Valid/
#
# The posterior-grid file is not produced by the three-city
# Figure9d.R script. A separate all-location assimilation step is
# required before Figure S5 can be regenerated.
#
# Outputs:
#   Output/Supplementary/FigureS5.png
#   Output/Supplementary/FigureS5.pdf
#   Output/Intermediate/FigureS5/
#     FigureS5_cluster_assignments.csv
#     FigureS5_diagnostics.csv
# ============================================================

library(dplyr)
library(readr)
library(fda)

# ------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------

analysis_years <- 1368:1911
input_mode <- "auto"

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

# Figure 6 currently selects four clusters. The plotting code is
# dynamic, but this expected value catches stale Figure 6 output.
expected_number_of_clusters <- 4L

# Keep FALSE when Figure6_cluster_assignments.csv already stores
# the final labels displayed in Figure 6.
#
# Set TRUE only if Figure 6 displayed clusters 1 and 4 in swapped
# order without applying the same swap to its assignment CSV.
swap_cluster_1_and_4 <- FALSE

generated_files <- c(
  clusters = here::here(
    "Output",
    "Figure6",
    "Figure6_cluster_assignments.csv"
  ),
  reaches = here::here(
    "Output",
    "Intermediate",
    "REACHES",
    "reaches_kriging_lme_grid_mean.csv"
  ),
  lme = here::here(
    "Output",
    "Intermediate",
    "LME",
    "lme_ensemble_mean_1368_1911.csv"
  ),
  posterior = here::here(
    "Output",
    "Intermediate",
    "Assimilation",
    "assimilated_posterior_lme_grid_mean.csv"
  )
)

precomputed_files <- c(
  clusters = here::here(
    "Data",
    "REACHES",
    "precomputed",
    "Figure6_cluster_assignments.csv"
  ),
  reaches = here::here(
    "Data",
    "REACHES",
    "precomputed",
    "reaches_kriging_lme_grid_mean.csv"
  ),
  lme = here::here(
    "Data",
    "LME data",
    "precomputed",
    "lme_ensemble_mean_1368_1911.csv"
  ),
  posterior = here::here(
    "Data",
    "Valid",
    "assimilated_posterior_lme_grid_mean.csv"
  )
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

cluster_output_file <- file.path(
  intermediate_dir,
  "FigureS5_cluster_assignments.csv"
)

diagnostic_file <- file.path(
  intermediate_dir,
  "FigureS5_diagnostics.csv"
)

# ------------------------------------------------------------
# 2. Select one complete input set
# ------------------------------------------------------------

select_input_set <- function(
    input_mode,
    generated_files,
    precomputed_files) {

  generated_complete <- all(
    file.exists(
      generated_files
    )
  )

  precomputed_complete <- all(
    file.exists(
      precomputed_files
    )
  )

  if (input_mode == "generated") {
    if (!generated_complete) {
      missing_files <- generated_files[
        !file.exists(
          generated_files
        )
      ]

      stop(
        "The generated Figure S5 input set is incomplete. ",
        "Missing:\n  ",
        paste(
          missing_files,
          collapse = "\n  "
        ),
        "\nThe all-location posterior grid must be generated ",
        "before Figure S5 can run."
      )
    }

    return(
      generated_files
    )
  }

  if (input_mode == "precomputed") {
    if (!precomputed_complete) {
      missing_files <- precomputed_files[
        !file.exists(
          precomputed_files
        )
      ]

      stop(
        "The precomputed Figure S5 input set is incomplete. ",
        "Missing:\n  ",
        paste(
          missing_files,
          collapse = "\n  "
        )
      )
    }

    return(
      precomputed_files
    )
  }

  if (generated_complete) {
    return(
      generated_files
    )
  }

  if (precomputed_complete) {
    return(
      precomputed_files
    )
  }

  generated_missing <- generated_files[
    !file.exists(
      generated_files
    )
  ]

  precomputed_missing <- precomputed_files[
    !file.exists(
      precomputed_files
    )
  ]

  stop(
    "Neither a complete generated nor a complete precomputed ",
    "Figure S5 input set was found.\n",
    "Missing generated files:\n  ",
    paste(
      generated_missing,
      collapse = "\n  "
    ),
    "\nMissing precomputed files:\n  ",
    paste(
      precomputed_missing,
      collapse = "\n  "
    ),
    "\nIn particular, Figure9d.R currently produces only ",
    "three-city posterior series and cannot create the spatial ",
    "posterior input required by Figure S5."
  )
}

input_files <- select_input_set(
  input_mode = input_mode,
  generated_files = generated_files,
  precomputed_files = precomputed_files
)

message(
  "Figure S5 inputs selected from the ",
  if (
    identical(
      unname(
        input_files
      ),
      unname(
        generated_files
      )
    )
  ) {
    "generated"
  } else {
    "precomputed"
  },
  " input set."
)

# ------------------------------------------------------------
# 3. Grid-data helpers
# ------------------------------------------------------------

coordinate_key <- function(
    long,
    lat) {

  paste(
    sprintf(
      "%.8f",
      as.numeric(
        long
      )
    ),
    sprintf(
      "%.8f",
      as.numeric(
        lat
      )
    ),
    sep = "_"
  )
}

standardize_grid_data <- function(
    input_file,
    object_name) {

  grid_data <- readr::read_csv(
    input_file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  if (
    !"lat" %in%
      names(
        grid_data
      ) &&
      "lati" %in%
        names(
          grid_data
        )
  ) {
    grid_data <- grid_data %>%
      rename(
        lat = lati
      )
  }

  required_coordinate_columns <- c(
    "long",
    "lat"
  )

  missing_coordinate_columns <- setdiff(
    required_coordinate_columns,
    names(
      grid_data
    )
  )

  if (length(
    missing_coordinate_columns
  ) > 0L) {
    stop(
      object_name,
      " is missing coordinate columns: ",
      paste(
        missing_coordinate_columns,
        collapse = ", "
      )
    )
  }

  year_columns <- grep(
    "^[Xx]?[0-9]{4}$",
    names(
      grid_data
    ),
    value = TRUE
  )

  if (length(
    year_columns
  ) == 0L) {
    stop(
      object_name,
      " contains no annual columns."
    )
  }

  year_values <- as.integer(
    sub(
      "^[Xx]",
      "",
      year_columns
    )
  )

  names(
    grid_data
  )[
    match(
      year_columns,
      names(
        grid_data
      )
    )
  ] <- paste0(
    "x",
    year_values
  )

  grid_data <- grid_data %>%
    mutate(
      long = as.numeric(
        long
      ),
      lat = as.numeric(
        lat
      ),
      location_id = if (
        "location_id" %in%
          names(
            grid_data
          )
      ) {
        as.integer(
          location_id
        )
      } else {
        NA_integer_
      },
      cell_id = coordinate_key(
        long,
        lat
      )
    )

  if (
    any(!is.finite(
      grid_data$long
    )) ||
      any(!is.finite(
        grid_data$lat
      )) ||
      anyDuplicated(
        grid_data$cell_id
      )
  ) {
    stop(
      object_name,
      " contains invalid or duplicated coordinates."
    )
  }

  grid_data
}

extract_years <- function(
    grid_data) {

  year_columns <- grep(
    "^x[0-9]{4}$",
    names(
      grid_data
    ),
    value = TRUE
  )

  as.integer(
    sub(
      "^x",
      "",
      year_columns
    )
  )
}

center_rows <- function(
    matrix_data) {

  row_means <- rowMeans(
    matrix_data
  )

  sweep(
    matrix_data,
    1,
    row_means,
    FUN = "-"
  )
}

# ------------------------------------------------------------
# 4. Read Figure 6 assignments and the three spatial products
# ------------------------------------------------------------

cluster_assignments <- readr::read_csv(
  unname(
    input_files[["clusters"]]
  ),
  show_col_types = FALSE
)

if (
  !"lat" %in%
    names(
      cluster_assignments
    ) &&
    "lati" %in%
      names(
        cluster_assignments
      )
) {
  cluster_assignments <- cluster_assignments %>%
    rename(
      lat = lati
    )
}

required_cluster_columns <- c(
  "long",
  "lat",
  "figure6_cluster"
)

missing_cluster_columns <- setdiff(
  required_cluster_columns,
  names(
    cluster_assignments
  )
)

if (length(
  missing_cluster_columns
) > 0L) {
  stop(
    "Figure6_cluster_assignments.csv is missing: ",
    paste(
      missing_cluster_columns,
      collapse = ", "
    )
  )
}

cluster_assignments <- cluster_assignments %>%
  transmute(
    location_id = if (
      "location_id" %in%
        names(
          cluster_assignments
        )
    ) {
      as.integer(
        location_id
      )
    } else {
      NA_integer_
    },
    long = as.numeric(
      long
    ),
    lat = as.numeric(
      lat
    ),
    original_mclust_component = if (
      "original_mclust_component" %in%
        names(
          cluster_assignments
        )
    ) {
      as.integer(
        original_mclust_component
      )
    } else {
      NA_integer_
    },
    cluster = as.integer(
      figure6_cluster
    ),
    cell_id = coordinate_key(
      long,
      lat
    )
  )

if (swap_cluster_1_and_4) {
  cluster_assignments <- cluster_assignments %>%
    mutate(
      cluster = dplyr::recode(
        cluster,
        `1` = 4L,
        `4` = 1L,
        .default = cluster
      )
    )
}

if (
  anyNA(
    cluster_assignments[
      ,
      c(
        "long",
        "lat",
        "cluster",
        "cell_id"
      )
    ]
  ) ||
    anyDuplicated(
      cluster_assignments$cell_id
    )
) {
  stop(
    "The Figure 6 cluster assignments contain invalid or ",
    "duplicated locations."
  )
}

cluster_levels <- sort(
  unique(
    cluster_assignments$cluster
  )
)

number_of_clusters <- length(
  cluster_levels
)

if (!identical(
  cluster_levels,
  seq_len(
    number_of_clusters
  )
)) {
  stop(
    "Figure 6 cluster labels must be consecutive integers ",
    "starting at 1."
  )
}

if (
  number_of_clusters !=
    expected_number_of_clusters
) {
  warning(
    "Expected ",
    expected_number_of_clusters,
    " Figure 6 clusters, but found ",
    number_of_clusters,
    ". The layout will use the available cluster count."
  )
}

reaches_grid <- standardize_grid_data(
  unname(
    input_files[["reaches"]]
  ),
  "REACHES-at-LME-grid field"
)

lme_grid <- standardize_grid_data(
  unname(
    input_files[["lme"]]
  ),
  "LME ensemble-mean field"
)

posterior_grid <- standardize_grid_data(
  unname(
    input_files[["posterior"]]
  ),
  "Assimilated posterior field"
)

# ------------------------------------------------------------
# 5. Align all products to the exact Figure 6 locations
# ------------------------------------------------------------

align_grid_to_clusters <- function(
    grid_data,
    cluster_assignments,
    object_name) {

  # Prefer location_id when both objects contain complete,
  # unique identifiers. Otherwise align by coordinates.
  use_location_id <-
    all(
      !is.na(
        cluster_assignments$location_id
      )
    ) &&
    all(
      !is.na(
        grid_data$location_id
      )
    ) &&
    !anyDuplicated(
      cluster_assignments$location_id
    ) &&
    !anyDuplicated(
      grid_data$location_id
    )

  if (use_location_id) {
    match_rows <- match(
      cluster_assignments$location_id,
      grid_data$location_id
    )
  } else {
    match_rows <- match(
      cluster_assignments$cell_id,
      grid_data$cell_id
    )
  }

  if (anyNA(
    match_rows
  )) {
    stop(
      sum(
        is.na(
          match_rows
        )
      ),
      " Figure 6 locations are absent from ",
      object_name,
      "."
    )
  }

  aligned <- grid_data[
    match_rows,
    ,
    drop = FALSE
  ]

  if (!identical(
    aligned$cell_id,
    cluster_assignments$cell_id
  )) {
    stop(
      object_name,
      " could not be aligned exactly to the Figure 6 ",
      "coordinates."
    )
  }

  aligned
}

reaches_aligned <- align_grid_to_clusters(
  reaches_grid,
  cluster_assignments,
  "the REACHES field"
)

lme_aligned <- align_grid_to_clusters(
  lme_grid,
  cluster_assignments,
  "the LME field"
)

posterior_aligned <- align_grid_to_clusters(
  posterior_grid,
  cluster_assignments,
  "the posterior field"
)

# ------------------------------------------------------------
# 6. Determine the common REACHES event-year grid
# ------------------------------------------------------------

event_years <- sort(
  Reduce(
    intersect,
    list(
      extract_years(
        reaches_aligned
      ),
      extract_years(
        lme_aligned
      ),
      extract_years(
        posterior_aligned
      ),
      analysis_years
    )
  )
)

if (length(
  event_years
) < 2L) {
  stop(
    "Fewer than two common event years were found across the ",
    "three spatial products."
  )
}

year_columns <- paste0(
  "x",
  event_years
)

Y_reaches <- as.matrix(
  reaches_aligned[
    ,
    year_columns,
    drop = FALSE
  ]
)

Y_lme <- as.matrix(
  lme_aligned[
    ,
    year_columns,
    drop = FALSE
  ]
)

Y_posterior <- as.matrix(
  posterior_aligned[
    ,
    year_columns,
    drop = FALSE
  ]
)

storage.mode(
  Y_reaches
) <- "double"

storage.mode(
  Y_lme
) <- "double"

storage.mode(
  Y_posterior
) <- "double"

if (
  any(!is.finite(
    Y_reaches
  )) ||
    any(!is.finite(
      Y_lme
    )) ||
    any(!is.finite(
      Y_posterior
    ))
) {
  stop(
    "At least one aligned trajectory contains missing or ",
    "non-finite values over the common event years."
  )
}

Y_lme_centered <- center_rows(
  Y_lme
)

Y_posterior_centered <- center_rows(
  Y_posterior
)

# ------------------------------------------------------------
# 7. Smooth the three sets of trajectories
# ------------------------------------------------------------

make_evaluated_functions <- function(
    trajectory_matrix,
    years,
    number_of_basis_functions) {

  basis <- fda::create.bspline.basis(
    rangeval = range(
      years
    ),
    nbasis = number_of_basis_functions
  )

  functional_data <- fda::Data2fd(
    argvals = years,
    y = t(
      trajectory_matrix
    ),
    basisobj = basis
  )

  t(
    fda::eval.fd(
      years,
      functional_data
    )
  )
}

reaches_evaluated <- make_evaluated_functions(
  trajectory_matrix = Y_reaches,
  years = event_years,
  number_of_basis_functions = 20
)

lme_evaluated <- make_evaluated_functions(
  trajectory_matrix = Y_lme_centered,
  years = event_years,
  number_of_basis_functions = 30
)

posterior_evaluated <- make_evaluated_functions(
  trajectory_matrix = Y_posterior_centered,
  years = event_years,
  number_of_basis_functions = 30
)

# ------------------------------------------------------------
# 8. Draw Figure S5
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
    title = "REACHES kriged\nindex"
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
    title = "LME ensemble mean\n(cellwise centered)"
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
    title = "Assimilated posterior\n(cellwise centered)"
  )
)

draw_functional_boxplot <- function(
    evaluated_matrix,
    selected_rows,
    y_limits,
    y_ticks,
    show_x_labels) {

  if (length(
    selected_rows
  ) < 2L) {
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
      number_of_clusters,
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

  for (
    cluster_id in cluster_levels
  ) {

    selected_rows <- which(
      cluster_assignments$cluster ==
        cluster_id
    )

    for (
      column_index in seq_along(
        settings_order
      )
    ) {

      setting <- column_settings[
        [
          settings_order[
            column_index
          ]
        ]
      ]

      draw_functional_boxplot(
        evaluated_matrix =
          setting$matrix,
        selected_rows = selected_rows,
        y_limits = setting$limits,
        y_ticks = setting$ticks,
        show_x_labels =
          cluster_id ==
          max(
            cluster_levels
          )
      )

      if (
        cluster_id ==
          min(
            cluster_levels
          )
      ) {
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

figure_height <- 2.45 *
  number_of_clusters +
  0.5

png(
  filename = figure_png,
  width = 10,
  height = figure_height,
  units = "in",
  res = 300
)

draw_complete_figure()

dev.off()

pdf(
  file = figure_pdf,
  width = 10,
  height = figure_height,
  onefile = TRUE
)

draw_complete_figure()

dev.off()

# ------------------------------------------------------------
# 9. Save assignments and diagnostics
# ------------------------------------------------------------

cluster_output <- cluster_assignments %>%
  mutate(
    number_of_event_years =
      length(
        event_years
      )
  )

readr::write_csv(
  cluster_output,
  cluster_output_file
)

cluster_counts <- as.numeric(
  table(
    factor(
      cluster_assignments$cluster,
      levels = cluster_levels
    )
  )
)

diagnostics <- data.frame(
  quantity = c(
    "number_of_event_years",
    "number_of_figure6_locations",
    "number_of_clusters",
    paste0(
      "cluster_",
      cluster_levels,
      "_locations"
    ),
    "swap_cluster_1_and_4"
  ),
  value = c(
    as.character(
      length(
        event_years
      )
    ),
    as.character(
      nrow(
        cluster_assignments
      )
    ),
    as.character(
      number_of_clusters
    ),
    as.character(
      cluster_counts
    ),
    as.character(
      swap_cluster_1_and_4
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
  "Saved Figure S5 cluster assignments: ",
  cluster_output_file
)
