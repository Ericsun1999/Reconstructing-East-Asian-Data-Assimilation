here::i_am("Code/Figure9abc.R")

# ============================================================
# Generate the prior-parameter plots for Beijing, Shanghai, and
# Hong Kong.
#
# Supported input modes:
#   1. "generated":
#      Output/Intermediate/Prior/
#        mtB.csv, muB.csv, rtB.csv
#        mtS.csv, muS.csv, rtS.csv
#        mtH.csv, muH.csv, rtH.csv
#
#   2. "precomputed":
#      Data/par/
#        mtB.csv, muB.csv, rtB.csv
#        mtS.csv, muS.csv, rtS.csv
#        mtH.csv, muH.csv, rtH.csv
#
#   3. "auto" (default):
#      Use the complete generated set when available; otherwise
#      use the complete precomputed set.
#
# Outputs:
#   Beijing:
#     Output/Figure9/Figure9a.png
#     Output/Figure9/Figure9b.png
#     Output/Figure9/Figure9c.png
#
#   Shanghai:
#     Output/Supplementary/FigureS4a.png
#     Output/Supplementary/FigureS4b.png
#     Output/Supplementary/FigureS4c.png
#
#   Hong Kong:
#     Output/Supplementary/FigureS5a.png
#     Output/Supplementary/FigureS5b.png
#     Output/Supplementary/FigureS5c.png
#
#   Numerical plotting data:
#     Output/Figure9/Figure9abc_parameter_plot_data.csv
# ============================================================

library(dplyr)
library(readr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Input and output configuration
# ------------------------------------------------------------

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

generated_parameter_dir <- here::here(
  "Output",
  "Intermediate",
  "Prior"
)

precomputed_parameter_dir <- here::here(
  "Data",
  "par"
)

required_parameter_files <- c(
  "mtB.csv",
  "muB.csv",
  "rtB.csv",
  "mtS.csv",
  "muS.csv",
  "rtS.csv",
  "mtH.csv",
  "muH.csv",
  "rtH.csv"
)

is_complete_parameter_directory <- function(
    directory,
    required_files) {

  all(
    file.exists(
      file.path(
        directory,
        required_files
      )
    )
  )
}

select_parameter_directory <- function(
    input_mode,
    generated_directory,
    precomputed_directory,
    required_files) {

  generated_complete <-
    is_complete_parameter_directory(
      generated_directory,
      required_files
    )

  precomputed_complete <-
    is_complete_parameter_directory(
      precomputed_directory,
      required_files
    )

  if (input_mode == "generated") {
    if (!generated_complete) {
      missing_files <- file.path(
        generated_directory,
        required_files
      )

      missing_files <- missing_files[
        !file.exists(
          missing_files
        )
      ]

      stop(
        "The generated parameter set is incomplete. Missing:\n  ",
        paste(
          missing_files,
          collapse = "\n  "
        )
      )
    }

    return(
      generated_directory
    )
  }

  if (input_mode == "precomputed") {
    if (!precomputed_complete) {
      missing_files <- file.path(
        precomputed_directory,
        required_files
      )

      missing_files <- missing_files[
        !file.exists(
          missing_files
        )
      ]

      stop(
        "The precomputed parameter set is incomplete. Missing:\n  ",
        paste(
          missing_files,
          collapse = "\n  "
        )
      )
    }

    return(
      precomputed_directory
    )
  }

  if (generated_complete) {
    return(
      generated_directory
    )
  }

  if (precomputed_complete) {
    return(
      precomputed_directory
    )
  }

  stop(
    "Neither a complete generated nor a complete precomputed ",
    "parameter set was found."
  )
}

parameter_dir <- select_parameter_directory(
  input_mode = input_mode,
  generated_directory = generated_parameter_dir,
  precomputed_directory = precomputed_parameter_dir,
  required_files = required_parameter_files
)

message(
  "Figure 9(a)--(c) parameter inputs selected from: ",
  parameter_dir
)

figure9_dir <- here::here(
  "Output",
  "Figure9"
)

supplementary_dir <- here::here(
  "Output",
  "Supplementary"
)

dir.create(
  figure9_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  supplementary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

city_config <- list(
  Beijing = list(
    code = "B",
    output_dir = figure9_dir,
    figure_prefix = "Figure9"
  ),
  Shanghai = list(
    code = "S",
    output_dir = supplementary_dir,
    figure_prefix = "FigureS4"
  ),
  HongKong = list(
    code = "H",
    output_dir = supplementary_dir,
    figure_prefix = "FigureS5"
  )
)

parameter_config <- list(
  M = list(
    file_prefix = "mt",
    panel_suffix = "a",
    y_label = expression(M[t]),
    expected_years = 1368:1910,
    height = 4
  ),
  mu = list(
    file_prefix = "mu",
    panel_suffix = "b",
    y_label = expression(mu[t]),
    expected_years = 1368:1911,
    height = 4
  ),
  r2 = list(
    file_prefix = "rt",
    panel_suffix = "c",
    y_label = expression(r[t]^2),
    expected_years = 1368:1911,
    height = 3
  )
)

line_colors <- c(
  "ML" = "steelblue",
  "Penalized ML" = "darkred"
)

line_types <- c(
  "ML" = "dotted",
  "Penalized ML" = "solid"
)

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

read_parameter_file <- function(
    filename,
    city_name,
    parameter_name,
    expected_years) {

  if (!file.exists(
    filename
  )) {
    stop(
      "Missing ",
      parameter_name,
      " file for ",
      city_name,
      ": ",
      filename
    )
  }

  parameter_data <- readr::read_csv(
    filename,
    show_col_types = FALSE
  )

  required_columns <- c(
    "year",
    "coefficient",
    "value"
  )

  missing_columns <- setdiff(
    required_columns,
    names(
      parameter_data
    )
  )

  if (length(
    missing_columns
  ) > 0L) {
    stop(
      filename,
      " is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      "."
    )
  }

  parameter_data <- parameter_data %>%
    transmute(
      year = as.integer(
        year
      ),
      coefficient = as.character(
        coefficient
      ),
      value = as.numeric(
        value
      )
    )

  if (
    anyNA(
      parameter_data
    ) ||
      any(
        !is.finite(
          parameter_data$value
        )
      )
  ) {
    stop(
      "Non-finite or missing entries were found in ",
      filename,
      "."
    )
  }

  expected_coefficient_labels <- c(
    "ML",
    "Penalized ML"
  )

  unexpected_labels <- setdiff(
    unique(
      parameter_data$coefficient
    ),
    expected_coefficient_labels
  )

  if (length(
    unexpected_labels
  ) > 0L) {
    stop(
      "Unexpected coefficient labels in ",
      filename,
      ": ",
      paste(
        unexpected_labels,
        collapse = ", "
      ),
      "."
    )
  }

  duplicate_rows <- parameter_data %>%
    count(
      year,
      coefficient,
      name = "number_of_rows"
    ) %>%
    filter(
      number_of_rows != 1L
    )

  if (nrow(
    duplicate_rows
  ) > 0L) {
    stop(
      "Each year-coefficient combination must occur exactly ",
      "once in ",
      filename,
      "."
    )
  }

  observed_years_by_coefficient <- parameter_data %>%
    group_by(
      coefficient
    ) %>%
    summarise(
      year_set = list(
        sort(
          unique(
            year
          )
        )
      ),
      .groups = "drop"
    )

  if (
    nrow(
      observed_years_by_coefficient
    ) != 2L ||
      any(
        !vapply(
          observed_years_by_coefficient$year_set,
          identical,
          logical(1),
          as.integer(
            expected_years
          )
        )
      )
  ) {
    stop(
      filename,
      " does not contain the expected years for both ML and ",
      "Penalized ML."
    )
  }

  parameter_data %>%
    mutate(
      coefficient = factor(
        coefficient,
        levels = c(
          "Penalized ML",
          "ML"
        )
      ),
      city = city_name,
      parameter = parameter_name
    ) %>%
    arrange(
      coefficient,
      year
    )
}

make_parameter_plot <- function(
    data,
    y_label) {

  ggplot(
    data = data,
    aes(
      x = year,
      y = value,
      colour = coefficient,
      linetype = coefficient
    )
  ) +
    geom_line(
      linewidth = 0.5
    ) +
    scale_colour_manual(
      values = line_colors,
      breaks = c(
        "Penalized ML",
        "ML"
      )
    ) +
    scale_linetype_manual(
      values = line_types,
      breaks = c(
        "Penalized ML",
        "ML"
      )
    ) +
    labs(
      x = "Year",
      y = y_label,
      colour = NULL,
      linetype = NULL
    ) +
    theme(
      text = element_text(
        size = 12
      ),
      plot.title = element_text(
        hjust = 0.5
      ),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

save_parameter_plot <- function(
    plot_object,
    output_file,
    width,
    height) {

  ggsave(
    filename = output_file,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )

  message(
    "Saved: ",
    output_file
  )
}

# ------------------------------------------------------------
# 3. Generate all nine plots
# ------------------------------------------------------------

figure9_plots <- vector(
  "list",
  length(
    city_config
  )
)

names(
  figure9_plots
) <- names(
  city_config
)

all_parameter_plot_data <- list()

for (
  city_name in names(
    city_config
  )
) {

  city_settings <- city_config[[city_name]]

  city_code <- city_settings$code

  figure9_plots[[city_name]] <- list()

  for (
    parameter_name in names(
      parameter_config
    )
  ) {

    parameter_settings <- parameter_config[[parameter_name]]

    parameter_file <- file.path(
      parameter_dir,
      paste0(
        parameter_settings$file_prefix,
        city_code,
        ".csv"
      )
    )

    parameter_data <- read_parameter_file(
      filename = parameter_file,
      city_name = city_name,
      parameter_name = parameter_name,
      expected_years =
        parameter_settings$expected_years
    )

    all_parameter_plot_data[[paste(city_name, parameter_name, sep = "_")]] <- parameter_data

    parameter_plot <- make_parameter_plot(
      data = parameter_data,
      y_label = parameter_settings$y_label
    )

    figure9_plots[[city_name]][[parameter_settings$panel_suffix]] <- parameter_plot

    output_file <- file.path(
      city_settings$output_dir,
      paste0(
        city_settings$figure_prefix,
        parameter_settings$panel_suffix,
        ".png"
      )
    )

    save_parameter_plot(
      plot_object = parameter_plot,
      output_file = output_file,
      width = 6,
      height = parameter_settings$height
    )
  }
}

combined_parameter_plot_data <- bind_rows(
  all_parameter_plot_data
)

readr::write_csv(
  combined_parameter_plot_data,
  file.path(
    figure9_dir,
    "Figure9abc_parameter_plot_data.csv"
  )
)

message(
  "Completed Figure 9, Figure S4, and Figure S5."
)
