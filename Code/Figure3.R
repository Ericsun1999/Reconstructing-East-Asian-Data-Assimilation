here::i_am("Code/Figure3.R")

# ============================================================
# LME yearly temperature time series for Figure 3
#
# This script automatically generates:
#   Figure3(a).png: Beijing
#   Figure3(b).png: Shanghai
#   Figure3(c).png: Hong Kong
#
# Outputs are saved under:
#   Output/Figure3/
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# 1. Define city-specific inputs and panel labels
# ------------------------------------------------------------

city_config <- list(
  Beijing = list(
    input_file = here::here("Data", "LME data", "d3.csv"),
    panel_label = "a"
  ),
  Shanghai = list(
    input_file = here::here("Data", "LME data", "d2.csv"),
    panel_label = "b"
  ),
  HongKong = list(
    input_file = here::here("Data", "LME data", "d1.csv"),
    panel_label = "c"
  )
)

figure3_output_dir <- here::here("Output", "Figure3")

dir.create(
  figure3_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Generate one Figure 3 panel
# ------------------------------------------------------------

make_figure3 <- function(
    city_name,
    input_file,
    panel_label,
    output_dir = figure3_output_dir) {

  if (!file.exists(input_file)) {
    stop(
      "Input file for ",
      city_name,
      " was not found: ",
      input_file
    )
  }

  city_data <- read.csv(
    input_file,
    row.names = 1,
    check.names = FALSE
  )

  if (ncol(city_data) < 3L) {
    stop(
      "The input data for ",
      city_name,
      " must contain at least two metadata columns ",
      "followed by yearly temperature columns."
    )
  }

  # The first two columns contain location information.
  temperature_data <- city_data[
    ,
    -c(1, 2),
    drop = FALSE
  ]

  if (!all(vapply(temperature_data, is.numeric, logical(1)))) {
    stop(
      "One or more yearly temperature columns for ",
      city_name,
      " are not numeric."
    )
  }

  # Rows correspond to 13 LME simulations.
  # Columns correspond to years 1350--1949.
  temperature_matrix <- as.matrix(temperature_data) - 273

  expected_years <- 1350:1949

  if (nrow(temperature_matrix) != 13L) {
    stop(
      "Expected 13 LME simulations for ",
      city_name,
      ", but found ",
      nrow(temperature_matrix),
      "."
    )
  }

  if (ncol(temperature_matrix) != length(expected_years)) {
    stop(
      "Expected ",
      length(expected_years),
      " yearly values for ",
      city_name,
      " covering 1350--1949, but found ",
      ncol(temperature_matrix),
      "."
    )
  }

  # Convert the simulation matrix to long format.
  simulation_data <- data.frame(
    Year = rep(
      expected_years,
      times = nrow(temperature_matrix)
    ),
    Simulation = factor(
      rep(
        seq_len(nrow(temperature_matrix)),
        each = length(expected_years)
      )
    ),
    Temperature = as.vector(t(temperature_matrix))
  )

  average_data <- data.frame(
    Year = expected_years,
    Temperature = colMeans(
      temperature_matrix,
      na.rm = TRUE
    )
  )

  # ----------------------------------------------------------
  # Create the city-specific plot
  # ----------------------------------------------------------

  p_figure3 <- ggplot() +
    geom_line(
      data = simulation_data,
      aes(
        x = Year,
        y = Temperature,
        group = Simulation
      ),
      colour = "grey50",
      linewidth = 0.35
    ) +
    geom_line(
      data = average_data,
      aes(
        x = Year,
        y = Temperature
      ),
      colour = "deepskyblue",
      linewidth = 0.8
    ) +
    xlim(1368, 1911) +
    labs(
      x = "year",
      y = "temperature"
    ) +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )

  output_file <- file.path(
    output_dir,
    paste0("Figure3(", panel_label, ").png")
  )

  ggsave(
    filename = output_file,
    plot = p_figure3,
    width = 6,
    height = 3,
    units = "in",
    dpi = 300
  )

  message(
    "Figure 3 panel ",
    panel_label,
    " for ",
    city_name,
    " saved to: ",
    output_file
  )

  invisible(output_file)
}

# ------------------------------------------------------------
# 3. Generate Figure 3(a)--(c)
# ------------------------------------------------------------

figure3_output_files <- vapply(
  names(city_config),
  function(city_name) {
    config <- city_config[[city_name]]

    make_figure3(
      city_name = city_name,
      input_file = config$input_file,
      panel_label = config$panel_label
    )
  },
  character(1)
)

message(
  "All Figure 3 panels were saved to: ",
  figure3_output_dir
)
