here::i_am("Code/Figure3.R")

# ============================================================
# LME annual temperature time series for Figure 3
#
# Supported input modes:
#   1. "generated":
#      Output/Intermediate/LME/
#        lme_city3_annual_1368_1911.csv
#
#   2. "precomputed":
#      Data/LME data/precomputed/
#        lme_city3_annual_1368_1911.csv
#
#   3. "auto" (default):
#      Use the generated file when available; otherwise use the
#      precomputed copy distributed with the repository.
#
# Outputs:
#   Output/Figure3/Figure3(a).png  Beijing
#   Output/Figure3/Figure3(b).png  Shanghai
#   Output/Figure3/Figure3(c).png  Hong Kong
#   Output/Figure3/Figure3_plot_data.csv
# ============================================================

library(ggplot2)
library(dplyr)
library(readr)

# ------------------------------------------------------------
# 1. Select the input source
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

input_files <- c(
  generated = here::here(
    "Output",
    "Intermediate",
    "LME",
    "lme_city3_annual_1368_1911.csv"
  ),
  precomputed = here::here(
    "Data",
    "LME data",
    "precomputed",
    "lme_city3_annual_1368_1911.csv"
  )
)

select_input_file <- function(
    input_mode,
    input_files) {

  file_available <- file.exists(
    input_files
  )

  if (input_mode == "generated") {
    if (!file_available[["generated"]]) {
      stop(
        "The generated Figure 3 input was not found:\n  ",
        input_files[["generated"]],
        "\nRun Code/DataPreparation/prepare_lme_annual.R first."
      )
    }

    return(
      unname(
        input_files[["generated"]]
      )
    )
  }

  if (input_mode == "precomputed") {
    if (!file_available[["precomputed"]]) {
      stop(
        "The precomputed Figure 3 input was not found:\n  ",
        input_files[["precomputed"]]
      )
    }

    return(
      unname(
        input_files[["precomputed"]]
      )
    )
  }

  # For "auto", prefer a freshly generated file. If it is not
  # available, fall back to the precomputed repository copy.
  if (file_available[["generated"]]) {
    return(
      unname(
        input_files[["generated"]]
      )
    )
  }

  if (file_available[["precomputed"]]) {
    return(
      unname(
        input_files[["precomputed"]]
      )
    )
  }

  stop(
    "Neither Figure 3 input file was found.\n",
    "Generated location:\n  ",
    input_files[["generated"]],
    "\nPrecomputed location:\n  ",
    input_files[["precomputed"]],
    "\nRun Code/DataPreparation/prepare_lme_annual.R or ",
    "add the precomputed file to Data/LME data/precomputed/."
  )
}

figure3_input_file <- select_input_file(
  input_mode = input_mode,
  input_files = input_files
)

message(
  "Figure 3 input selected: ",
  figure3_input_file
)

# ------------------------------------------------------------
# 2. Read and validate the common three-city file
# ------------------------------------------------------------

city_data <- readr::read_csv(
  figure3_input_file,
  show_col_types = FALSE
)

required_columns <- c(
  "city",
  "long",
  "lati",
  "member",
  "year",
  "temperature_kelvin"
)

missing_columns <- setdiff(
  required_columns,
  names(
    city_data
  )
)

if (length(missing_columns) > 0L) {
  stop(
    "The Figure 3 input is missing the following columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

city_data <- city_data %>%
  transmute(
    city = as.character(city),
    long = as.numeric(long),
    lati = as.numeric(lati),
    member = as.character(member),
    year = as.integer(year),
    temperature_kelvin =
      as.numeric(
        temperature_kelvin
      )
  ) %>%
  filter(
    year >= 1368L,
    year <= 1911L
  )

if (
  anyNA(
    city_data[
      ,
      required_columns
    ]
  )
) {
  stop(
    "The Figure 3 input contains missing values in required ",
    "columns."
  )
}

expected_cities <- c(
  "Beijing",
  "Shanghai",
  "HongKong"
)

missing_cities <- setdiff(
  expected_cities,
  unique(
    city_data$city
  )
)

if (length(missing_cities) > 0L) {
  stop(
    "The Figure 3 input is missing the following cities: ",
    paste(
      missing_cities,
      collapse = ", "
    )
  )
}

duplicate_rows <- city_data %>%
  count(
    city,
    member,
    year,
    name = "number_of_rows"
  ) %>%
  filter(
    number_of_rows != 1L
  )

if (nrow(duplicate_rows) > 0L) {
  stop(
    "Each city-member-year combination must occur exactly once."
  )
}

member_counts <- city_data %>%
  distinct(
    city,
    member
  ) %>%
  count(
    city,
    name = "number_of_members"
  )

if (any(
  member_counts$number_of_members != 13L
)) {
  stop(
    "Each city must contain exactly 13 LME ensemble members."
  )
}

year_counts <- city_data %>%
  distinct(
    city,
    year
  ) %>%
  count(
    city,
    name = "number_of_years"
  )

expected_number_of_years <- length(
  1368:1911
)

if (any(
  year_counts$number_of_years !=
    expected_number_of_years
)) {
  stop(
    "Each city must contain all annual values from 1368 to 1911."
  )
}

# Convert Kelvin to degrees Celsius using the exact physical
# offset rather than the previous approximation of 273.
plot_data <- city_data %>%
  mutate(
    temperature_celsius =
      temperature_kelvin -
      273.15,
    member = factor(
      member,
      levels = paste0(
        "a",
        1:13
      )
    )
  ) %>%
  arrange(
    city,
    member,
    year
  )

ensemble_mean_data <- plot_data %>%
  group_by(
    city,
    year
  ) %>%
  summarise(
    temperature_celsius =
      mean(
        temperature_celsius
      ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 3. Define panel order and output directory
# ------------------------------------------------------------

city_config <- list(
  Beijing = list(
    panel_label = "a",
    display_name = "Beijing"
  ),
  Shanghai = list(
    panel_label = "b",
    display_name = "Shanghai"
  ),
  HongKong = list(
    panel_label = "c",
    display_name = "Hong Kong"
  )
)

figure3_output_dir <- here::here(
  "Output",
  "Figure3"
)

dir.create(
  figure3_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save the exact values used in the figure.
figure3_plot_data_file <- file.path(
  figure3_output_dir,
  "Figure3_plot_data.csv"
)

readr::write_csv(
  plot_data %>%
    select(
      city,
      long,
      lati,
      member,
      year,
      temperature_kelvin,
      temperature_celsius
    ),
  figure3_plot_data_file
)

# ------------------------------------------------------------
# 4. Generate one Figure 3 panel
# ------------------------------------------------------------

make_figure3 <- function(
    city_name,
    panel_label,
    display_name,
    plot_data,
    ensemble_mean_data,
    output_dir) {

  simulation_city_data <- plot_data %>%
    filter(
      city == city_name
    )

  mean_city_data <- ensemble_mean_data %>%
    filter(
      city == city_name
    )

  if (nrow(simulation_city_data) == 0L) {
    stop(
      "No LME data were found for ",
      city_name,
      "."
    )
  }

  p_figure3 <- ggplot() +
    geom_line(
      data = simulation_city_data,
      aes(
        x = year,
        y = temperature_celsius,
        group = member
      ),
      colour = "grey50",
      linewidth = 0.35,
      alpha = 0.9
    ) +
    geom_line(
      data = mean_city_data,
      aes(
        x = year,
        y = temperature_celsius
      ),
      colour = "deepskyblue3",
      linewidth = 0.85
    ) +
    scale_x_continuous(
      limits = c(
        1368,
        1911
      ),
      breaks = c(
        1400,
        1500,
        1600,
        1700,
        1800,
        1900
      ),
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    labs(
      title = display_name,
      x = "Year",
      y = expression(
        "Temperature (" *
          degree *
          "C)"
      )
    ) +
    theme_gray(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5
      ),
      legend.position = "none"
    )

  output_file <- file.path(
    output_dir,
    paste0(
      "Figure3(",
      panel_label,
      ").png"
    )
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
    display_name,
    " saved to: ",
    output_file
  )

  invisible(
    output_file
  )
}

# ------------------------------------------------------------
# 5. Generate Figure 3(a)--(c)
# ------------------------------------------------------------

figure3_output_files <- vapply(
  names(
    city_config
  ),
  function(city_name) {

    config <- city_config[[city_name]]

    make_figure3(
      city_name = city_name,
      panel_label =
        config$panel_label,
      display_name =
        config$display_name,
      plot_data = plot_data,
      ensemble_mean_data =
        ensemble_mean_data,
      output_dir =
        figure3_output_dir
    )
  },
  character(1)
)

message(
  "All Figure 3 panels were saved to: ",
  figure3_output_dir
)

message(
  "Figure 3 plotting data were saved to: ",
  figure3_plot_data_file
)
