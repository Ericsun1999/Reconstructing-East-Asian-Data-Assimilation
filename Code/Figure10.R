here::i_am("Code/Figure10.R")

# ============================================================
# Generate Figure 10(a)--(c):
#   Figure 10(a): Beijing
#   Figure 10(b): Shanghai
#   Figure 10(c): Hong Kong
#
# Required inputs:
#   Data/GHCNv4.xlsx
#   Data/Valid/tempBv5.csv
#   Data/Valid/tempSv5.csv
#   Data/Valid/tempHv5.csv
#
# Outputs:
#   Output/Figure10/Figure10a.png
#   Output/Figure10/Figure10b.png
#   Output/Figure10/Figure10c.png
# ============================================================

library(here)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------


ghcn_file <- here::here("Data", "GHCNv4.xlsx")

validation_dir <- here::here(
  "Data",
  "Valid"
)

output_dir <- here::here(
  "Output",
  "Figure10"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Read the GHCN workbook
#
# The workbook contains repeating value/flag columns. This
# preserves the original selection of year plus 12 monthly
# temperature columns.
# ------------------------------------------------------------

ghcn_column_types <- c(
  "skip", "numeric", "skip",
  "numeric", "skip", "skip",
  "numeric", "skip", "skip",
  "numeric", "skip", "skip",
  "numeric", "skip", "skip",
  "numeric", "skip", "skip",
  "numeric", "skip", "skip",
  "numeric", "skip", "skip",
  "numeric", "skip", "skip",
  "numeric", "skip", "skip",
  "numeric", "skip", "skip",
  "numeric", "skip", "skip",
  "numeric", "skip", "skip"
)

ghcn_all <- readxl::read_excel(
  ghcn_file,
  col_types = ghcn_column_types
)

if (ncol(ghcn_all) != 13L) {
  stop(
    "Expected the selected GHCN data to contain one year ",
    "column and 12 monthly columns, but found ",
    ncol(ghcn_all),
    " columns."
  )
}

names(ghcn_all) <- c(
  "year",
  paste0(
    "month",
    1:12
  )
)

# ------------------------------------------------------------
# 3. City configuration
# ------------------------------------------------------------

city_config <- list(
  Beijing = list(
    row_indices = 1:161,
    remove_local_rows = c(103, 104),
    validation_file = file.path(
      validation_dir,
      "tempBv5.csv"
    ),
    output_file = file.path(
      output_dir,
      "Figure10a.png"
    ),
    x_limits = c(10, 13),
    y_limits = c(10, 13),
    label_x = 10.05,
    label_y = c(12.90, 12.65, 12.40)
  ),
  Shanghai = list(
    row_indices = 443:503,
    remove_local_rows = integer(0),
    validation_file = file.path(
      validation_dir,
      "tempSv5.csv"
    ),
    output_file = file.path(
      output_dir,
      "Figure10b.png"
    ),
    x_limits = c(15, 18),
    y_limits = c(15, 18),
    label_x = 15.05,
    label_y = c(17.88, 17.63, 17.38)
  ),
  HongKong = list(
    row_indices = 839:895,
    remove_local_rows = integer(0),
    validation_file = file.path(
      validation_dir,
      "tempHv5.csv"
    ),
    output_file = file.path(
      output_dir,
      "Figure10c.png"
    ),
    x_limits = c(21.5, 24.5),
    y_limits = c(21.5, 24.5),
    label_x = 21.58,
    label_y = c(24.35, 24.10, 23.85)
  )
)

method_levels <- c(
  "Assimilated",
  "LME",
  "REACHES"
)

method_colors <- c(
  "Assimilated" = "#F8766D",
  "LME" = "#00BA38",
  "REACHES" = "#619CFF"
)

method_shapes <- c(
  "Assimilated" = 16,
  "LME" = 17,
  "REACHES" = 15
)

base_text_size <- 23
correlation_text_size <- 7
legend_text_size <- 18

# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------

prepare_ghcn_city <- function(
    ghcn_data,
    row_indices,
    remove_local_rows,
    city_name) {

  if (max(row_indices) > nrow(ghcn_data)) {
    stop(
      "The configured GHCN row range for ",
      city_name,
      " exceeds the number of rows in GHCNv4.xlsx."
    )
  }

  city_data <- ghcn_data[
    row_indices,
    ,
    drop = FALSE
  ]

  if (length(remove_local_rows) > 0L) {

    if (max(remove_local_rows) > nrow(city_data)) {
      stop(
        "A configured row removal for ",
        city_name,
        " exceeds the selected city block."
      )
    }

    city_data <- city_data[
      -remove_local_rows,
      ,
      drop = FALSE
    ]
  }

  monthly_values <- as.matrix(
    city_data[
      ,
      paste0(
        "month",
        1:12
      ),
      drop = FALSE
    ]
  )

  storage.mode(monthly_values) <- "double"

  # Preserve the original rule: retain only years with all 12
  # monthly temperatures available.
  complete_months <- rowSums(
    is.na(monthly_values)
  ) == 0L

  annual_temperature <- rowMeans(
    monthly_values,
    na.rm = TRUE
  ) /
    100

  result <- data.frame(
    year = as.numeric(
      city_data$year[
        complete_months
      ]
    ),
    GHCN = annual_temperature[
      complete_months
    ]
  )

  if (
    any(!is.finite(result$year)) ||
      any(!is.finite(result$GHCN))
  ) {
    stop(
      "Non-finite annual GHCN values were produced for ",
      city_name,
      "."
    )
  }

  result
}


read_validation_data <- function(
    filename,
    city_name) {

  if (!file.exists(filename)) {
    stop(
      "Validation file for ",
      city_name,
      " was not found: ",
      filename,
      "\nRun Code/Figure9d.R first."
    )
  }

  data <- read.csv(
    filename,
    check.names = FALSE
  )

  required_columns <- c(
    "year",
    "predicted",
    "REACHES",
    "LME"
  )

  if (!all(required_columns %in% names(data))) {
    stop(
      filename,
      " must contain: ",
      paste(
        required_columns,
        collapse = ", "
      ),
      "."
    )
  }

  data$year <- as.numeric(
    data$year
  )

  data
}


prepare_plot_data <- function(
    validation_data,
    ghcn_data,
    city_name) {

  joined_data <- dplyr::inner_join(
    validation_data,
    ghcn_data,
    by = "year"
  )

  if (nrow(joined_data) == 0L) {
    stop(
      "No common years were found between GHCN and the ",
      "validation data for ",
      city_name,
      "."
    )
  }

  long_data <- joined_data |>
    dplyr::select(
      year,
      GHCN,
      predicted,
      REACHES,
      LME
    ) |>
    tidyr::pivot_longer(
      cols = c(
        predicted,
        LME,
        REACHES
      ),
      names_to = "Method",
      values_to = "Estimate"
    ) |>
    dplyr::mutate(
      Method = dplyr::recode(
        Method,
        predicted = "Assimilated",
        LME = "LME",
        REACHES = "REACHES"
      ),
      Method = factor(
        Method,
        levels = method_levels
      )
    )

  long_data
}


make_correlation_labels <- function(
    plot_data,
    label_x,
    label_y) {

  correlation_data <- plot_data |>
    dplyr::group_by(
      Method
    ) |>
    dplyr::summarise(
      cor_val = cor(
        GHCN,
        Estimate,
        use = "complete.obs"
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      Method
    )

  # Reapply the intended order explicitly.
  correlation_data$Method <- factor(
    correlation_data$Method,
    levels = method_levels
  )

  correlation_data <- correlation_data |>
    dplyr::arrange(
      Method
    ) |>
    dplyr::mutate(
      label = paste0(
        Method,
        ": Cor = ",
        round(
          cor_val,
          2
        )
      ),
      x = label_x,
      y = label_y
    )

  correlation_data
}


make_scatter_plot <- function(
    plot_data,
    correlation_labels,
    x_limits,
    y_limits) {

  ggplot(
    plot_data,
    aes(
      x = GHCN,
      y = Estimate,
      color = Method,
      shape = Method
    )
  ) +
    geom_point(
      size = 3.5,
      alpha = 0.85
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "red",
      linewidth = 0.8
    ) +
    geom_text(
      data = correlation_labels,
      aes(
        x = x,
        y = y,
        label = label,
        color = Method
      ),
      inherit.aes = FALSE,
      show.legend = FALSE,
      hjust = 0,
      size = correlation_text_size
    ) +
    scale_color_manual(
      values = method_colors,
      breaks = method_levels
    ) +
    scale_shape_manual(
      values = method_shapes,
      breaks = method_levels
    ) +
    coord_equal(
      xlim = x_limits,
      ylim = y_limits
    ) +
    labs(
      x = "GHCN",
      y = "Estimate"
    ) +
    theme_bw() +
    theme(
      text = element_text(
        size = base_text_size
      ),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(
        size = legend_text_size
      )
    )
}

# ------------------------------------------------------------
# 5. Generate Figure 10(a)--(c)
# ------------------------------------------------------------

figure10_plots <- vector(
  "list",
  length(city_config)
)

names(figure10_plots) <- names(
  city_config
)

correlation_tables <- vector(
  "list",
  length(city_config)
)

names(correlation_tables) <- names(
  city_config
)

for (city_name in names(city_config)) {

  config <- city_config[[city_name]]

  ghcn_city <- prepare_ghcn_city(
    ghcn_data = ghcn_all,
    row_indices = config$row_indices,
    remove_local_rows = config$remove_local_rows,
    city_name = city_name
  )

  validation_data <- read_validation_data(
    filename = config$validation_file,
    city_name = city_name
  )

  plot_data <- prepare_plot_data(
    validation_data = validation_data,
    ghcn_data = ghcn_city,
    city_name = city_name
  )

  correlation_labels <- make_correlation_labels(
    plot_data = plot_data,
    label_x = config$label_x,
    label_y = config$label_y
  )

  plot_object <- make_scatter_plot(
    plot_data = plot_data,
    correlation_labels = correlation_labels,
    x_limits = config$x_limits,
    y_limits = config$y_limits
  )

  figure10_plots[[city_name]] <- plot_object
  correlation_tables[[city_name]] <- correlation_labels

  ggsave(
    filename = config$output_file,
    plot = plot_object,
    width = 6,
    height = 6,
    units = "in",
    dpi = 300
  )

  message(
    "Saved ",
    city_name,
    ": ",
    config$output_file
  )
}

correlation_summary <- dplyr::bind_rows(
  lapply(
    names(correlation_tables),
    function(city_name) {

      data.frame(
        city = city_name,
        method = as.character(
          correlation_tables[[city_name]]$Method
        ),
        correlation = correlation_tables[
          [city_name]
        ]$cor_val
      )
    }
  )
)

write.csv(
  correlation_summary,
  file.path(
    output_dir,
    "Figure10_correlations.csv"
  ),
  row.names = FALSE
)

message(
  "Completed Figure 10(a), Figure 10(b), and Figure 10(c)."
)
