here::i_am("Code/Supplementary/FigureS2.R")

# ============================================================
# Generate Figure S2(a)--(c):
#   Figure S2(a): Beijing
#   Figure S2(b): Shanghai
#   Figure S2(c): Hong Kong
#
# Required inputs:
#   Data/GHCNv4.xlsx
#   Data/Valid/tempBv5.csv
#   Data/Valid/tempSv5.csv
#   Data/Valid/tempHv5.csv
#
# Outputs:
#   Output/Supplementary/FigureS2a.png
#   Output/Supplementary/FigureS2b.png
#   Output/Supplementary/FigureS2c.png
# ============================================================

library(here)
library(ggplot2)
library(readxl)
library(dplyr)

ghcn_file <- here::here("Data", "GHCNv4.xlsx")
validation_dir <- here::here("Data", "Valid")
output_dir <- here::here("Output", "Supplementary")

if (!file.exists(ghcn_file)) {
  stop("GHCN file was not found: ", ghcn_file)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
    "Expected one year column and 12 monthly columns, but found ",
    ncol(ghcn_all),
    " columns."
  )
}

names(ghcn_all) <- c("year", paste0("month", 1:12))

city_config <- list(
  Beijing = list(
    row_indices = 1:161,
    remove_local_rows = c(103, 104),
    start_year_exclusive = 1837,
    validation_file = file.path(validation_dir, "tempBv5.csv"),
    output_file = file.path(output_dir, "FigureS2a.png")
  ),
  Shanghai = list(
    row_indices = 443:503,
    remove_local_rows = integer(0),
    start_year_exclusive = 1847,
    validation_file = file.path(validation_dir, "tempSv5.csv"),
    output_file = file.path(output_dir, "FigureS2b.png")
  ),
  HongKong = list(
    row_indices = 839:895,
    remove_local_rows = integer(0),
    start_year_exclusive = 1854,
    validation_file = file.path(validation_dir, "tempHv5.csv"),
    output_file = file.path(output_dir, "FigureS2c.png")
  )
)

series_colors <- c(
  "Assimilated" = "black",
  "REACHES" = "firebrick",
  "LME" = "deepskyblue",
  "GHCN" = "green3"
)

series_linetypes <- c(
  "Assimilated" = "solid",
  "REACHES" = "dashed",
  "LME" = "dotted"
)

legend_order <- c("Assimilated", "GHCN", "LME", "REACHES")

prepare_ghcn_city <- function(
    ghcn_data,
    row_indices,
    remove_local_rows,
    city_name) {

  if (max(row_indices) > nrow(ghcn_data)) {
    stop(
      "The configured GHCN row range for ",
      city_name,
      " exceeds the workbook size."
    )
  }

  city_data <- ghcn_data[row_indices, , drop = FALSE]

  if (length(remove_local_rows) > 0L) {
    if (max(remove_local_rows) > nrow(city_data)) {
      stop(
        "A configured row removal for ",
        city_name,
        " exceeds the selected block."
      )
    }

    city_data <- city_data[-remove_local_rows, , drop = FALSE]
  }

  monthly_values <- as.matrix(
    city_data[, paste0("month", 1:12), drop = FALSE]
  )

  storage.mode(monthly_values) <- "double"

  complete_year <- rowSums(is.na(monthly_values)) == 0L

  annual_temperature <- rowMeans(
    monthly_values,
    na.rm = TRUE
  ) / 100

  result <- data.frame(
    year = as.numeric(city_data$year[complete_year]),
    GHCN = annual_temperature[complete_year]
  )

  if (
    any(!is.finite(result$year)) ||
      any(!is.finite(result$GHCN))
  ) {
    stop(
      "Non-finite GHCN annual values were produced for ",
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
    "sigmasmooth2",
    "alpha",
    "beta",
    "LME"
  )

  if (!all(required_columns %in% names(data))) {
    stop(
      filename,
      " must contain: ",
      paste(required_columns, collapse = ", "),
      "."
    )
  }

  data[, required_columns] <- lapply(
    data[, required_columns, drop = FALSE],
    as.numeric
  )

  if (
    any(!is.finite(data$year)) ||
      any(!is.finite(data$predicted)) ||
      any(!is.finite(data$sigmasmooth2)) ||
      any(!is.finite(data$alpha)) ||
      any(!is.finite(data$beta)) ||
      any(!is.finite(data$LME))
  ) {
    stop("Non-finite values were found in ", filename, ".")
  }

  if (any(data$sigmasmooth2 < 0)) {
    stop("Negative smoothing variances were found in ", filename, ".")
  }

  data
}


prepare_plot_data <- function(
    validation_data,
    ghcn_data,
    start_year_exclusive) {

  # Assimilated, REACHES, LME, and the uncertainty ribbon are
  # available annually, so retain every available validation
  # year after the city-specific starting year.
  validation_data <- validation_data %>%
    mutate(
      REACHES_transformed = predicted * beta + alpha,
      se = sqrt(sigmasmooth2),
      lo = predicted - 1.96 * se,
      hi = predicted + 1.96 * se
    ) %>%
    filter(
      year > start_year_exclusive
    )

  # GHCN is plotted only for years with all 12 monthly values.
  # Do not restrict the annual model series to these GHCN years.
  ghcn_data <- ghcn_data %>%
    filter(
      year > start_year_exclusive
    )

  list(
    validation = validation_data,
    ghcn = ghcn_data
  )
}


make_time_series_plot <- function(
    validation_data,
    ghcn_data) {

  ggplot() +
    geom_ribbon(
      data = validation_data,
      aes(
        x = year,
        ymin = lo,
        ymax = hi
      ),
      fill = "grey30",
      alpha = 0.25
    ) +
    geom_line(
      data = validation_data,
      aes(
        x = year,
        y = predicted,
        color = "Assimilated",
        linetype = "Assimilated"
      ),
      linewidth = 0.7,
      alpha = 0.75
    ) +
    geom_line(
      data = validation_data,
      aes(
        x = year,
        y = REACHES_transformed,
        color = "REACHES",
        linetype = "REACHES"
      ),
      linewidth = 0.7,
      alpha = 0.75
    ) +
    geom_line(
      data = validation_data,
      aes(
        x = year,
        y = LME,
        color = "LME",
        linetype = "LME"
      ),
      linewidth = 0.7,
      alpha = 0.75
    ) +
    geom_point(
      data = ghcn_data,
      aes(
        x = year,
        y = GHCN,
        color = "GHCN"
      ),
      size = 0.7
    ) +
    scale_color_manual(
      values = series_colors,
      breaks = legend_order
    ) +
    scale_linetype_manual(
      values = series_linetypes
    ) +
    xlab("year") +
    ylab("temperature") +
    theme(
      text = element_text(size = 12),
      legend.position = "bottom",
      legend.key.height = grid::unit(-0.5, "cm"),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(
      linetype = "none",
      color = guide_legend(
        title = NULL,
        override.aes = list(
          linetype = c(
            "solid",
            "blank",
            "dotted",
            "dashed"
          ),
          shape = c(
            NA,
            16,
            NA,
            NA
          )
        )
      )
    )
}

figureS2_plots <- vector(
  "list",
  length(city_config)
)

names(figureS2_plots) <- names(city_config)

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
    start_year_exclusive = config$start_year_exclusive
  )

  if (nrow(plot_data$validation) == 0L) {
    stop(
      "No annual model results after ",
      config$start_year_exclusive,
      " were found for ",
      city_name,
      "."
    )
  }

  if (nrow(plot_data$ghcn) == 0L) {
    warning(
      "No complete GHCN years after ",
      config$start_year_exclusive,
      " were found for ",
      city_name,
      ". The model series will still be plotted."
    )
  }

  plot_object <- make_time_series_plot(
    validation_data = plot_data$validation,
    ghcn_data = plot_data$ghcn
  )

  figureS2_plots[[city_name]] <- plot_object

  ggsave(
    filename = config$output_file,
    plot = plot_object,
    width = 6,
    height = 3,
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

message("Completed Figure S2(a)--(c).")
