here::i_am("Code/Figure9abc.R")

# ============================================================
# Generate the parameter plots for Beijing, Shanghai, and
# Hong Kong.
#
# Inputs:
#   Data/par/mtB.csv
#   Data/par/muB.csv
#   Data/par/rtB.csv
#   Data/par/mtS.csv
#   Data/par/muS.csv
#   Data/par/rtS.csv
#   Data/par/mtH.csv
#   Data/par/muH.csv
#   Data/par/rtH.csv
#
# Outputs:
#   Beijing:
#     Output/Figure9/Figure9a.png
#     Output/Figure9/Figure9b.png
#     Output/Figure9/Figure9c.png
#
#   Shanghai:
#     Output/Supplementary/FigureS3a.png
#     Output/Supplementary/FigureS3b.png
#     Output/Supplementary/FigureS3c.png
#
#   Hong Kong:
#     Output/Supplementary/FigureS4a.png
#     Output/Supplementary/FigureS4b.png
#     Output/Supplementary/FigureS4c.png
# ============================================================

library(here)
library(ggplot2)

# ------------------------------------------------------------
# 1. Paths and city-specific figure names
# ------------------------------------------------------------

parameter_dir <- here::here(
  "Data",
  "par"
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
    figure_prefix = "FigureS3"
  ),
  HongKong = list(
    code = "H",
    output_dir = supplementary_dir,
    figure_prefix = "FigureS4"
  )
)

# Preserve the original line appearance.
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
    parameter_name) {

  if (!file.exists(filename)) {
    stop(
      "Missing ",
      parameter_name,
      " file for ",
      city_name,
      ": ",
      filename
    )
  }

  data <- read.csv(
    filename,
    check.names = FALSE
  )

  required_columns <- c(
    "year",
    "coefficient",
    "value"
  )

  if (!all(required_columns %in% names(data))) {
    stop(
      filename,
      " must contain the columns: ",
      paste(required_columns, collapse = ", "),
      "."
    )
  }

  data$year <- as.numeric(
    data$year
  )

  data$value <- as.numeric(
    data$value
  )

  data$coefficient <- factor(
    data$coefficient,
    levels = c(
      "ML",
      "Penalized ML"
    )
  )

  if (
    any(!is.finite(data$year)) ||
      any(!is.finite(data$value))
  ) {
    stop(
      "Non-finite year or value entries were found in ",
      filename,
      "."
    )
  }

  if (anyNA(data$coefficient)) {
    stop(
      "Unexpected coefficient labels were found in ",
      filename,
      ". Expected only 'ML' and 'Penalized ML'."
    )
  }

  data
}


make_parameter_plot <- function(
    data,
    y_label) {

  ggplot(
    data = data,
    aes(
      x = year,
      y = value,
      color = coefficient,
      linetype = coefficient
    )
  ) +
    geom_line() +
    scale_color_manual(
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
    xlab("year") +
    ylab(y_label) +
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
  length(city_config)
)

names(figure9_plots) <- names(
  city_config
)

for (city_name in names(city_config)) {

  config <- city_config[[city_name]]
  city_code <- config$code
  figure_prefix <- config$figure_prefix
  output_dir <- config$output_dir

  df_M <- read_parameter_file(
    filename = file.path(
      parameter_dir,
      paste0(
        "mt",
        city_code,
        ".csv"
      )
    ),
    city_name = city_name,
    parameter_name = "M"
  )

  df_mu <- read_parameter_file(
    filename = file.path(
      parameter_dir,
      paste0(
        "mu",
        city_code,
        ".csv"
      )
    ),
    city_name = city_name,
    parameter_name = "mu"
  )

  df_r2 <- read_parameter_file(
    filename = file.path(
      parameter_dir,
      paste0(
        "rt",
        city_code,
        ".csv"
      )
    ),
    city_name = city_name,
    parameter_name = "r2"
  )

  plot_M <- make_parameter_plot(
    data = df_M,
    y_label = expression(M)
  )

  plot_mu <- make_parameter_plot(
    data = df_mu,
    y_label = expression(mu)
  )

  plot_r2 <- make_parameter_plot(
    data = df_r2,
    y_label = expression(r^2)
  )

  figure9_plots[[city_name]] <- list(
    a = plot_M,
    b = plot_mu,
    c = plot_r2
  )

  save_parameter_plot(
    plot_object = plot_M,
    output_file = file.path(
      output_dir,
      paste0(
        figure_prefix,
        "a.png"
      )
    ),
    width = 6,
    height = 4
  )

  save_parameter_plot(
    plot_object = plot_mu,
    output_file = file.path(
      output_dir,
      paste0(
        figure_prefix,
        "b.png"
      )
    ),
    width = 6,
    height = 4
  )

  save_parameter_plot(
    plot_object = plot_r2,
    output_file = file.path(
      output_dir,
      paste0(
        figure_prefix,
        "c.png"
      )
    ),
    width = 6,
    height = 3
  )
}

message(
  "Completed Figure 9, Figure S3, and Figure S4."
)
