here::i_am("Code/Figure9abc.R")

# ============================================================
# Generate Figure 9(a)--(c) automatically for Beijing,
# Shanghai, and Hong Kong.
#
# Required inputs produced by Code/Prior.R:
#   Output/Intermediate/Prior/mtB.csv
#   Output/Intermediate/Prior/muB.csv
#   Output/Intermediate/Prior/rtB.csv
#   Output/Intermediate/Prior/mtS.csv
#   Output/Intermediate/Prior/muS.csv
#   Output/Intermediate/Prior/rtS.csv
#   Output/Intermediate/Prior/mtH.csv
#   Output/Intermediate/Prior/muH.csv
#   Output/Intermediate/Prior/rtH.csv
#
# Outputs:
#   Output/Figure9/Beijing/Figure9a.png
#   Output/Figure9/Beijing/Figure9b.png
#   Output/Figure9/Beijing/Figure9c.png
#   Output/Figure9/Shanghai/Figure9a.png
#   ...
# ============================================================

library(here)
library(ggplot2)

# ------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------

prior_dir <- here::here(
  "Output",
  "Intermediate",
  "Prior"
)

figure_dir <- here::here(
  "Output",
  "Figure9"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

city_config <- list(
  Beijing = list(code = "B"),
  Shanghai = list(code = "S"),
  HongKong = list(code = "H")
)

# Preserve the appearance of the original plots:
#   ML            = blue dotted line
#   Penalized ML  = dark-red solid line
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
      filename,
      "\nRun Code/Prior.R first."
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
    width = 6,
    height = 4) {

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
# 3. Generate Figure 9(a)--(c) for every city
# ------------------------------------------------------------

figure9_plots <- vector(
  "list",
  length(city_config)
)

names(figure9_plots) <- names(
  city_config
)

for (city_name in names(city_config)) {

  city_code <- city_config[[city_name]]$code

  city_output_dir <- file.path(
    figure_dir,
    city_name
  )

  dir.create(
    city_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  df_M <- read_parameter_file(
    filename = file.path(
      prior_dir,
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
      prior_dir,
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
      prior_dir,
      paste0(
        "rt",
        city_code,
        ".csv"
      )
    ),
    city_name = city_name,
    parameter_name = "r2"
  )

  # Figure 9(a): estimated M_t.
  plot_M <- make_parameter_plot(
    data = df_M,
    y_label = expression(M)
  )

  # Figure 9(b): estimated mu_t.
  plot_mu <- make_parameter_plot(
    data = df_mu,
    y_label = expression(mu)
  )

  # Figure 9(c): estimated r_t^2.
  plot_r2 <- make_parameter_plot(
    data = df_r2,
    y_label = expression(r^2)
  )

  figure9_plots[[city_name]] <- list(
    Figure9a = plot_M,
    Figure9b = plot_mu,
    Figure9c = plot_r2
  )

  save_parameter_plot(
    plot_object = plot_M,
    output_file = file.path(
      city_output_dir,
      "Figure9a.png"
    ),
    width = 6,
    height = 4
  )

  save_parameter_plot(
    plot_object = plot_mu,
    output_file = file.path(
      city_output_dir,
      "Figure9b.png"
    ),
    width = 6,
    height = 4
  )

  save_parameter_plot(
    plot_object = plot_r2,
    output_file = file.path(
      city_output_dir,
      "Figure9c.png"
    ),
    width = 6,
    height = 3
  )
}

message(
  "Figure 9(a)--(c) completed for: ",
  paste(
    names(city_config),
    collapse = ", "
  )
)
