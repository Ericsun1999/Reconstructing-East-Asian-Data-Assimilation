here::i_am("Code/Figure7-8.R")

# ============================================================
# Quantile-mapping calibration for Figures 7 and 8
#
# Supported input modes:
#   1. "generated":
#      Read outputs created by prepare_lme_annual.R and
#      Get_tempe_all_data.R.
#
#   2. "precomputed":
#      Read equivalent precomputed files distributed under Data/.
#
#   3. "auto" (default):
#      Use the complete generated input set when available;
#      otherwise use the complete precomputed input set.
#
# Outputs:
#   Figure 7(a)--(d): Beijing
#   Figure 8(a): Beijing
#   Figure 8(b): Beijing
#   Figure 8(c): Shanghai
#   Figure 8(d): Hong Kong
#   Figure7-8_quantile_mapping_values.csv
#
# Figure 7(e) is generated separately by Code/Figure7e.R.
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(np)

# ------------------------------------------------------------
# 1. Select one internally consistent input set
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

generated_files <- c(
  reaches_mean = here::here(
    "Output",
    "Intermediate",
    "REACHES",
    "reaches_kriging_city3_mean.csv"
  ),
  reaches_sd = here::here(
    "Output",
    "Intermediate",
    "REACHES",
    "reaches_kriging_city3_sd.csv"
  ),
  lme_city3 = here::here(
    "Output",
    "Intermediate",
    "LME",
    "lme_city3_annual_1368_1911.csv"
  )
)

precomputed_files <- c(
  reaches_mean = here::here(
    "Data",
    "REACHES",
    "precomputed",
    "reaches_kriging_city3_mean.csv"
  ),
  reaches_sd = here::here(
    "Data",
    "REACHES",
    "precomputed",
    "reaches_kriging_city3_sd.csv"
  ),
  lme_city3 = here::here(
    "Data",
    "LME data",
    "precomputed",
    "lme_city3_annual_1368_1911.csv"
  )
)

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
        "The generated Figure 7--8 input set is incomplete. ",
        "Missing:\n  ",
        paste(
          missing_files,
          collapse = "\n  "
        )
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
        "The precomputed Figure 7--8 input set is incomplete. ",
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

  stop(
    "Neither a complete generated nor a complete precomputed ",
    "Figure 7--8 input set was found.\n",
    "Run prepare_lme_annual.R and Get_tempe_all_data.R, or ",
    "provide all three precomputed files under Data/."
  )
}

input_files <- select_input_set(
  input_mode = input_mode,
  generated_files = generated_files,
  precomputed_files = precomputed_files
)

message(
  "Figure 7--8 inputs selected from: ",
  dirname(
    unname(
      input_files[["reaches_mean"]]
    )
  )
)

# ------------------------------------------------------------
# 2. City definitions
# ------------------------------------------------------------

city_locations <- tibble::tibble(
  city = c(
    "HongKong",
    "Shanghai",
    "Beijing"
  ),
  long = c(
    114.167,
    121.433,
    116.283
  ),
  lat = c(
    22.333,
    31.167,
    39.933
  )
)


coordinate_key <- function(
    long,
    lat) {

  sprintf(
    "%.2f_%.2f",
    as.numeric(
      long
    ),
    as.numeric(
      lat
    )
  )
}

normalize_city_name <- function(city) {

  city <- gsub(
    "[[:space:]_-]",
    "",
    as.character(
      city
    )
  )

  city <- tolower(
    city
  )

  dplyr::recode(
    city,
    hongkong = "HongKong",
    shanghai = "Shanghai",
    beijing = "Beijing",
    .default = NA_character_
  )
}

city_locations <- city_locations %>%
  mutate(
    coordinate_key = coordinate_key(
      long,
      lat
    )
  )

# ------------------------------------------------------------
# 3. Read the kriged REACHES city files
# ------------------------------------------------------------

read_reaches_city_wide <- function(
    input_file,
    value_name,
    city_locations) {

  input_data <- readr::read_csv(
    input_file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  if (
    !"lat" %in%
      names(
        input_data
      ) &&
      "lati" %in%
        names(
          input_data
        )
  ) {
    input_data <- input_data %>%
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
      input_data
    )
  )

  if (length(missing_coordinate_columns) > 0L) {
    stop(
      "The file ",
      input_file,
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
      input_data
    ),
    value = TRUE
  )

  if (length(year_columns) == 0L) {
    stop(
      "No annual columns were found in: ",
      input_file
    )
  }

  long_data <- input_data %>%
    transmute(
      long = as.numeric(
        long
      ),
      lat = as.numeric(
        lat
      ),
      across(
        all_of(
          year_columns
        ),
        as.numeric
      )
    ) %>%
    mutate(
      coordinate_key = coordinate_key(
        long,
        lat
      )
    ) %>%
    left_join(
      city_locations %>%
        dplyr::select(
          city,
          coordinate_key
        ),
      by = "coordinate_key"
    ) %>%
    pivot_longer(
      cols = all_of(
        year_columns
      ),
      names_to = "year_column",
      values_to = value_name
    ) %>%
    mutate(
      year = as.integer(
        sub(
          "^[Xx]",
          "",
          year_column
        )
      )
    ) %>%
    dplyr::select(
      city,
      long,
      lat,
      year,
      all_of(
        value_name
      )
    ) %>%
    arrange(
      city,
      year
    )

  if (anyNA(long_data$city)) {
    stop(
      "At least one REACHES city row could not be matched ",
      "to Hong Kong, Shanghai, or Beijing by coordinates."
    )
  }

  duplicate_city_years <- long_data %>%
    count(
      city,
      year,
      name = "number_of_rows"
    ) %>%
    filter(
      number_of_rows != 1L
    )

  if (nrow(duplicate_city_years) > 0L) {
    stop(
      "Each city-year must occur exactly once in: ",
      input_file
    )
  }

  long_data
}

reaches_mean <- read_reaches_city_wide(
  input_file = unname(
    input_files[["reaches_mean"]]
  ),
  value_name = "reaches_index_mean",
  city_locations = city_locations
)

reaches_sd <- read_reaches_city_wide(
  input_file = unname(
    input_files[["reaches_sd"]]
  ),
  value_name = "reaches_index_sd",
  city_locations = city_locations
)

reaches_city_data <- reaches_mean %>%
  inner_join(
    reaches_sd %>%
      dplyr::select(
        city,
        year,
        reaches_index_sd
      ),
    by = c(
      "city",
      "year"
    )
  ) %>%
  arrange(
    city,
    year
  )

if (
  nrow(reaches_city_data) !=
    nrow(reaches_mean) ||
    nrow(reaches_city_data) !=
      nrow(reaches_sd)
) {
  stop(
    "The REACHES mean and standard-deviation files do not ",
    "contain the same city-year combinations."
  )
}

if (
  any(!is.finite(
    reaches_city_data$reaches_index_mean
  )) ||
    any(!is.finite(
      reaches_city_data$reaches_index_sd
    )) ||
    any(
      reaches_city_data$reaches_index_sd <
        0
    )
) {
  stop(
    "The REACHES city files contain invalid means or ",
    "standard deviations."
  )
}

# ------------------------------------------------------------
# 4. Read the long-format LME city file
# ------------------------------------------------------------

lme_city_data <- readr::read_csv(
  unname(
    input_files[["lme_city3"]]
  ),
  show_col_types = FALSE
)

required_lme_columns <- c(
  "city",
  "long",
  "lati",
  "member",
  "year",
  "temperature_kelvin"
)

missing_lme_columns <- setdiff(
  required_lme_columns,
  names(
    lme_city_data
  )
)

if (length(missing_lme_columns) > 0L) {
  stop(
    "The LME city file is missing columns: ",
    paste(
      missing_lme_columns,
      collapse = ", "
    )
  )
}

lme_city_data <- lme_city_data %>%
  transmute(
    city = normalize_city_name(
      city
    ),
    long = as.numeric(
      long
    ),
    lat = as.numeric(
      lati
    ),
    member = as.character(
      member
    ),
    year = as.integer(
      year
    ),
    temperature_kelvin = as.numeric(
      temperature_kelvin
    ),
    temperature_celsius =
      temperature_kelvin -
      273.15
  ) %>%
  filter(
    year >= 1368L,
    year <= 1911L
  ) %>%
  arrange(
    city,
    member,
    year
  )

if (
  anyNA(
    lme_city_data
  ) ||
    any(!is.finite(
      lme_city_data$temperature_celsius
    ))
) {
  stop(
    "The LME city file contains invalid or missing values."
  )
}

duplicate_lme_rows <- lme_city_data %>%
  count(
    city,
    member,
    year,
    name = "number_of_rows"
  ) %>%
  filter(
    number_of_rows != 1L
  )

if (nrow(duplicate_lme_rows) > 0L) {
  stop(
    "Each LME city-member-year combination must occur ",
    "exactly once."
  )
}

member_counts <- lme_city_data %>%
  distinct(
    city,
    member
  ) %>%
  count(
    city,
    name = "number_of_members"
  )

if (
  nrow(member_counts) != 3L ||
    any(
      member_counts$number_of_members !=
        13L
    )
) {
  stop(
    "Each city must contain exactly 13 LME members."
  )
}

lme_ensemble_mean <- lme_city_data %>%
  group_by(
    city,
    year
  ) %>%
  summarise(
    lme_mean_celsius = mean(
      temperature_celsius
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 5. Quantile-mapping functions
# ------------------------------------------------------------

Fx_hat <- function(
    q,
    fx_bw) {

  fitted_distribution <- np::npudist(
    bws = fx_bw,
    edat = data.frame(
      x = as.numeric(
        q
      )
    )
  )

  as.numeric(
    fitted(
      fitted_distribution
    )
  )
}

FY_hat <- function(
    y,
    yhat,
    nu) {

  yhat <- as.numeric(
    yhat
  )

  nu <- as.numeric(
    nu
  )

  if (length(yhat) != length(nu)) {
    stop(
      "yhat and nu must have the same length."
    )
  }

  nu_safe <- pmax(
    nu,
    1e-8
  )

  vapply(
    as.numeric(
      y
    ),
    function(current_y) {
      mean(
        pnorm(
          (
            current_y -
              yhat
          ) /
            nu_safe
        )
      )
    },
    numeric(1)
  )
}

Fx_inv <- function(
    u,
    x,
    fx_bw) {

  x <- as.numeric(
    x
  )

  u <- pmin(
    pmax(
      as.numeric(
        u
      ),
      1e-8
    ),
    1 -
      1e-8
  )

  x_sd <- stats::sd(
    x
  )

  if (
    !is.finite(
      x_sd
    ) ||
      x_sd <= 0
  ) {
    stop(
      "The LME calibration sample has zero or invalid ",
      "standard deviation."
    )
  }

  lower_bound <- min(
    x
  ) -
    8 *
      x_sd

  upper_bound <- max(
    x
  ) +
    8 *
      x_sd

  lower_probability <- Fx_hat(
    lower_bound,
    fx_bw
  )

  upper_probability <- Fx_hat(
    upper_bound,
    fx_bw
  )

  if (
    lower_probability >
      min(
        u
      ) ||
      upper_probability <
        max(
          u
        )
  ) {
    stop(
      "The numerical interval does not bracket all requested ",
      "LME distribution probabilities."
    )
  }

  vapply(
    u,
    function(current_u) {

      root_function <- function(q) {
        Fx_hat(
          q,
          fx_bw
        ) -
          current_u
      }

      uniroot(
        root_function,
        interval = c(
          lower_bound,
          upper_bound
        ),
        tol = 1e-6
      )$root
    },
    numeric(1)
  )
}

qmapping <- function(
    yhat,
    std,
    x,
    xhat,
    lower_xhat = NULL,
    upper_xhat = NULL,
    ymax = 1,
    ymin = -2) {

  yhat <- as.numeric(
    yhat
  )

  std <- as.numeric(
    std
  )

  x <- as.numeric(
    x
  )

  xhat <- as.numeric(
    xhat
  )

  fx_bw <- np::npudistbw(
    dat = data.frame(
      x = x
    )
  )

  transform_index_to_temperature <- function(index_value) {

    index_probability <- FY_hat(
      index_value,
      yhat,
      std
    )

    Fx_inv(
      u = index_probability,
      x = x,
      fx_bw = fx_bw
    )
  }

  y_seq <- seq(
    ymin,
    ymax,
    length.out = 150
  )

  mapped_curve <- transform_index_to_temperature(
    y_seq
  )

  mapped_xhat <- transform_index_to_temperature(
    xhat
  )

  mapped_lower <- NULL
  mapped_upper <- NULL

  if (
    !is.null(
      lower_xhat
    ) &&
      !is.null(
        upper_xhat
      )
  ) {
    mapped_lower_raw <-
      transform_index_to_temperature(
        lower_xhat
      )

    mapped_upper_raw <-
      transform_index_to_temperature(
        upper_xhat
      )

    mapped_lower <- pmin(
      mapped_lower_raw,
      mapped_upper_raw
    )

    mapped_upper <- pmax(
      mapped_lower_raw,
      mapped_upper_raw
    )
  }

  list(
    y_seq = y_seq,
    mapped_curve = mapped_curve,
    mapped_xhat = mapped_xhat,
    mapped_lower = mapped_lower,
    mapped_upper = mapped_upper,
    fx_bandwidth = fx_bw
  )
}

# ------------------------------------------------------------
# 6. Prepare one city
# ------------------------------------------------------------

prepare_city <- function(
    city_name,
    reaches_city_data,
    lme_city_data,
    lme_ensemble_mean) {

  reaches_data <- reaches_city_data %>%
    filter(
      city == city_name
    ) %>%
    arrange(
      year
    )

  lme_member_data <- lme_city_data %>%
    filter(
      city == city_name
    ) %>%
    arrange(
      member,
      year
    )

  lme_mean_data <- lme_ensemble_mean %>%
    filter(
      city == city_name
    ) %>%
    arrange(
      year
    )

  if (
    nrow(reaches_data) == 0L ||
      nrow(lme_member_data) == 0L
  ) {
    stop(
      "Missing REACHES or LME data for ",
      city_name,
      "."
    )
  }

  quantile_mapping_result <- qmapping(
    yhat =
      reaches_data$reaches_index_mean,
    std =
      reaches_data$reaches_index_sd,
    x =
      lme_member_data$temperature_celsius,
    xhat =
      reaches_data$reaches_index_mean,
    lower_xhat =
      reaches_data$reaches_index_mean -
      reaches_data$reaches_index_sd,
    upper_xhat =
      reaches_data$reaches_index_mean +
      reaches_data$reaches_index_sd,
    ymax = 1,
    ymin = -2
  )

  mapped_reaches <- reaches_data %>%
    mutate(
      temperature_celsius =
        quantile_mapping_result$mapped_xhat,
      lower_temperature_celsius =
        quantile_mapping_result$mapped_lower,
      upper_temperature_celsius =
        quantile_mapping_result$mapped_upper
    )

  comparison_data <- mapped_reaches %>%
    dplyr::select(
      city,
      year,
      reaches_temperature_celsius =
        temperature_celsius
    ) %>%
    inner_join(
      lme_mean_data %>%
        dplyr::select(
          city,
          year,
          lme_temperature_celsius =
            lme_mean_celsius
        ),
      by = c(
        "city",
        "year"
      )
    ) %>%
    arrange(
      year
    )

  if (
    nrow(comparison_data) !=
      nrow(mapped_reaches)
  ) {
    stop(
      "At least one REACHES event year for ",
      city_name,
      " could not be matched to an LME year."
    )
  }

  list(
    city_name = city_name,
    reaches_data = reaches_data,
    lme_member_data = lme_member_data,
    lme_mean_data = lme_mean_data,
    quantile_mapping =
      quantile_mapping_result,
    mapped_reaches = mapped_reaches,
    comparison_data = comparison_data
  )
}

city_names <- c(
  "Beijing",
  "Shanghai",
  "HongKong"
)

city_results <- setNames(
  lapply(
    city_names,
    prepare_city,
    reaches_city_data =
      reaches_city_data,
    lme_city_data =
      lme_city_data,
    lme_ensemble_mean =
      lme_ensemble_mean
  ),
  city_names
)

beijing <- city_results[["Beijing"]]
shanghai <- city_results[["Shanghai"]]
hongkong <- city_results[["HongKong"]]

# ------------------------------------------------------------
# 7. Output directory and graphics helper
# ------------------------------------------------------------

output_dir <- here::here(
  "Output",
  "Figure7-8"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

save_png <- function(
    plot_object,
    filename,
    width,
    height,
    res = 300) {

  output_file <- file.path(
    output_dir,
    filename
  )

  ggsave(
    filename = output_file,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = res
  )

  invisible(
    output_file
  )
}

# Save exact numerical values used in Figures 7--8.
quantile_mapping_output <- bind_rows(
  lapply(
    city_results,
    function(city_result) {

      city_result$mapped_reaches %>%
        dplyr::select(
          city,
          year,
          reaches_index_mean,
          reaches_index_sd,
          mapped_temperature_celsius =
            temperature_celsius,
          mapped_lower_temperature_celsius =
            lower_temperature_celsius,
          mapped_upper_temperature_celsius =
            upper_temperature_celsius
        ) %>%
        left_join(
          city_result$lme_mean_data %>%
            dplyr::select(
              city,
              year,
              lme_ensemble_mean_celsius =
                lme_mean_celsius
            ),
          by = c(
            "city",
            "year"
          )
        )
    }
  )
)

readr::write_csv(
  quantile_mapping_output,
  file.path(
    output_dir,
    "Figure7-8_quantile_mapping_values.csv"
  )
)

# ------------------------------------------------------------
# 8. Figure 7(a)--(d): Beijing
# ------------------------------------------------------------

figure7a_data <- beijing$reaches_data

save_png(
  plot_object =
    ggplot(
      figure7a_data,
      aes(
        x = year,
        y = reaches_index_mean
      )
    ) +
    geom_line(
      colour = "darkorchid1"
    ) +
    labs(
      x = "Year",
      y = "Level"
    ) +
    theme(
      text = element_text(
        size = 19
      ),
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5
      )
    ),
  filename = "Figure7a.png",
  width = 6,
  height = 4
)

figure7b_data <- data.frame(
  index = beijing$quantile_mapping$y_seq,
  temperature_celsius =
    beijing$quantile_mapping$mapped_curve
)

save_png(
  plot_object =
    ggplot(
      figure7b_data,
      aes(
        x = index,
        y = temperature_celsius
      )
    ) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE
    ) +
    labs(
      x = "Index",
      y = expression(
        "Temperature (" *
          degree *
          "C)"
      )
    ) +
    theme(
      text = element_text(
        size = 19
      ),
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5
      )
    ),
  filename = "Figure7b.png",
  width = 6,
  height = 4
)

save_png(
  plot_object =
    ggplot(
      figure7a_data,
      aes(
        x = reaches_index_mean
      )
    ) +
    geom_histogram(
      aes(
        y = after_stat(
          density
        )
      ),
      breaks = seq(
        -2.5,
        1.5,
        by = 0.2
      ),
      fill = "darkorchid1",
      colour = "white"
    ) +
    geom_density(
      colour = "black",
      linewidth = 0.4
    ) +
    labs(
      x = "Index",
      y = "Density"
    ) +
    theme(
      text = element_text(
        size = 18
      ),
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5
      )
    ),
  filename = "Figure7c.png",
  width = 6,
  height = 3.5
)

figure7d_data <- beijing$lme_member_data

save_png(
  plot_object =
    ggplot(
      figure7d_data,
      aes(
        x = temperature_celsius
      )
    ) +
    geom_histogram(
      aes(
        y = after_stat(
          density
        )
      ),
      breaks = seq(
        5,
        9,
        by = 0.1
      ),
      colour = "white",
      fill = "deepskyblue"
    ) +
    geom_density(
      colour = "black",
      linewidth = 0.4
    ) +
    coord_cartesian(
      xlim = c(
        5,
        9
      )
    ) +
    labs(
      x = expression(
        "Temperature (" *
          degree *
          "C)"
      ),
      y = "Density"
    ) +
    theme(
      text = element_text(
        size = 18
      ),
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5
      )
    ),
  filename = "Figure7d.png",
  width = 6,
  height = 3.5
)

# ------------------------------------------------------------
# 9. Figure 8(a): Beijing mapped reconstruction
# ------------------------------------------------------------

figure8a_data <- beijing$mapped_reaches

save_png(
  plot_object =
    ggplot(
      figure8a_data,
      aes(
        x = year,
        y = temperature_celsius
      )
    ) +
    geom_ribbon(
      aes(
        ymin =
          lower_temperature_celsius,
        ymax =
          upper_temperature_celsius
      ),
      alpha = 0.5,
      fill = "grey30"
    ) +
    geom_line(
      colour = "red"
    ) +
    labs(
      x = "Year",
      y = expression(
        "Temperature (" *
          degree *
          "C)"
      )
    ) +
    theme(
      text = element_text(
        size = 11
      ),
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5
      )
    ),
  filename = "Figure8a.png",
  width = 6,
  height = 3
)

# ------------------------------------------------------------
# 10. Figure 8(b)--(d): mapped REACHES and LME ensemble mean
# ------------------------------------------------------------

make_figure8_city_plot <- function(
    city_result) {

  comparison_data <-
    city_result$comparison_data

  ggplot() +
    geom_line(
      data = comparison_data,
      aes(
        x = year,
        y =
          reaches_temperature_celsius
      ),
      colour = "palevioletred1",
      alpha = 0.75
    ) +
    geom_smooth(
      data = comparison_data,
      aes(
        x = year,
        y =
          reaches_temperature_celsius
      ),
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 1.1,
      colour = "firebrick",
      linetype = "dashed",
      span = 0.75
    ) +
    geom_line(
      data = comparison_data,
      aes(
        x = year,
        y =
          lme_temperature_celsius
      ),
      colour = "skyblue1",
      alpha = 0.75
    ) +
    geom_smooth(
      data = comparison_data,
      aes(
        x = year,
        y =
          lme_temperature_celsius
      ),
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 1.1,
      colour = "deepskyblue",
      linetype = "dashed"
    ) +
    labs(
      x = "Year",
      y = expression(
        "Temperature (" *
          degree *
          "C)"
      )
    ) +
    theme(
      text = element_text(
        size = 12
      ),
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5
      )
    )
}

save_png(
  plot_object =
    make_figure8_city_plot(
      beijing
    ),
  filename = "Figure8b.png",
  width = 6,
  height = 3
)

save_png(
  plot_object =
    make_figure8_city_plot(
      shanghai
    ),
  filename = "Figure8c.png",
  width = 6,
  height = 3
)

save_png(
  plot_object =
    make_figure8_city_plot(
      hongkong
    ),
  filename = "Figure8d.png",
  width = 6,
  height = 3
)

message(
  "All Figure 7--8 panels and numerical outputs were saved to: ",
  output_dir
)
