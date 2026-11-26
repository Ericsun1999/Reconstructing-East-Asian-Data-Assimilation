here::i_am("Code/Figure2.R")

# ============================================================
# REACHES descriptive plots for Figure 2
#
# This script produces:
#   1. Figure2(a).png: annual counts of REACHES records
#   2. Figure2(b).png: empirical frequencies of temperature levels
#   3. Figure2(c).png: spatial distribution of REACHES records
#
# Outputs are saved under:
#   Output/Figure2/
# ============================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

if (!requireNamespace(
  "mapproj",
  quietly = TRUE
)) {
  stop(
    "The package 'mapproj' is required for coord_map(). ",
    "Install it with renv::install('mapproj')."
  )
}

# ------------------------------------------------------------
# 1. Load REACHES data
# ------------------------------------------------------------

temperature_input_mode <- Sys.getenv(
  "REACHES_TEMPERATURE_INPUT", unset = "original"
)
if (!temperature_input_mode %in% c("original", "augmented")) {
  stop("REACHES_TEMPERATURE_INPUT must be original or augmented.")
}
if (temperature_input_mode == "augmented") {
  prepared_file <- here::here(
    "Output", "Intermediate", "REACHES",
    "temperature_site_year_augmented.csv"
  )
  if (!file.exists(prepared_file)) {
    stop("Prepare the inferred-zero site-years before Figure2.R.")
  }
  temperature <- readr::read_csv(
    prepared_file, show_col_types = FALSE
  )
  if (!all(c("level", "year", "long", "lat",
             "is_inferred_zero") %in% names(temperature))) {
    stop("The inferred-zero dataset is incomplete.")
  }
  temperature <- temperature %>%
    transmute(
      level = as.numeric(level), year = as.integer(year),
      long = as.numeric(long), lat = as.numeric(lat),
      is_inferred_zero = as.logical(is_inferred_zero)
    )
} else {
  temperature <- read_excel(
    here::here("Data", "temperature index value.v1.xlsx"),
    col_types = c(
      "skip", "skip", "numeric", "numeric",
      "skip", "skip", "skip", "skip",
      "skip", "numeric", "numeric",
      "skip", "skip"
    )
  )
  colnames(temperature) <- c("level", "year", "long", "lat")
  temperature$is_inferred_zero <- FALSE
}
# Main reconstruction and descriptive displays use 1368–1911.
temperature <- temperature %>%
  filter(between(year, 1368L, 1911L),
         level %in% c(-2, -1, 0, 1),
         is.finite(long), is.finite(lat))
if (anyNA(temperature$is_inferred_zero)) {
  stop("Missing inferred-zero provenance in Figure 2 data.")
}

# ------------------------------------------------------------
# 2. Figure 2(b): empirical frequencies of temperature levels
# ------------------------------------------------------------

p_figure2b <- ggplot(
  temperature,
  aes(x = level)
) +
  geom_histogram(
    bins = 9
  ) +
  theme(
    text = element_text(size = 15)
  )

# ------------------------------------------------------------
# 3. Figure 2(a): annual counts of REACHES records
# ------------------------------------------------------------

year_count <- temperature %>%
  count(year)

p_figure2a <- ggplot(
  year_count,
  aes(x = year, y = n)
) +
  geom_col(
    fill = "red",
    color = "red"
  ) +
  geom_smooth(
    method = "loess",
    span = 0.3,
    color = "black",
    linewidth = 1,
    se = FALSE
  ) +
  theme(
    text = element_text(size = 15)
  ) +
  labs(
    y = "count"
  )

# ------------------------------------------------------------
# 4. Figure 2(c): spatial distribution of REACHES records
# ------------------------------------------------------------

# Preserve the grouped data structure used in the original code.
temp2 <- temperature %>%
  group_by(year, long, lat)

plot_originREACHES_map <- function(
    data,
    text_size = 15,
    legend_height_cm = 2.5) {

  ggplot(
    data,
    aes(long, lat)
  ) +
    borders(
      database = "world",
      xlim = c(95, 126),
      ylim = c(19, 45),
      fill = NA,
      colour = "grey30"
    ) +
    # Draw inferred zeros first so their outlines cannot obscure
    # documentary temperature observations from other years.
    geom_point(
      data = function(d) d[d$is_inferred_zero, , drop = FALSE],
      shape = 21, fill = "grey94", colour = "grey55",
      size = 1.6, stroke = 0.25, alpha = 0.75,
      show.legend = FALSE
    ) +
    geom_point(
      data = function(d) d[!d$is_inferred_zero & d$level != 0, ,
                           drop = FALSE],
      aes(colour = level),
      size = 1, alpha = 0.9
    ) +
    # Observed index-zero records must remain distinguishable from
    # precipitation-inferred zeros: a dark outline marks true zeros.
    geom_point(
      data = function(d) d[!d$is_inferred_zero & d$level == 0, ,
                           drop = FALSE],
      shape = 21, fill = "white", colour = "grey15",
      stroke = 0.4, size = 1.3, show.legend = FALSE
    ) +
    coord_map(
      xlim = c(98, 124.5),
      ylim = c(19, 42.5)
    ) +
    scale_colour_gradientn(
      colours = rev(
        brewer.pal(
          n = 9,
          name = "RdBu"
        )
      ),
      limits = c(-2, 2),
      na.value = "transparent",
      guide = "colourbar"
    ) +
    theme(
      text = element_text(size = text_size),
      legend.position = c(1.12, 0.61),
      legend.title = element_blank(),
      legend.key.height = grid::unit(
        legend_height_cm,
        "cm"
      )
    )
}

p_figure2c <- plot_originREACHES_map(temp2)

# ------------------------------------------------------------
# 5. Save Figure 2 panels
# ------------------------------------------------------------

figure2_output_dir <- here::here(
  "Output",
  "Figure2"
)

dir.create(
  figure2_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

ggsave(
  filename = file.path(
    figure2_output_dir,
    "Figure2(a).png"
  ),
  plot = p_figure2a,
  width = 4,
  height = 2.5,
  units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(
    figure2_output_dir,
    "Figure2(b).png"
  ),
  plot = p_figure2b,
  width = 4,
  height = 2.5,
  units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(
    figure2_output_dir,
    "Figure2(c).png"
  ),
  plot = p_figure2c,
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)

message(
  "All Figure 2 panels were saved to: ",
  figure2_output_dir
)

# Display figures when running interactively in RStudio.
if (interactive()) {
  print(p_figure2a)
  print(p_figure2b)
  print(p_figure2c)
}
