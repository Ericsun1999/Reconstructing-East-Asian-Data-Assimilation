here::i_am("Code/Supplementary/FigureS2.R")

# ============================================================
# Figure S2(a)--(e):
# Maps of temperature event levels for different periods
#
# Outputs:
#   Output/Supplementary/FigureS2a.png
#   Output/Supplementary/FigureS2b.png
#   Output/Supplementary/FigureS2c.png
#   Output/Supplementary/FigureS2d.png
#   Output/Supplementary/FigureS2e.png
# ============================================================

library(here)
library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(maps)

# ------------------------------------------------------------
# 1. Read REACHES temperature-event data
# ------------------------------------------------------------
temperature_input_mode <- Sys.getenv(
  "REACHES_TEMPERATURE_INPUT", unset = "original"
)
if (!temperature_input_mode %in% c("original", "augmented")) {
  stop("REACHES_TEMPERATURE_INPUT must be original or augmented.")
}
if (temperature_input_mode == "augmented") {
  input_file <- here::here(
    "Output", "Intermediate", "REACHES",
    "temperature_site_year_augmented.csv"
  )
  temperature <- readr::read_csv(input_file, show_col_types = FALSE)
  if (!all(c("level", "year", "long", "lat", "is_inferred_zero") %in%
           names(temperature))) {
    stop("Incomplete inferred-zero source for Supplementary Figure S2.")
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
      "skip", "skip", "skip", "skip", "skip",
      "numeric", "numeric", "skip", "skip"
    )
  )
  colnames(temperature) <- c("level", "year", "long", "lat")
  temperature$is_inferred_zero <- FALSE
}
temperature <- temperature %>% filter(
  between(year, 1368L, 1911L),
  level %in% c(-2, -1, 0, 1),
  is.finite(long), is.finite(lat)
)

# Keep the same structure as the original code
temp2 <- temperature %>%
  group_by(year, long, lat)

# ------------------------------------------------------------
# 2. Split into the five periods used in Figure S2
# ------------------------------------------------------------
period_data <- list(
  FigureS2a = temp2 %>% filter(year < 1501),
  FigureS2b = temp2 %>% filter(year > 1500 & year < 1601),
  FigureS2c = temp2 %>% filter(year > 1600 & year < 1701),
  FigureS2d = temp2 %>% filter(year > 1700 & year < 1801),
  FigureS2e = temp2 %>% filter(year > 1800 & year <= 1911)
)

# ------------------------------------------------------------
# 3. Plotting function
# ------------------------------------------------------------
plot_originREACHES_map <- function(
    data,
    text_size = 15,
    legend_height_cm = 2.5) {

  ggplot(data, aes(long, lat)) +
    borders(
      database = "world",
      xlim = c(95, 126),
      ylim = c(19, 45),
      fill = NA,
      colour = "grey30"
    ) +
    geom_point(
      aes(colour = level),
      size = 1
    ) +
    geom_point(
      data = function(d) d[d$is_inferred_zero, , drop = FALSE],
      shape = 21, fill = NA, colour = "grey40",
      size = 1, stroke = 0.35, show.legend = FALSE
    ) +
    coord_map(
      xlim = c(98, 124.5),
      ylim = c(19, 42.5)
    ) +
    scale_colour_gradientn(
      colours = rev(brewer.pal(n = 9, name = "RdBu")),
      limits = c(-2, 2),
      na.value = "transparent",
      guide = "colourbar"
    ) +
    theme(
      text = element_text(size = text_size),
      legend.position = c(1.12, 0.61),
      legend.title = element_blank(),
      legend.key.height = grid::unit(legend_height_cm, "cm")
    )
}

# ------------------------------------------------------------
# 4. Output directory
# ------------------------------------------------------------
output_dir <- here::here("Output", "Supplementary")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 5. Generate and save Figure S2(a)--(e)
# ------------------------------------------------------------
figureS2_plots <- vector("list", length(period_data))
names(figureS2_plots) <- names(period_data)

for (figure_name in names(period_data)) {
  plot_object <- plot_originREACHES_map(period_data[[figure_name]])

  figureS2_plots[[figure_name]] <- plot_object

  ggsave(
    filename = file.path(output_dir, paste0(figure_name, ".png")),
    plot = plot_object,
    width = 5,
    height = 4,
    units = "in",
    dpi = 300
  )

  message("Saved: ", file.path(output_dir, paste0(figure_name, ".png")))
}

message("Completed Figure S2(a)--(e).")
