here::i_am("Code/Supplementary/FigureS1.R")

# ============================================================
# Figure S1(a)--(e):
# Maps of temperature event levels for different periods
#
# Outputs:
#   Output/Supplementary/FigureS1a.png
#   Output/Supplementary/FigureS1b.png
#   Output/Supplementary/FigureS1c.png
#   Output/Supplementary/FigureS1d.png
#   Output/Supplementary/FigureS1e.png
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
temperature <- read_excel(
  here::here("Data", "temperature index value.v1.xlsx"),
  col_types = c(
    "skip", "skip", "numeric", "numeric",
    "skip", "skip", "skip", "skip", "skip",
    "numeric", "numeric", "skip", "skip"
  )
)

colnames(temperature) <- c("level", "year", "long", "lat")

# Keep the same structure as the original code
temp2 <- temperature %>%
  group_by(year, long, lat)

# ------------------------------------------------------------
# 2. Split into the five periods used in Figure S1
# ------------------------------------------------------------
period_data <- list(
  FigureS1a = temp2 %>% filter(year < 1501),
  FigureS1b = temp2 %>% filter(year > 1500 & year < 1601),
  FigureS1c = temp2 %>% filter(year > 1600 & year < 1701),
  FigureS1d = temp2 %>% filter(year > 1700 & year < 1801),
  FigureS1e = temp2 %>% filter(year > 1800 & year <= 1911)
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
# 5. Generate and save Figure S1(a)--(e)
# ------------------------------------------------------------
figureS1_plots <- vector("list", length(period_data))
names(figureS1_plots) <- names(period_data)

for (figure_name in names(period_data)) {
  plot_object <- plot_originREACHES_map(period_data[[figure_name]])

  figureS1_plots[[figure_name]] <- plot_object

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

message("Completed Figure S1(a)--(e).")
