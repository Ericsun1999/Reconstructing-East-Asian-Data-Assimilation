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

# ------------------------------------------------------------
# 1. Load REACHES data
# ------------------------------------------------------------

temperature <- read_excel(
  here::here(
    "Data",
    "temperature index value.v1.xlsx"
  ),
  col_types = c(
    "skip", "skip", "numeric", "numeric",
    "skip", "skip", "skip", "skip",
    "skip", "numeric", "numeric",
    "skip", "skip"
  )
)

colnames(temperature) <- c(
  "level",
  "year",
  "long",
  "lat"
)

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
    geom_point(
      aes(colour = level),
      cex = 1
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
