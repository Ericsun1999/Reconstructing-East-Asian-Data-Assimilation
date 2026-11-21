here::i_am("Code/Analysis/diagnose_precipitation_augmentation.R")

# Diagnostics for the optional documentary-availability temperature input.
# All annual precipitation records, including dry-weather and drought
# reports, remain eligible to support an inferred zero. Event categories
# below are descriptive, non-exclusive, and do not control imputation.
# Input source workbooks and the augmented input are never modified.

library(here)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

augmented_file <- here::here(
  "Output", "Intermediate", "REACHES",
  "temperature_site_year_augmented.csv"
)
temperature_file <- here::here("Data", "temperature index value.v1.xlsx")
precipitation_file <- path.expand(Sys.getenv(
  "REACHES_PRECIPITATION_FILE",
  unset = here::here("Data", "precipitation index value.xlsx")
))
output_dir <- here::here(
  "Output", "Intermediate", "Precipitation", "AugmentationDiagnostics"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (p in c(augmented_file, temperature_file, precipitation_file)) {
  if (!file.exists(p)) stop("Required input file not found: ", p)
}

write_result <- function(data, filename) {
  readr::write_csv(data, file.path(output_dir, filename))
}

assign_block <- function(year) {
  dplyr::case_when(
    year <= 1500L ~ "1368-1500",
    year <= 1600L ~ "1501-1600",
    year <= 1700L ~ "1601-1700",
    year <= 1800L ~ "1701-1800",
    TRUE ~ "1801-1911"
  )
}

mode_first <- function(values) {
  frequencies <- table(values)
  modes <- as.numeric(names(frequencies)[frequencies == max(frequencies)])
  values[which(values %in% modes)[1L]]
}

augmented <- readr::read_csv(
  augmented_file, show_col_types = FALSE
)
required <- c(
  "year", "long", "lat", "level", "source_type",
  "is_inferred_zero", "number_of_precipitation_records"
)
if (!all(required %in% names(augmented))) {
  stop("Missing augmented data columns: ",
       paste(setdiff(required, names(augmented)), collapse = ", "))
}
augmented <- augmented %>%
  mutate(
    year = as.integer(year),
    long = round(as.numeric(long), 8),
    lat = round(as.numeric(lat), 8),
    level = as.integer(level),
    is_inferred_zero = as.logical(is_inferred_zero),
    block = assign_block(year)
  )
key <- c("year", "long", "lat")
if (anyDuplicated(augmented[key])) {
  duplicate_keys <- augmented %>%
    count(across(all_of(key))) %>%
    filter(n > 1L)
  write_result(duplicate_keys, "augmentation_duplicate_coordinates.csv")
  stop(
    "Duplicate augmented site-years after rounding coordinates to eight ",
    "decimal places. Rerun ",
    "Code/Analysis/prepare_precipitation_augmented_reaches.R ",
    "to rebuild the normalized dataset. See ",
    file.path(output_dir, "augmentation_duplicate_coordinates.csv")
  )
}
if (any(is.na(augmented$is_inferred_zero))) {
  stop("Missing inferred-zero flags.")
}
if (any(
  augmented$is_inferred_zero &
    (augmented$level != 0L |
       augmented$source_type != "inferred_zero_precipitation")
)) {
  stop("An inferred temperature record is not a flagged zero.")
}
if (any(
  !augmented$is_inferred_zero &
    augmented$source_type != "observed_temperature"
)) stop("Invalid observed-temperature source flag.")

# Independently reconstruct the original site-year temperature category
# (modal value, with source-row order breaking ties).
temperature_source <- readxl::read_excel(
  temperature_file, sheet = "全年溫度分級", col_types = "text"
) %>%
  transmute(
    year = suppressWarnings(as.integer(.data[["資料判斷後，所屬年份"]])),
    long = round(suppressWarnings(as.numeric(.data[["空間-經度"]])), 8),
    lat = round(suppressWarnings(as.numeric(.data[["空間-緯度"]])), 8),
    level = suppressWarnings(as.integer(.data[["溫度等級"]]))
  ) %>%
  filter(
    between(year, 1368L, 1911L),
    level %in% c(-2L, -1L, 0L, 1L),
    is.finite(long), is.finite(lat)
  ) %>%
  group_by(year, long, lat) %>%
  summarise(level = mode_first(level), .groups = "drop")
recorded <- augmented %>%
  filter(!is_inferred_zero) %>%
  select(all_of(key), level)
recorded_missing_from_source <- anti_join(
  recorded, temperature_source, by = c(key, "level")
)
source_missing_from_recorded <- anti_join(
  temperature_source, recorded, by = c(key, "level")
)
if (nrow(recorded) != nrow(temperature_source) ||
    nrow(recorded_missing_from_source) > 0L ||
    nrow(source_missing_from_recorded) > 0L) {
  write_result(
    recorded_missing_from_source,
    "augmentation_recorded_not_in_source.csv"
  )
  write_result(
    source_missing_from_recorded,
    "augmentation_source_not_in_recorded.csv"
  )
  stop(
    "The recorded temperature site-years differ from the source workbook: ",
    nrow(recorded), " recorded, ",
    nrow(temperature_source), " original, ",
    nrow(recorded_missing_from_source), " unmatched recorded, ",
    nrow(source_missing_from_recorded), " unmatched source. ",
    "Inspect the augmentation_recorded_not_in_source.csv and ",
    "augmentation_source_not_in_recorded.csv diagnostic files."
  )
}

precipitation <- readxl::read_excel(
  precipitation_file, sheet = "全年降水分級", col_types = "text"
) %>%
  transmute(
    year = suppressWarnings(as.integer(.data[["資料判斷後，所屬年份"]])),
    long = round(suppressWarnings(as.numeric(.data[["空間-經度"]])), 8),
    lat = round(suppressWarnings(as.numeric(.data[["空間-緯度"]])), 8),
    event_codes = coalesce(.data[["事件編碼"]], ""),
    description = coalesce(.data[["事件-事件敘述"]], "")
  ) %>%
  filter(
    between(year, 1368L, 1911L),
    is.finite(long), is.finite(lat)
  )

# Prefixes were checked against examples in the source workbook.
# The 10111 family also includes prayers for snow; do not interpret
# this group as proof of snow absence for every record.
tokens <- strsplit(precipitation$event_codes, ";", fixed = TRUE)
has_prefix <- function(codes, prefix) {
  any(startsWith(trimws(codes), prefix))
}
no_rain_code <- vapply(tokens, has_prefix, logical(1L), prefix = "10011")
snow_absence_code <- vapply(
  tokens, has_prefix, logical(1L), prefix = "10111"
)
hail_code <- vapply(tokens, has_prefix, logical(1L), prefix = "1021")
drought_code <- vapply(tokens, has_prefix, logical(1L), prefix = "300")
# Positive and negative codes may coexist within one record, so
# exclusions must be applied token by token, not across the record.
rain_event_code <- vapply(tokens, function(codes) {
  codes <- trimws(codes)
  any(startsWith(codes, "100") & !startsWith(codes, "10011"))
}, logical(1L))
snow_event_code <- vapply(tokens, function(codes) {
  codes <- trimws(codes)
  any(startsWith(codes, "101") & !startsWith(codes, "10111"))
}, logical(1L))
precipitation <- precipitation %>%
  mutate(
    rain_event = rain_event_code,
    no_rain = no_rain_code,
    drought = drought_code | grepl("旱", description),
    snow_event = snow_event_code,
    snow_absence_or_prayer = snow_absence_code,
    hail_event = hail_code
  )

precipitation_site_year <- precipitation %>%
  group_by(year, long, lat) %>%
  summarise(
    n_precipitation_records = n(),
    rain_event = any(rain_event),
    no_rain = any(no_rain),
    drought = any(drought),
    snow_event = any(snow_event),
    snow_absence_or_prayer = any(snow_absence_or_prayer),
    hail_event = any(hail_event),
    .groups = "drop"
  )

inferred <- augmented %>%
  filter(is_inferred_zero) %>%
  left_join(precipitation_site_year, by = key)

if (any(is.na(inferred$n_precipitation_records))) {
  stop("Some inferred zero site-years have no annual precipitation source.")
}
if (any(
  inferred$number_of_precipitation_records !=
    inferred$n_precipitation_records
)) stop("Source precipitation record counts do not match the input.")
if (nrow(anti_join(
  inferred, precipitation_site_year, by = key
)) != 0L) stop("An inferred zero lacks a precipitation source.")

# Temporal summaries, including the full 544-year calendar sequence.
yearly <- augmented %>%
  group_by(year) %>%
  summarise(
    recorded_temperature = sum(!is_inferred_zero),
    recorded_zero = sum(!is_inferred_zero & level == 0L),
    inferred_zero = sum(is_inferred_zero),
    augmented_total = n(),
    .groups = "drop"
  ) %>%
  right_join(tibble(year = 1368:1911), by = "year") %>%
  mutate(
    across(c(recorded_temperature, recorded_zero,
             inferred_zero, augmented_total),
           ~ replace_na(.x, 0L)),
    block = assign_block(year),
    inferred_share = if_else(
      augmented_total > 0L, inferred_zero / augmented_total, NA_real_
    )
  ) %>%
  arrange(year)
write_result(yearly, "augmentation_counts_by_year.csv")

new_years <- yearly %>%
  filter(recorded_temperature == 0L, inferred_zero > 0L)
write_result(new_years, "augmentation_new_informative_years.csv")

by_block <- yearly %>%
  group_by(block) %>%
  summarise(
    recorded_temperature_site_years = sum(recorded_temperature),
    originally_recorded_zero = sum(recorded_zero),
    inferred_zero_site_years = sum(inferred_zero),
    augmented_site_years = sum(augmented_total),
    years_with_recorded_temperature = sum(recorded_temperature > 0),
    years_with_augmented_temperature = sum(augmented_total > 0),
    recorded_zero_proportion = if_else(
      sum(recorded_temperature) > 0,
      sum(recorded_zero) / sum(recorded_temperature), NA_real_
    ),
    augmented_zero_proportion = if_else(
      sum(augmented_total) > 0,
      (sum(recorded_zero) + sum(inferred_zero)) / sum(augmented_total),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  arrange(match(block, c(
    "1368-1500", "1501-1600", "1601-1700",
    "1701-1800", "1801-1911"
  )))
write_result(by_block, "augmentation_counts_by_block.csv")

# Compare the entire ordinal distribution, not only the zero share.
level_distribution <- bind_rows(
  augmented %>%
    filter(!is_inferred_zero) %>%
    count(block, level, name = "site_years") %>%
    mutate(dataset = "recorded_only"),
  augmented %>%
    count(block, level, name = "site_years") %>%
    mutate(dataset = "with_inferred_zero")
) %>%
  group_by(dataset, block) %>%
  mutate(
    total_site_years = sum(site_years),
    proportion = site_years / total_site_years
  ) %>%
  ungroup() %>%
  arrange(block, dataset, level)
write_result(
  level_distribution,
  "augmentation_temperature_level_by_block.csv"
)
total_level_distribution <- bind_rows(
  augmented %>%
    filter(!is_inferred_zero) %>%
    count(level, name = "site_years") %>%
    mutate(dataset = "recorded_only"),
  augmented %>%
    count(level, name = "site_years") %>%
    mutate(dataset = "with_inferred_zero")
) %>%
  group_by(dataset) %>%
  mutate(
    total_site_years = sum(site_years),
    proportion = site_years / total_site_years
  ) %>%
  ungroup()
write_result(
  total_level_distribution,
  "augmentation_temperature_level_overall.csv"
)


# Non-exclusive event categories: a site-year may belong to multiple groups.
category_fields <- c(
  "rain_event", "no_rain", "drought", "snow_event",
  "snow_absence_or_prayer", "hail_event"
)
inferred$other_precipitation_event <- (
  rowSums(as.matrix(inferred[, category_fields])) == 0L
)
category_fields <- c(category_fields, "other_precipitation_event")
categories <- inferred %>%
  select(year, long, lat, block, all_of(category_fields)) %>%
  pivot_longer(
    cols = all_of(category_fields),
    names_to = "event_category", values_to = "category_present"
  ) %>%
  filter(category_present)
category_totals <- categories %>%
  count(event_category, name = "site_years") %>%
  mutate(
    proportion_of_inferred = site_years / nrow(inferred),
    categories_can_overlap = TRUE
  ) %>%
  arrange(desc(site_years))
write_result(category_totals, "augmentation_event_category_counts.csv")
write_result(
  categories %>%
    count(block, event_category, name = "site_years") %>%
    arrange(block, desc(site_years)),
  "augmentation_event_category_by_block.csv"
)

category_patterns <- apply(
  as.matrix(inferred[, category_fields]),
  1L,
  function(values) paste(category_fields[as.logical(values)], collapse = "+")
)
write_result(
  tibble(category_combination = category_patterns) %>%
    count(category_combination, sort = TRUE, name = "site_years"),
  "augmentation_event_category_combinations.csv"
)
write_result(
  inferred %>%
    count(number_of_precipitation_records, name = "site_years") %>%
    arrange(number_of_precipitation_records),
  "augmentation_source_record_counts.csv"
)

# Two-degree site-year density cells; each record remains a separate
# annual observation, rather than a distinct geographical location.
spatial <- augmented %>%
  mutate(
    longitude_cell_center = floor(long / 2) * 2 + 1,
    latitude_cell_center = floor(lat / 2) * 2 + 1
  ) %>%
  count(
    source_type, longitude_cell_center, latitude_cell_center,
    name = "site_years"
  )
write_result(spatial, "augmentation_spatial_2deg_counts.csv")
write_result(
  augmented %>%
    mutate(
      longitude_cell_center = floor(long / 2) * 2 + 1,
      latitude_cell_center = floor(lat / 2) * 2 + 1
    ) %>%
    count(
      block, source_type, longitude_cell_center, latitude_cell_center,
      name = "site_years"
    ),
  "augmentation_spatial_2deg_by_block.csv"
)

p_year <- yearly %>%
  select(year, recorded_temperature, inferred_zero) %>%
  pivot_longer(
    cols = c(recorded_temperature, inferred_zero),
    names_to = "record_type", values_to = "site_years"
  ) %>%
  ggplot(aes(x = year, y = site_years, fill = record_type)) +
  geom_col(width = 0.95) +
  labs(
    title = "Documentary temperature coverage before and after augmentation",
    x = "Year", y = "Site-year count", fill = "Source"
  ) +
  theme_minimal(base_size = 11)
ggsave(
  file.path(output_dir, "augmentation_counts_by_year.png"),
  p_year, width = 12, height = 5, dpi = 180, bg = "white"
)

# Use the same geographic frame as Figure2(c): grey world borders,
# geographic projection, and point markers. Each point represents a
# two-degree cell center, with its color encoding the site-year count.
# A sequential palette is used because counts do not have a natural zero
# midpoint (unlike the signed temperature index shown in Figure2(c)).
if (!requireNamespace("mapproj", quietly = TRUE)) {
  stop(
    "The package 'mapproj' is required for the Figure 2 map projection. ",
    "Install it with renv::install('mapproj')."
  )
}
if (!requireNamespace("maps", quietly = TRUE)) {
  stop(
    "The package 'maps' is required for world borders. ",
    "Install it with renv::install('maps')."
  )
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  stop(
    "The package 'RColorBrewer' is required for the map palette. ",
    "Install it with renv::install('RColorBrewer')."
  )
}

spatial_plot_data <- spatial %>%
  mutate(
    source_type = factor(
      source_type,
      levels = c("observed_temperature", "inferred_zero_precipitation"),
      labels = c("Recorded temperature", "Inferred zero")
    )
  )

p_spatial <- ggplot() +
  borders(
    database = "world",
    xlim = c(95, 126),
    ylim = c(19, 45),
    fill = NA,
    colour = "grey30",
    linewidth = 0.35
  ) +
  geom_point(
    data = spatial_plot_data,
    aes(
      x = longitude_cell_center,
      y = latitude_cell_center,
      colour = site_years
    ),
    size = 3,
    alpha = 0.9
  ) +
  facet_wrap(~ source_type, ncol = 2) +
  coord_map(
    xlim = c(98, 124.5),
    ylim = c(19, 42.5),
    expand = FALSE
  ) +
  scale_colour_gradientn(
    colours = RColorBrewer::brewer.pal(n = 9, name = "Blues"),
    trans = "sqrt",
    name = "Site-years\nper 2-degree cell"
  ) +
  labs(
    title = "Geographic coverage of documentary site-years",
    subtitle = "Aggregated to 2-degree cells, 1368-1911",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "grey95", colour = "grey80"),
    strip.text = element_text(face = "bold"),
    legend.key.height = grid::unit(3, "cm"),
    plot.title = element_text(face = "bold")
  )
ggsave(
  file.path(output_dir, "augmentation_spatial_2deg.png"),
  p_spatial, width = 11, height = 6, dpi = 300, bg = "white"
)

p_levels <- total_level_distribution %>%
  ggplot(aes(
    x = dataset, y = proportion, fill = factor(level)
  )) +
  geom_col(width = 0.6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Temperature-index distribution before and after augmentation",
    x = "Documentary input", y = "Proportion",
    fill = "Temperature index"
  ) +
  theme_minimal(base_size = 11)
ggsave(
  file.path(output_dir, "augmentation_temperature_level_shares.png"),
  p_levels, width = 9, height = 5, dpi = 180, bg = "white"
)

p_events <- category_totals %>%
  ggplot(aes(
    x = reorder(event_category, site_years),
    y = site_years
  )) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Precipitation event categories supporting inferred zeros",
    subtitle = "Non-exclusive event-code categories; totals may overlap",
    x = "Event category", y = "Inferred-zero site-years"
  ) +
  theme_minimal(base_size = 11)
ggsave(
  file.path(output_dir, "augmentation_event_category_counts.png"),
  p_events, width = 9, height = 5, dpi = 180, bg = "white"
)

metadata <- tibble(
  observed_site_years = nrow(recorded),
  inferred_zero_site_years = nrow(inferred),
  augmented_site_years = nrow(augmented),
  years_with_recorded_temperature = sum(yearly$recorded_temperature > 0),
  years_with_augmented_temperature = sum(yearly$augmented_total > 0),
  newly_informative_years = nrow(new_years),
  observed_index_zero_count = sum(recorded$level == 0L),
  inferred_supported_by_one_precipitation_record =
    sum(inferred$number_of_precipitation_records == 1L),
  event_category_counts_are_nonexclusive = TRUE,
  precipitation_source = precipitation_file
)
write_result(metadata, "augmentation_diagnostic_metadata.csv")
print(metadata)
print(by_block)
print(category_totals)
message("Augmentation diagnostics saved to: ", output_dir)
