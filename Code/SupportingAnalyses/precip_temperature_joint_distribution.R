here::i_am("Code/SupportingAnalyses/precip_temperature_joint_distribution.R")

# Annual REACHES precipitation/temperature co-occurrence analysis.
# All category rules below are provisional until checked against
# the REACHES event-code dictionary. No zero values are imputed.

library(here)
library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

precip_file <- path.expand(Sys.getenv(
  "REACHES_PRECIPITATION_FILE",
  unset = here::here("Data", "precipitation index value.xlsx")
))
temp_file <- here::here("Data", "temperature index value.v1.xlsx")
output_dir <- here::here("Output", "Intermediate", "Precipitation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(precip_file)) stop("Precipitation file not found: ", precip_file)
if (!file.exists(temp_file)) stop("Temperature file not found: ", temp_file)

# The annual worksheets must not be combined with the summer/winter worksheets.
precip <- readxl::read_excel(
  precip_file, sheet = "全年降水分級", col_types = "text"
) %>%
  transmute(
    year = as.integer(.data[["資料判斷後，所屬年份"]]),
    long = as.numeric(.data[["空間-經度"]]),
    lat = as.numeric(.data[["空間-緯度"]]),
    precipitation_index = as.numeric(.data[["降水等級"]]),
    text_id = .data[["資料編號-文本ID"]],
    event_id = .data[["資料編號-事件編碼"]],
    event_codes = coalesce(.data[["事件編碼"]], ""),
    description = .data[["事件-事件敘述"]]
  ) %>%
  filter(between(year, 1368L, 1911L), is.finite(long), is.finite(lat))

temperature <- readxl::read_excel(
  temp_file, sheet = "全年溫度分級", col_types = "text"
) %>%
  transmute(
    year = as.integer(.data[["資料判斷後，所屬年份"]]),
    long = as.numeric(.data[["空間-經度"]]),
    lat = as.numeric(.data[["空間-緯度"]]),
    Z = as.integer(.data[["溫度等級"]]),
    source_order = row_number()
  ) %>%
  filter(
    between(year, 1368L, 1911L),
    Z %in% c(-2L, -1L, 0L, 1L),
    is.finite(long), is.finite(lat)
  )

# Match the existing calibration pipeline: modal index at each site-year;
# if tied, take the first modal category in the original row order.
mode_first <- function(x) {
  frequencies <- table(x)
  modes <- as.integer(names(frequencies)[frequencies == max(frequencies)])
  x[which(x %in% modes)[1L]]
}

temperature_site_year <- temperature %>%
  group_by(year, long, lat) %>%
  summarise(Z = mode_first(Z), n_temperature_records = n(), .groups = "drop")

# Preliminary event-code grouping based on inspection of annual records:
#   100* = rain-related; 10011* includes non-rain/prayer-related codes;
#   101* = snow-related; 10111* includes non-snow codes;
#   1021* = hail-related; 1022*/1023* also include other ice phenomena.
# These are working coding rules, not a substitute for the event-code dictionary.
code_flags <- function(codes) {
  tokens <- trimws(codes)
  c(
    rain = any(startsWith(tokens, "100") & !startsWith(tokens, "10011")),
    snow = any(startsWith(tokens, "101") & !startsWith(tokens, "10111")),
    hail = any(startsWith(tokens, "1021"))
  )
}

token_lists <- strsplit(precip$event_codes, ";", fixed = TRUE)
category_matrix <- t(vapply(token_lists, code_flags, logical(3L)))
precip <- dplyr::bind_cols(precip, as.data.frame(category_matrix))

# One row per location/year. A location/year may have several categories,
# but repeated documentary entries are never counted multiple times per category.
precipitation_site_year <- precip %>%
  group_by(year, long, lat) %>%
  summarise(
    rain = any(rain),
    snow = any(snow),
    hail = any(hail),
    n_precipitation_records = n(),
    source_text_ids = paste(unique(text_id[!is.na(text_id)]), collapse = ";"),
    source_event_ids = paste(unique(event_id[!is.na(event_id)]), collapse = ";"),
    .groups = "drop"
  )

precipitation_categories <- precipitation_site_year %>%
  pivot_longer(c(rain, snow, hail), names_to = "C", values_to = "present") %>%
  filter(present) %>%
  select(-present)

# This is the sample R = 1: precipitation AND temperature at the same site/year.
matched <- precipitation_categories %>%
  inner_join(temperature_site_year, by = c("year", "long", "lat"))

blocks <- c("1368-1500", "1501-1600", "1601-1700", "1701-1800", "1801-1911")
assign_block <- function(year) {
  case_when(
    year <= 1500L ~ blocks[1L],
    year <= 1600L ~ blocks[2L],
    year <= 1700L ~ blocks[3L],
    year <= 1800L ~ blocks[4L],
    TRUE ~ blocks[5L]
  )
}

# P(Z = k | C = c, R = 1) and the denominator N_c for every category.
conditional <- matched %>%
  count(C, Z, name = "n") %>%
  complete(C = c("rain", "snow", "hail"),
           Z = c(-2L, -1L, 0L, 1L), fill = list(n = 0L)) %>%
  group_by(C) %>%
  mutate(N = sum(n), probability = if_else(N > 0, n / N, NA_real_)) %>%
  ungroup() %>%
  arrange(match(C, c("rain", "snow", "hail")), Z)

conditional_by_block <- matched %>%
  mutate(block = assign_block(year)) %>%
  count(block, C, Z, name = "n") %>%
  complete(
    block = blocks, C = c("rain", "snow", "hail"),
    Z = c(-2L, -1L, 0L, 1L), fill = list(n = 0L)
  ) %>%
  group_by(block, C) %>%
  mutate(N = sum(n), probability = if_else(N > 0, n / N, NA_real_)) %>%
  ungroup() %>%
  arrange(match(block, blocks), match(C, c("rain", "snow", "hail")), Z)

# R = 0 candidate site-years for later review; do NOT impute zero here.
candidates <- precipitation_categories %>%
  anti_join(temperature_site_year, by = c("year", "long", "lat")) %>%
  mutate(block = assign_block(year))

sample_sizes <- precipitation_categories %>%
  left_join(temperature_site_year, by = c("year", "long", "lat")) %>%
  mutate(has_temperature = !is.na(Z), block = assign_block(year)) %>%
  group_by(block, C) %>%
  summarise(
    precipitation_site_years = n(),
    with_temperature = sum(has_temperature),
    without_temperature = sum(!has_temperature),
    .groups = "drop"
  ) %>%
  arrange(match(block, blocks), C)

# Preserve original identifiers and descriptions for code-category auditing.
code_audit <- precip %>%
  select(year, long, lat, text_id, event_id, event_codes, description,
         rain, snow, hail) %>%
  pivot_longer(c(rain, snow, hail), names_to = "C", values_to = "present") %>%
  filter(present) %>%
  group_by(C) %>%
  slice_head(n = 30L) %>%
  ungroup() %>%
  select(-present)

readr::write_csv(conditional,
  file.path(output_dir, "precip_temperature_conditional.csv"))
readr::write_csv(conditional_by_block,
  file.path(output_dir, "precip_temperature_conditional_by_block.csv"))
readr::write_csv(sample_sizes,
  file.path(output_dir, "precip_temperature_site_year_counts.csv"))
readr::write_csv(candidates,
  file.path(output_dir, "precip_no_temperature_candidates.csv"))
readr::write_csv(code_audit,
  file.path(output_dir, "precip_event_code_audit_sample.csv"))

plot <- ggplot(conditional, aes(x = C, y = probability, fill = factor(Z))) +
  geom_col(width = 0.72) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Documentary precipitation category", y = "Proportion",
       fill = "Temperature index", title = "Temperature-index distribution by precipitation category",
       subtitle = "Same site and year; preliminary event-code classification") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "precip_temperature_conditional.png"),
       plot = plot, width = 8, height = 5, dpi = 160, bg = "white")

print(conditional)
message("Outputs saved to: ", output_dir)
