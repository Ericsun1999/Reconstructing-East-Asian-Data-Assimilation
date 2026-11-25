here::i_am("Code/DataPreparation/prepare_precipitation_augmented_reaches.R")

# Construct a site-year temperature-index dataset supported by
# annual documentary precipitation records.
#
# Existing temperature indices are retained. A site-year with an
# annual precipitation record but no temperature record receives
# a working index of zero. The source flag distinguishes this
# assumption from an observed temperature category.
#
# The original spreadsheets are read-only. All generated products
# are saved under Output/Intermediate/REACHES/.

library(here)
library(readxl)
library(readr)
library(dplyr)

first_nonempty_unique <- function(x) {
  x <- trimws(as.character(x))
  paste(unique(x[!is.na(x) & nzchar(x)]), collapse = ";")
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

precip_file <- path.expand(Sys.getenv(
  "REACHES_PRECIPITATION_FILE",
  unset = here::here("Data", "precipitation index value.xlsx")
))
temperature_file <- here::here("Data", "temperature index value.v1.xlsx")
output_dir <- here::here("Output", "Intermediate", "REACHES")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (file in c(precip_file, temperature_file)) {
  if (!file.exists(file)) stop("Required input file not found: ", file)
}

# Read the annual worksheets only; do not combine annual and
# seasonal summaries.
# Coordinates are rounded to eight decimal places in both inputs
# before site-year matching. This removes negligible floating-point
# discrepancies that would otherwise create duplicate site-years.
temp_sheet <- readxl::read_excel(
  temperature_file, sheet = "全年溫度分級", col_types = "text"
)
precip_sheet <- readxl::read_excel(
  precip_file, sheet = "全年降水分級", col_types = "text"
)

required_temp <- c(
  "溫度等級", "資料判斷後，所屬年份", "空間-經度",
  "空間-緯度", "資料編號-文本ID", "資料編號-事件編碼"
)
required_precip <- c(
  "資料判斷後，所屬年份", "空間-經度", "空間-緯度",
  "資料編號-文本ID", "資料編號-事件編碼", "事件編碼"
)
if (!all(required_temp %in% names(temp_sheet))) {
  stop("Missing temperature columns: ",
       paste(setdiff(required_temp, names(temp_sheet)), collapse = ", "))
}
if (!all(required_precip %in% names(precip_sheet))) {
  stop("Missing precipitation columns: ",
       paste(setdiff(required_precip, names(precip_sheet)), collapse = ", "))
}

temperature <- temp_sheet %>%
  transmute(
    source_order = row_number(),
    year = suppressWarnings(as.integer(.data[["資料判斷後，所屬年份"]])),
    long = round(suppressWarnings(as.numeric(.data[["空間-經度"]])), 8),
    lat = round(suppressWarnings(as.numeric(.data[["空間-緯度"]])), 8),
    level = suppressWarnings(as.integer(.data[["溫度等級"]])),
    temperature_text_id = .data[["資料編號-文本ID"]],
    temperature_event_id = .data[["資料編號-事件編碼"]]
  ) %>%
  filter(
    between(year, 1368L, 1911L),
    level %in% c(-2L, -1L, 0L, 1L),
    is.finite(long), is.finite(lat)
  )

# Site-year keys only determine eligibility for one inferred zero.
# All original temperature records, including repeated site-years and
# recorded zeros, remain separate observations.
temperature_site_year <- temperature %>%
  distinct(year, long, lat)

precipitation <- precip_sheet %>%
  transmute(
    year = suppressWarnings(as.integer(.data[["資料判斷後，所屬年份"]])),
    long = round(suppressWarnings(as.numeric(.data[["空間-經度"]])), 8),
    lat = round(suppressWarnings(as.numeric(.data[["空間-緯度"]])), 8),
    precipitation_text_id = .data[["資料編號-文本ID"]],
    precipitation_event_id = .data[["資料編號-事件編碼"]],
    precipitation_event_codes = .data[["事件編碼"]]
  ) %>%
  filter(between(year, 1368L, 1911L), is.finite(long), is.finite(lat))

# Every recorded annual precipitation event is treated as evidence
# of documentary availability, regardless of its index value or
# tentative rain/snow/hail classification. Original event codes
# are retained for later event-type sensitivity checks.
precipitation_site_year <- precipitation %>%
  group_by(year, long, lat) %>%
  summarise(
    number_of_precipitation_records = n(),
    precipitation_text_ids = first_nonempty_unique(precipitation_text_id),
    precipitation_event_ids = first_nonempty_unique(precipitation_event_id),
    precipitation_event_codes = first_nonempty_unique(precipitation_event_codes),
    .groups = "drop"
  )

# Only site-years with a documentary precipitation record and no
# temperature observation may receive an inferred zero.
inferred_site_year <- precipitation_site_year %>%
  anti_join(
    temperature_site_year,
    by = c("year", "long", "lat")
  )

augmented <- bind_rows(
  temperature %>%
    transmute(
      source_order, year, long, lat, level,
      source_type = "observed_temperature",
      is_inferred_zero = FALSE,
      number_of_source_records = 1L,
      number_of_temperature_records = 1L,
      number_of_precipitation_records = NA_integer_,
      supporting_event_types = "temperature",
      temperature_text_ids = as.character(temperature_text_id),
      temperature_event_ids = as.character(temperature_event_id),
      precipitation_text_ids = NA_character_,
      precipitation_event_ids = NA_character_,
      precipitation_event_codes = NA_character_
    ),
  inferred_site_year %>%
    transmute(
      source_order = NA_integer_,
      year, long, lat, level = 0L,
      source_type = "inferred_zero_precipitation",
      is_inferred_zero = TRUE,
      number_of_source_records = 0L,
      number_of_temperature_records = 0L,
      number_of_precipitation_records,
      supporting_event_types = "annual_precipitation",
      temperature_text_ids = NA_character_,
      temperature_event_ids = NA_character_,
      precipitation_text_ids, precipitation_event_ids,
      precipitation_event_codes
    )
) %>%
  arrange(year, long, lat, is_inferred_zero, source_order)

# A recorded temperature of zero blocks imputation just like any other
# observed category. Duplicate originals are allowed; inferred zeros are not.
key_columns <- c("year", "long", "lat")
if (anyDuplicated(temperature_site_year[key_columns]) ||
    anyDuplicated(precipitation_site_year[key_columns]) ||
    anyDuplicated(inferred_site_year[key_columns])) {
  stop("The source site-year key indexes must be unique.")
}
inferred_check <- augmented %>% filter(is_inferred_zero)
if (anyDuplicated(inferred_check[key_columns]) ||
    nrow(semi_join(inferred_check, temperature_site_year,
                   by = key_columns)) != 0L) {
  stop("Inferred zeros must be unique and have no recorded temperature.")
}
observed_check <- augmented %>%
  filter(!is_inferred_zero) %>%
  select(source_order, year, long, lat, level) %>%
  arrange(source_order)
original_check <- temperature %>%
  select(source_order, year, long, lat, level) %>%
  arrange(source_order)
if (!identical(as.data.frame(observed_check),
               as.data.frame(original_check))) {
  stop("Original temperature observations were dropped or altered.")
}
if (nrow(augmented) != nrow(temperature) + nrow(inferred_site_year) ||
    any(augmented$is_inferred_zero & augmented$level != 0L) ||
    any(augmented$is_inferred_zero &
        augmented$number_of_source_records != 0L) ||
    any(!augmented$is_inferred_zero &
        augmented$number_of_source_records != 1L)) {
  stop("Augmented record counts or inferred-zero values are invalid.")
}

output_file <- file.path(output_dir, "temperature_site_year_augmented.csv")
readr::write_csv(augmented, output_file)

inferred_audit <- augmented %>% filter(is_inferred_zero)
readr::write_csv(
  inferred_audit,
  file.path(output_dir, "temperature_inferred_zero_audit.csv")
)

summary_by_block <- augmented %>%
  mutate(
    block = assign_block(year),
    is_observed = !is_inferred_zero
  ) %>%
  group_by(block) %>%
  summarise(
    observed_temperature_records = sum(is_observed),
    observed_temperature_site_years =
      n_distinct(year[is_observed], long[is_observed], lat[is_observed]),
    inferred_zero_site_years = sum(is_inferred_zero),
    augmented_temperature_records = n(),
    augmented_site_years = n_distinct(year, long, lat),
    informative_years = n_distinct(year),
    .groups = "drop"
  ) %>%
  arrange(match(
    block, c("1368-1500", "1501-1600", "1601-1700",
             "1701-1800", "1801-1911")
  ))
readr::write_csv(
  summary_by_block,
  file.path(output_dir, "temperature_augmented_by_block.csv")
)

yearly_summary <- augmented %>%
  group_by(year) %>%
  summarise(
    recorded_temperature = sum(!is_inferred_zero),
    inferred_zero = sum(is_inferred_zero),
    total = n(),
    .groups = "drop"
  ) %>%
  right_join(tibble::tibble(year = 1368:1911), by = "year") %>%
  mutate(
    recorded_temperature = coalesce(recorded_temperature, 0L),
    inferred_zero = coalesce(inferred_zero, 0L),
    total = coalesce(total, 0L)
  ) %>%
  arrange(year)
readr::write_csv(
  yearly_summary,
  file.path(output_dir, "temperature_augmented_by_year.csv")
)

metadata <- tibble::tibble(
  source = "annual_precipitation_documentary_availability",
  temperature_file = temperature_file,
  precipitation_file = precip_file,
  observed_temperature_records = nrow(temperature),
  observed_temperature_site_years = nrow(temperature_site_year),
  precipitation_site_years = nrow(precipitation_site_year),
  inferred_zero_site_years = nrow(inferred_site_year),
  augmented_temperature_records = nrow(augmented),
  augmented_site_years = nrow(temperature_site_year) +
    nrow(inferred_site_year),
  observed_informative_years = n_distinct(temperature_site_year$year),
  augmented_informative_years = n_distinct(augmented$year)
)
readr::write_csv(
  metadata,
  file.path(output_dir, "temperature_augmented_metadata.csv")
)

print(metadata)
print(summary_by_block)
message("Augmented REACHES data saved to: ", output_file)
message("The original temperature spreadsheet was not modified.")
