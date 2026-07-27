library(ncdf4)

# ============================================================
# 1. path
# ============================================================

data_dir <- path.expand("./LME nc data")

output_dir <- file.path(
  data_dir,
  "processed_1350_1949"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 19 longitude × 14 latitude = 266 locations
lon_idx <- 40:58
lat_idx <- 58:71

stopifnot(length(lon_idx) == 19L)
stopifnot(length(lat_idx) == 14L)

# 1350-01 - 1949-12， 600 years、7200 months
month_names <- sprintf(
  "T_%04d_%02d",
  rep(1350:1949, each = 12),
  rep(1:12, times = 600)
)

stopifnot(length(month_names) == 7200L)


# ============================================================
# 2. 
# ============================================================

get_member_files <- function(member_number) {

  member <- sprintf("%03d", member_number)

  pre1850_file <- file.path(
    data_dir,
    sprintf(
      paste0(
        "b.e11.BLMTRC5CN.f19_g16.%s.",
        "cam.h0.TREFHT.085001-184912.nc"
      ),
      member
    )
  )

  post1850_file <- file.path(
    data_dir,
    sprintf(
      paste0(
        "b.e11.BLMTRC5CN.f19_g16.%s.",
        "cam.h0.TREFHT.185001-200512.nc"
      ),
      member
    )
  )

  list(
    pre1850 = pre1850_file,
    post1850 = post1850_file
  )
}


# ============================================================
# 3. 
# ============================================================

read_trefht_subset <- function(
    nc_file,
    lon_idx,
    lat_idx,
    time_start,
    time_count,
    convert_to_celsius = TRUE) {

  if (!file.exists(nc_file)) {
    stop(
      "no NetCDF file：\n",
      nc_file
    )
  }

  nc <- nc_open(nc_file)
  on.exit(nc_close(nc), add = TRUE)

  if (!"TREFHT" %in% names(nc$var)) {
    stop(
      "no TREFHT：\n",
      nc_file
    )
  }

  dim_names <- vapply(
    nc$var$TREFHT$dim,
    function(x) x$name,
    character(1)
  )


  lon_all <- as.numeric(
    ncvar_get(nc, "lon")
  )

  lat_all <- as.numeric(
    ncvar_get(nc, "lat")
  )

  n_time_available <- nc$var$TREFHT$dim[[3]]$len
  time_end <- time_start + time_count - 1L


  temperature <- ncvar_get(
    nc,
    "TREFHT",
    start = c(
      min(lon_idx),
      min(lat_idx),
      time_start
    ),
    count = c(
      length(lon_idx),
      length(lat_idx),
      time_count
    ),
    collapse_degen = FALSE
  )

  units_result <- ncatt_get(
    nc,
    "TREFHT",
    "units"
  )

  if (isTRUE(units_result$hasatt)) {
    units <- as.character(units_result$value)
  } else {
    units <- NA_character_
  }

  if (convert_to_celsius) {

    units_lower <- tolower(trimws(units))

    if (units_lower %in% c("k", "kelvin")) {

      temperature <- temperature - 273.15

    } else if (
      units_lower %in%
        c(
          "c",
          "degc",
          "degree_celsius",
          "degrees_celsius"
        )
    ) {


    } else {

      stop(
        paste0(
          "Can't identify TREFHT unit：",
          units
        )
      )
    }
  }

  list(
    temperature = temperature,
    lon = as.numeric(lon_all[lon_idx]),
    lat = as.numeric(lat_all[lat_idx]),
    original_units = units
  )
}


# ============================================================
# 4.
# ============================================================

process_lme_member <- function(
    member_number,
    save_csv = TRUE,
    save_rds = TRUE,
    convert_to_celsius = TRUE) {


  member_number <- as.integer(member_number)
  member <- sprintf("%03d", member_number)

  files <- get_member_files(member_number)

  cat("\n")
  cat("========================================\n")
  cat("Processing member", member, "\n")
  cat("========================================\n")

  cat(
    "Pre-1850 file:\n",
    files$pre1850,
    "\n\n"
  )

  cat(
    "Post-1850 file:\n",
    files$post1850,
    "\n\n"
  )


  pre1850 <- read_trefht_subset(
    nc_file = files$pre1850,
    lon_idx = lon_idx,
    lat_idx = lat_idx,
    time_start = 6001L,
    time_count = 6000L,
    convert_to_celsius = convert_to_celsius
  )


  post1850 <- read_trefht_subset(
    nc_file = files$post1850,
    lon_idx = lon_idx,
    lat_idx = lat_idx,
    time_start = 1L,
    time_count = 1200L,
    convert_to_celsius = convert_to_celsius
  )


  lon_selected <- as.numeric(pre1850$lon)
  lat_selected <- as.numeric(pre1850$lat)

  n_lon <- length(lon_selected)
  n_lat <- length(lat_selected)

  stopifnot(n_lon == 19L)
  stopifnot(n_lat == 14L)



  full_temperature <- array(
    NA_real_,
    dim = c(
      n_lon,
      n_lat,
      7200L
    )
  )

  full_temperature[, , 1:6000] <-
    pre1850$temperature

  full_temperature[, , 6001:7200] <-
    post1850$temperature


  coordinate_lon <- rep(
    lon_selected,
    times = n_lat
  )

  coordinate_lat <- rep(
    lat_selected,
    each = n_lon
  )

  stopifnot(length(coordinate_lon) == 266L)
  stopifnot(length(coordinate_lat) == 266L)

  # location × time
  temperature_matrix <- matrix(
    full_temperature,
    nrow = n_lon * n_lat,
    ncol = 7200L
  )

  colnames(temperature_matrix) <- month_names

  output <- data.frame(
    lon = coordinate_lon,
    lat = coordinate_lat,
    temperature_matrix,
    check.names = FALSE
  )


  stopifnot(nrow(output) == 266L)
  stopifnot(ncol(output) == 7202L)

  stopifnot(names(output)[1] == "lon")
  stopifnot(names(output)[2] == "lat")
  stopifnot(names(output)[3] == "T_1350_01")
  stopifnot(names(output)[7202] == "T_1949_12")


  stopifnot(
    max(
      abs(
        as.numeric(output$lon) -
          coordinate_lon
      ),
      na.rm = TRUE
    ) < 1e-10
  )

  stopifnot(
    max(
      abs(
        as.numeric(output$lat) -
          coordinate_lat
      ),
      na.rm = TRUE
    ) < 1e-10
  )

  stopifnot(
    isTRUE(
      all.equal(
        as.numeric(output[1, 3]),
        as.numeric(full_temperature[1, 1, 1]),
        tolerance = 1e-10,
        check.attributes = FALSE
      )
    )
  )

  stopifnot(
    isTRUE(
      all.equal(
        as.numeric(output[n_lon, 3]),
        as.numeric(
          full_temperature[n_lon, 1, 1]
        ),
        tolerance = 1e-10,
        check.attributes = FALSE
      )
    )
  )

  stopifnot(
    isTRUE(
      all.equal(
        as.numeric(output[n_lon + 1L, 3]),
        as.numeric(
          full_temperature[1, 2, 1]
        ),
        tolerance = 1e-10,
        check.attributes = FALSE
      )
    )
  )

  stopifnot(
    isTRUE(
      all.equal(
        as.numeric(output[266, 3]),
        as.numeric(
          full_temperature[n_lon, n_lat, 1]
        ),
        tolerance = 1e-10,
        check.attributes = FALSE
      )
    )
  )

  stopifnot(names(output)[6002] == "T_1849_12")
  stopifnot(names(output)[6003] == "T_1850_01")

  stopifnot(
    isTRUE(
      all.equal(
        as.numeric(output[1, 6002]),
        as.numeric(
          pre1850$temperature[1, 1, 6000]
        ),
        tolerance = 1e-10,
        check.attributes = FALSE
      )
    )
  )

  stopifnot(
    isTRUE(
      all.equal(
        as.numeric(output[1, 6003]),
        as.numeric(
          post1850$temperature[1, 1, 1]
        ),
        tolerance = 1e-10,
        check.attributes = FALSE
      )
    )
  )


  csv_file <- file.path(
    output_dir,
    sprintf("a%d.csv", member_number)
  )

  rds_file <- file.path(
    output_dir,
    sprintf("a%d.rds", member_number)
  )

  if (save_csv) {

    write.csv(
      output,
      file = csv_file,
      row.names = FALSE,
      quote = FALSE,
      na = ""
    )

    cat(
      "CSV saved:\n",
      csv_file,
      "\n"
    )
  }

  if (save_rds) {

    saveRDS(
      output,
      file = rds_file,
      compress = FALSE
    )

    cat(
      "RDS saved:\n",
      rds_file,
      "\n"
    )
  }

  cat(
    "Output dimension:",
    nrow(output),
    "rows ×",
    ncol(output),
    "columns\n"
  )

  cat(
    "Temperature period:",
    names(output)[3],
    "to",
    names(output)[ncol(output)],
    "\n"
  )

  cat(
    "Longitude range:",
    min(output$lon),
    "to",
    max(output$lon),
    "\n"
  )

  cat(
    "Latitude range:",
    min(output$lat),
    "to",
    max(output$lat),
    "\n"
  )

  cat(
    "Temperature range:",
    round(
      min(
        as.matrix(output[, 3:7202]),
        na.rm = TRUE
      ),
      3
    ),
    "to",
    round(
      max(
        as.matrix(output[, 3:7202]),
        na.rm = TRUE
      ),
      3
    ),
    if (convert_to_celsius) "°C" else "",
    "\n"
  )

  invisible(
    list(
      data = output,
      csv_file = csv_file,
      rds_file = rds_file
    )
  )
}

for (member_number in 1:13) {

  process_lme_member(
    member_number = member_number,
    save_csv = TRUE,
    save_rds = TRUE,
    convert_to_celsius = TRUE
  )

  gc()
}
