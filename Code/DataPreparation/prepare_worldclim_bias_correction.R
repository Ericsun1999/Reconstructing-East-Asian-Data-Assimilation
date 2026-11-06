here::i_am("Code/DataPreparation/prepare_worldclim_bias_correction.R")

# L1: WorldClim climatology correction for all 266 CESM-LME cells + 3 stations
# Keeps the historical repository unchanged; outputs only to Output/L1.
# Run from repository root:
# Rscript --vanilla Code/DataPreparation/L1_extract_modern_LME_WorldClim_full_grid_maskfix.R
suppressPackageStartupMessages({library(ncdf4); library(mgcv); library(terra)})

nc_dir <- path.expand(Sys.getenv('L1_LME_NC_DIR', '~/Downloads/LME nc data'))
wc_dir <- path.expand(Sys.getenv('L1_WORLDCLIM_DIR', '~/Downloads/wc2.1_30s_tavg'))
out_dir <- here::here('Output', 'L1')
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
years <- 1970:2000
lon_idx <- 40:58; lat_idx <- 58:71
nlon <- length(lon_idx); nlat <- length(lat_idx)
monthly_start <- (1970L - 1850L)*12L + 1L
monthly_count <- length(years)*12L
file_for <- function(j) file.path(nc_dir, sprintf(
  'b.e11.BLMTRC5CN.f19_g16.%03d.cam.h0.TREFHT.185001-200512.nc',j))

# Read native coordinates once; order MUST match lon-fastest NetCDF matrix layout.
f0 <- file_for(1L)
if (!file.exists(f0)) stop('Cannot find ',f0)
nc0 <- ncdf4::nc_open(f0)
lon <- as.numeric(ncdf4::ncvar_get(nc0,'lon'))[lon_idx]
lat <- as.numeric(ncdf4::ncvar_get(nc0,'lat'))[lat_idx]
ncdf4::nc_close(nc0)
coords <- data.frame(longitude=rep(lon,times=nlat),
                     latitude=rep(lat,each=nlon))
# Stable location IDs for 266 nodes. Use exact same coordinate grid in historical outputs.
grid <- data.frame(site=sprintf('grid_%03d',seq_len(nrow(coords))),
                   longitude=coords$longitude, latitude=coords$latitude,
                   type='grid')
stations <- data.frame(site=c('HongKong','Shanghai','Beijing'),
    longitude=c(114.167,121.433,116.283),
    latitude=c(22.333,31.167,39.933),type='station')
locations <- rbind(grid, stations)
stopifnot(nrow(grid)==266L, !anyDuplicated(locations$site))
write.csv(locations,file.path(out_dir,'L1_locations_all_269.csv'),row.names=FALSE)

# Extract monthly WorldClim values at each position. Do not construct a
# full-resolution global annual raster (30s rasters are very large).
wc_files <- file.path(wc_dir,sprintf('wc2.1_30s_tavg_%02d.tif',1:12))
if (!all(file.exists(wc_files))) stop('Missing WorldClim files: ',
  paste(wc_files[!file.exists(wc_files)],collapse=', '))
pts <- terra::vect(locations, geom=c('longitude','latitude'),crs='EPSG:4326')
wc_monthly <- vapply(wc_files, function(f) {
  r <- terra::rast(f)
  as.numeric(terra::extract(r,pts,method='simple')[,2])
},numeric(nrow(locations)))
# WorldClim land-temperature GeoTIFFs are NA over the ocean. Preserve
# those NAs; never interpret missing WorldClim as 0 degrees Celsius.
wc_valid <- rowSums(is.finite(wc_monthly)) == 12L
message('WorldClim complete locations: ',sum(wc_valid), '/',length(wc_valid))
message('WorldClim masked locations: ',sum(!wc_valid))
write.csv(locations[!wc_valid,],
  file.path(out_dir,'L1_worldclim_masked_locations.csv'),row.names=FALSE)
if (any(!wc_valid & locations$type == 'station')) {
  stop('WorldClim missing at a validation station; verify station coordinates/rasters')
}
if (any(wc_monthly[wc_valid,] < -80 | wc_monthly[wc_valid,] > 60))
  stop('Inspect WorldClim units/scaling')
colnames(wc_monthly) <- sprintf('WorldClim_month_%02d_C',1:12)
locations$WorldClim_1970_2000_C <- rowMeans(wc_monthly)
locations$worldclim_available <- wc_valid
write.csv(cbind(locations,wc_monthly),
  file.path(out_dir,'L1_worldclim_monthly_all_269.csv'),row.names=FALSE)

# Thin-plate regression-spline linear operator: same form as station pilot.
# The target includes every native grid cell and all validation stations.
train_coords <- data.frame(long=coords$longitude,lati=coords$latitude)
smooth <- mgcv::smoothCon(mgcv::s(long,lati,bs='tp',k=266),
  data=train_coords,absorb.cons=TRUE,scale.penalty=TRUE)[[1]]
train <- cbind(intercept=1,smooth$X)
target <- cbind(intercept=1,mgcv::PredictMat(smooth,
  data.frame(long=locations$longitude,lati=locations$latitude)))
penalty <- matrix(0,ncol(train),ncol(train))
penalty[-1,-1] <- smooth$S[[1]]
normal <- crossprod(train)+1e-8*penalty
if (rcond(normal)<1e-14) stop('Spline spatial operator poorly conditioned')
operator <- target %*% solve(normal,t(train))

# Modern CESM-LME annual temperatures for 13 members and 31 years.
# Annual weighting is arithmetic mean of 12 monthly means (same as helper).
rows <- vector('list',13L)
for(j in 1:13) {
  file <- file_for(j)
  if (!file.exists(file)) stop('Missing NetCDF: ',file)
  message('Reading TREFHT for a',j)
  nc <- ncdf4::nc_open(file)
  lon_j <- as.numeric(ncdf4::ncvar_get(nc,'lon'))[lon_idx]
  lat_j <- as.numeric(ncdf4::ncvar_get(nc,'lat'))[lat_idx]
  if (!isTRUE(all.equal(lon,lon_j,tolerance=1e-10)) ||
      !isTRUE(all.equal(lat,lat_j,tolerance=1e-10))) {
    ncdf4::nc_close(nc);stop('Coordinates differ: ',file)
  }
  u <- ncdf4::ncatt_get(nc,'TREFHT','units')
  if (!isTRUE(u$hasatt)) {ncdf4::nc_close(nc);stop('Missing units metadata')}
  arr <- ncdf4::ncvar_get(nc,'TREFHT',
    start=c(min(lon_idx),min(lat_idx),monthly_start),
    count=c(nlon,nlat,monthly_count),collapse_degen=FALSE)
  ncdf4::nc_close(nc)
  if (!identical(as.integer(dim(arr)),c(nlon,nlat,monthly_count)))
    stop('Unexpected TREFHT dimensions: ',file)
  unit <- tolower(trimws(as.character(u$value)))
  if (unit %in% c('k','kelvin')) arr <- arr-273.15
  else if (!unit %in% c('c','degc','degree_celsius','degrees_celsius'))
    stop('Unexpected temperature units: ',unit)
  if (any(!is.finite(arr))) stop('Nonfinite TREFHT: ',file)
  native <- matrix(arr,nrow=nrow(coords),ncol=monthly_count)
  annual <- vapply(seq_along(years),function(y)
    rowMeans(native[,((y-1)*12+1):(y*12),drop=FALSE]),
    numeric(nrow(coords)))
  target_annual <- operator %*% annual
  rows[[j]] <- data.frame(site=rep(locations$site,times=length(years)),
    member=paste0('a',j),year=rep(years,each=nrow(locations)),
    temp_C=as.vector(target_annual))
  rm(arr,native,annual,target_annual);gc(verbose=FALSE)
}
lme <- do.call(rbind,rows)
stopifnot(nrow(lme)==13L*31L*nrow(locations),all(is.finite(lme$temp_C)))
write.csv(lme,file.path(out_dir,'L1_lme_annual_1970_2000_all_269.csv'),row.names=FALSE)
clim <- aggregate(temp_C~site,lme,mean)
names(clim)[2] <- 'LME_1970_2000_C'
offsets <- merge(locations,clim,by='site',sort=FALSE)
offsets$d_C <- offsets$WorldClim_1970_2000_C-offsets$LME_1970_2000_C

# Extrapolate *offsets*, not WorldClim temperatures, to masked ocean cells.
# This is a sensitivity-oriented spatial fallback, NOT an observed WorldClim
# climatology over ocean. Use distances in km and label approximated values.
# For a land-only scientific product, you may instead exclude those cells.
offsets$correction_method <- ifelse(offsets$worldclim_available,
  'WorldClim_direct','IDW_from_land_offsets')
offsets$distance_to_land_km <- ifelse(offsets$worldclim_available,0,NA_real_)
source_idx <- which(offsets$worldclim_available & offsets$type=='grid')
if (length(source_idx) < 4L) stop('Fewer than 4 land-grid offsets available')
target_idx <- which(!offsets$worldclim_available)
for (i in target_idx) {
  # Great-circle distances from masked site to each valid land LME grid node.
  lon1 <- offsets$longitude[i]*pi/180
  lat1 <- offsets$latitude[i]*pi/180
  lon2 <- offsets$longitude[source_idx]*pi/180
  lat2 <- offsets$latitude[source_idx]*pi/180
  a <- sin((lat2-lat1)/2)^2 + cos(lat1)*cos(lat2)*sin((lon2-lon1)/2)^2
  km <- 6371*2*asin(pmin(1,sqrt(pmax(0,a))))
  nearest <- order(km)[seq_len(min(4L,length(km)))]
  distance <- km[nearest]
  w <- 1/pmax(distance,1)^2
  offsets$d_C[i] <- sum(w*offsets$d_C[source_idx[nearest]])/sum(w)
  offsets$distance_to_land_km[i] <- min(km)
}
stopifnot(all(is.finite(offsets$d_C)))
message('Masked grid offsets estimated by IDW: ',length(target_idx),
        '; farthest nearest-land grid distance (km): ',
        if (length(target_idx)) round(max(offsets$distance_to_land_km[target_idx]),1) else 0)
write.csv(offsets,file.path(out_dir,'L1_worldclim_offsets_all_269.csv'),row.names=FALSE)
write.csv(offsets[offsets$type=='grid',],
  file.path(out_dir,'L1_worldclim_offsets_grid266.csv'),row.names=FALSE)
write.csv(offsets[offsets$type=='station',],
  file.path(out_dir,'L1_worldclim_offsets_city3.csv'),row.names=FALSE)
print(offsets[offsets$type=='station',c('site','WorldClim_1970_2000_C',
  'LME_1970_2000_C','d_C')],row.names=FALSE)
message('Done: 266 grid cells + 3 validation stations; files saved in ',out_dir)
