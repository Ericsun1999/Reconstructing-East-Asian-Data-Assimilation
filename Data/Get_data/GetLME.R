library(ncdf4,tidync) # for NetCDF files
library(lubridate)
library(zoo)

data1 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.001.cam.h0.TREFHT.185001-200512.nc")
data2 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.002.cam.h0.TREFHT.085001-184912.nc")
data3 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.003.cam.h0.TREFHT.085001-184912.nc")
data4 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.004.cam.h0.TREFHT.085001-184912.nc")
data5 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.005.cam.h0.TREFHT.085001-184912.nc")
data6 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.006.cam.h0.TREFHT.185001-200512.nc")
data7 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.007.cam.h0.TREFHT.085001-184912.nc")
data8 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.008.cam.h0.TREFHT.085001-184912.nc")
data9 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.009.cam.h0.TREFHT.085001-184912.nc")
data10 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.010.cam.h0.TREFHT.085001-184912.nc")
data11 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.011.cam.h0.TREFHT.185001-200512.nc")
data12 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.012.cam.h0.TREFHT.085001-184912.nc")
data13 <- nc_open("./b.e11.BLMTRC5CN.f19_g16.013.cam.h0.TREFHT.085001-184912.nc")


TREFHT1<-ncvar_get(data1,"TREFHT")
TREFHT2<-ncvar_get(data2,"TREFHT")
TREFHT3<-ncvar_get(data3,"TREFHT")
TREFHT4<-ncvar_get(data4,"TREFHT")
TREFHT5<-ncvar_get(data5,"TREFHT")
TREFHT6<-ncvar_get(data6,"TREFHT")
TREFHT7<-ncvar_get(data7,"TREFHT")
TREFHT8<-ncvar_get(data8,"TREFHT")
TREFHT9<-ncvar_get(data9,"TREFHT")
TREFHT10<-ncvar_get(data10,"TREFHT")
TREFHT11<-ncvar_get(data11,"TREFHT")
TREFHT12<-ncvar_get(data12,"TREFHT")
TREFHT13<-ncvar_get(data13,"TREFHT")


#Data preprocess

time <- ncvar_get(data1,"time")
lon <- ncvar_get(data1,"lon")
lat <- ncvar_get(data1,"lat")

  long<-lon[40:58]
  lat1<-lat[58:71]
  lati<-replicate(19,lat1)
  lati<-array(t(lati))
  
  clean_LMEdata <- function(data){
    a1<-data[40:58,58:71,6001:12000]
    mat.a1 <- matrix(a1,19*14,1200)
    a1<-data.frame(mat.a1)
    a1<-cbind(lati,long,a1)
  }
  
  Data1 <- clean_LMEdata(TREFHT1)
  Data2 <- clean_LMEdata(TREFHT2)
  Data3 <- clean_LMEdata(TREFHT3)
  Data4 <- clean_LMEdata(TREFHT4)
  Data5 <- clean_LMEdata(TREFHT5)
  Data6 <- clean_LMEdata(TREFHT6)
  Data7 <- clean_LMEdata(TREFHT7)
  Data8 <- clean_LMEdata(TREFHT8)
  Data9 <- clean_LMEdata(TREFHT9)
  Data10 <- clean_LMEdata(TREFHT10)
  Data11 <- clean_LMEdata(TREFHT11)
  Data12 <- clean_LMEdata(TREFHT12)
  Data13 <- clean_LMEdata(TREFHT13)

#Get LME data

clean_data <- function(data){
  data %>%
  rowwise() %>% 
  mutate(
    annual_data = list(
      colMeans(matrix(c_across(3:ncol(.)), nrow = 12))
    )
  ) %>%
  ungroup() %>%
  dplyr::select(lati, long, annual_data) %>%  
  tidyr::unnest_wider(annual_data, names_sep = "_")  %>%
    rename_with(
      ~ as.character(1350:1949),
      -c(lati, long)
    )
}

data1 <- clean_data(Data1)  
data2 <- clean_data(Data2)  
data3 <- clean_data(Data3) 
data4 <- clean_data(Data4) 
data5 <- clean_data(Data5)  
data6 <- clean_data(Data6)
data7 <- clean_data(Data7)  
data8 <- clean_data(Data8) 
data9 <- clean_data(Data9)
data10 <- clean_data(Data10)
data11 <- clean_data(Data11) 
data12 <- clean_data(Data12)
data13 <- clean_data(Data13) 

pred_points <- data.frame(
  long = c(113.75, 121.25, 116.25),
  lati = c(22.25, 31.25, 39.25)
)

library(mgcv)

spatial_spline_predict_one_year <- function(df, year_col, pred_points, k = 60) {
  d <- data.frame(
    y    = as.numeric(df[[year_col]]),
    long = as.numeric(df$long),
    lati = as.numeric(df$lati)
  )
  d <- d[is.finite(d$y) & is.finite(d$long) & is.finite(d$lati), ]

  fit <- gam(
    y ~ s(long, lati, bs = "tp", k = k),
    data = d,
    method = "REML"
  )

  predict(fit, newdata = pred_points)
}

predict_one_series <- function(dataj, pred_points, k = 60) {
  year_cols <- 3:ncol(dataj)

  preds <- sapply(year_cols, function(col_id) {
    spatial_spline_predict_one_year(dataj, col_id, pred_points, k = k)
  })

  rownames(preds) <- paste0("pt", 1:nrow(pred_points))
  colnames(preds) <- colnames(dataj)[year_cols]
  preds
}

data_list <- list(data1, data2, data3, data4, data5, data6, data7,
                  data8, data9, data10, data11, data12, data13)

result_list <- lapply(data_list, predict_one_series,
                      pred_points = pred_points, k = 60)


library(dplyr)

extract_point_df <- function(result_list, point_id,
                             long, lati) {
  # point_id: "pt1", "pt2", or "pt3"

  pt_df <- lapply(seq_along(result_list), function(j) {
    result_list[[j]][point_id, ]
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()

  rownames(pt_df) <- paste0("data", seq_along(result_list))

  pt_df <- pt_df %>%
    dplyr::mutate(
      long   = long,
      lati   = lati,
      series = rownames(pt_df),
      point  = point_id,
      .before = 1
    )

  pt_df
}

pt_coords <- data.frame(
  point = c("pt1", "pt2", "pt3"),
  long  = c(113.75, 121.25, 116.25),
  lati  = c(22.25, 31.25, 39.25)
)


pt1_df <- extract_point_df(
  result_list, "pt1",
  long = 113.75, lati = 22.25
)

pt2_df <- extract_point_df(
  result_list, "pt2",
  long = 121.25, lati = 31.25
)

pt3_df <- extract_point_df(
  result_list, "pt3",
  long = 116.25, lati = 39.25
)

pt1_df1 <- pt1_df[,-c(3:4)]
pt2_df1 <- pt2_df[,-c(3:4)]
pt3_df1 <- pt3_df[,-c(3:4)]

write.csv(pt1_df1, "d1.csv")
write.csv(pt2_df1, "d2.csv")
write.csv(pt3_df1, "d3.csv")



#Get 'LME data/Figure6e'

library(sf)
library(rnaturalearth)
library(dplyr)
library(stringr)


world <- ne_countries(scale = "medium", returnclass = "sf")
targets <- world %>%
  filter(str_detect(admin, regex("China|Taiwan|Hong Kong|Macao|Macau", ignore_case = TRUE)))
china_tw_hk_mo <- st_union(targets)  
china_tw_hk_mo <- st_as_sf(data.frame(name = "China+Taiwan+HK+MO"), geometry = china_tw_hk_mo)
st_crs(china_tw_hk_mo) <- 4326

coast_km <- 80

china_tw_hk_mo_valid <- st_make_valid(china_tw_hk_mo)
poly_m <- st_transform(china_tw_hk_mo_valid, 3857)
poly_m_buf <- st_buffer(poly_m, dist = coast_km * 1000)
china_tw_hk_mo_buf <- st_transform(poly_m_buf, 4326)

pts_sf <- st_as_sf(data1, coords = c("long", "lati"), crs = 4326)
data1$in_china_tw_hk_mo_coast <- st_within(
  pts_sf,
  china_tw_hk_mo_buf,
  sparse = FALSE
)[,1]

arr <- (data1$in_china_tw_hk_mo_coast == "TRUE")
arr[260] <- arr[261] <- FALSE
dataa1 <- data1[(arr), -c(3:20,566:603)]
dataa2 <- data2[(arr), -c(3:20,566:602)]
dataa3 <- data3[(arr), -c(3:20,566:602)]
dataa4 <- data4[(arr), -c(3:20,566:602)]
dataa5 <- data5[(arr), -c(3:20,566:602)]
dataa6 <- data6[(arr), -c(3:20,566:602)]
dataa7 <- data7[(arr), -c(3:20,566:602)]
dataa8 <- data8[(arr), -c(3:20,566:602)]
dataa9 <- data9[(arr), -c(3:20,566:602)]
dataa10 <- data10[(arr), -c(3:20,566:602)]
dataa11 <- data11[(arr), -c(3:20,566:602)]
dataa12 <- data12[(arr), -c(3:20,566:602)]
dataa13 <- data13[(arr), -c(3:20,566:602)]

dataa <- rbind(dataa1, dataa2, dataa3, dataa4, dataa5, dataa6, dataa7, dataa8, dataa9, dataa10, dataa11, dataa12, dataa13)
write.csv(dataa, "b0.csv")

