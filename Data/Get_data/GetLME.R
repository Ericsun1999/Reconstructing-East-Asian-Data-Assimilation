library(ncdf4,tidync) # for NetCDF files
library(lubridate)
library(zoo)

data1 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.001.cam.h0.TREFHT.185001-200512.nc")
data2 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.002.cam.h0.TREFHT.085001-184912.nc")
data3 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.003.cam.h0.TREFHT.085001-184912.nc")
data4 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.004.cam.h0.TREFHT.085001-184912.nc")
data5 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.005.cam.h0.TREFHT.085001-184912.nc")
data6 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.006.cam.h0.TREFHT.185001-200512.nc")
data7 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.007.cam.h0.TREFHT.085001-184912.nc")
data8 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.008.cam.h0.TREFHT.085001-184912.nc")
data9 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.009.cam.h0.TREFHT.085001-184912.nc")
data10 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.010.cam.h0.TREFHT.085001-184912.nc")
data11 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.011.cam.h0.TREFHT.185001-200512.nc")
data12 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.012.cam.h0.TREFHT.085001-184912.nc")
data13 <- nc_open("~/Downloads/b.e11.BLMTRC5CN.f19_g16.013.cam.h0.TREFHT.085001-184912.nc")


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

data1 <- Data1[c(46,144,227),]
data2 <- Data2[c(46,144,227),]
data3 <- Data3[c(46,144,227),]
data4 <- Data4[c(46,144,227),]
data5 <- Data5[c(46,144,227),]
data6 <- Data6[c(46,144,227),]
data7 <- Data7[c(46,144,227),]
data8 <- Data8[c(46,144,227),]
data9 <- Data9[c(46,144,227),]
data10 <- Data10[c(46,144,227),]
data11 <- Data11[c(46,144,227),]
data12 <- Data12[c(46,144,227),]
data13 <- Data13[c(46,144,227),]

i=1
Dataa1<-rbind(data1[i,],data2[i,],data3[i,],data4[i,],data5[i,],data6[i,],data7[i,],data8[i,],data9[i,],data10[i,],data11[i,],data12[i,],data13[i,])

i=2
Dataa2<-rbind(data1[i,],data2[i,],data3[i,],data4[i,],data5[i,],data6[i,],data7[i,],data8[i,],data9[i,],data10[i,],data11[i,],data12[i,],data13[i,])

i=3
Dataa3<-rbind(data1[i,],data2[i,],data3[i,],data4[i,],data5[i,],data6[i,],data7[i,],data8[i,],data9[i,],data10[i,],data11[i,],data12[i,],data13[i,])

write.csv(Dataa1, "a1.csv")
write.csv(Dataa2, "a2.csv")
write.csv(Dataa3, "a3.csv")


#Get 'LME data/Figure9'
clean_Data <- function(data){
  data <- data[data$long < 125.5 & data$long >99 & data$lati > 19,]
  data <- data[-c(1,2,3,4,6,7,8,9,10,11,12,13,14,15,18,19,20,21,22,32,33,44,45,55,56,66,67,77,78,79,88,89,98,99,100,110,111,112,113,114,121,122,123,124,125,126,132,133,134,135,136,137,138,139,140,141,143 ),]
}

data1 <- clean_Data(Data1)
data2 <- clean_Data(Data2)
data3 <- clean_Data(Data3)
data4 <- clean_Data(Data4)
data5 <- clean_Data(Data5)
data6 <- clean_Data(Data6)
data7 <- clean_Data(Data7)
data8 <- clean_Data(Data8)
data9 <- clean_Data(Data9)
data10 <- clean_Data(Data10)
data11 <- clean_Data(Data11)
data12 <- clean_Data(Data12)
data13 <- clean_Data(Data13)


write.csv(data1, "b1.csv")
write.csv(data2, "b2.csv")
write.csv(data3, "b3.csv")
write.csv(data4, "b4.csv")
write.csv(data5, "b5.csv")
write.csv(data6, "b6.csv")
write.csv(data7, "b7.csv")
write.csv(data8, "b8.csv")
write.csv(data9, "b9.csv")
write.csv(data10, "b10.csv")
write.csv(data11, "b11.csv")
write.csv(data12, "b12.csv")
write.csv(data13, "b13.csv")

