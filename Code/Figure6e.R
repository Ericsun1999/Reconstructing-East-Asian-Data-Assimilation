#Load LME data
Data1 <- read.csv("./b1.csv", row.names=1)
Data2 <- read.csv("./b2.csv", row.names=1)
Data3 <- read.csv("./b3.csv", row.names=1)
Data4 <- read.csv("./b4.csv", row.names=1)
Data5 <- read.csv("./b5.csv", row.names=1)
Data6 <- read.csv("./b6.csv", row.names=1)
Data7 <- read.csv("./b7.csv", row.names=1)
Data8 <- read.csv("./b8.csv", row.names=1)
Data9 <- read.csv("./b9.csv", row.names=1)
Data10 <- read.csv("./b10.csv", row.names=1)
Data11 <- read.csv("./b11.csv", row.names=1)
Data12 <- read.csv("./b12.csv", row.names=1)
Data13 <- read.csv("./b13.csv", row.names=1)

clean_data <- function(data){
  data %>%
  rowwise() %>% 
  mutate(
    annual_data = list(
      colMeans(matrix(c_across(3:ncol(.)), nrow = 12))
    )
  ) %>%
  ungroup() %>%
  select(lati, long, annual_data) %>%  
  tidyr::unnest_wider(annual_data, names_sep = "_")  
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

haa<-as.data.frame(matrix(0, nrow = nrow(data1), ncol = ncol(data1)-2))


for (i in 1:nrow(data1)) {
  for (j in 1:(ncol(data1)-2)) {
    haa[i,j] <- mean(as.numeric(data1[i,j+2]), as.numeric(data2[i,j+2]),
                     as.numeric(data3[i,j+2]), as.numeric(data4[i,j+2]),
                     as.numeric(data5[i,j+2]), as.numeric(data6[i,j+2]),
                     as.numeric(data7[i,j+2]), as.numeric(data8[i,j+2]),
                     as.numeric(data9[i,j+2]), as.numeric(data10[i,j+2]),
                     as.numeric(data11[i,j+2]), as.numeric(data12[i,j+2]),
                     as.numeric(data13[i,j+2]) ) -273
  }
}

#Load REACHES data
library(readxl)
temperature<- read_excel("~/Downloads/DA/temperature index value.v1.xlsx",col_type = c("skip","skip","numeric","numeric","skip","skip","skip","skip","skip","numeric","numeric","skip","skip"))
colnames(temperature)<-c("level","year","long","lat")

# Take mode for events at duplicated locations and rearrange the data by year
library(dplyr)  
temp2 <-temperature%>% 
    group_by(year,long,lat) %>%  
    summarise(level = mean(level))

# Compute empirical variograms based on all the data at various years
  library(sp)
  library(spacetime)
  library(zoo)
  library(gstat)

  y2<- temp2[temp2$year>100 & temp2$year<=2000,]
  coord <- unique(round(as.data.frame(cbind(long=y2$long,lat=y2$lat)),4))
  coordinates(coord) <- ~ long + lat
  proj4string(coord) <- CRS('+proj=longlat +datum=WGS84')
  year <- sort(unique(y2$year))
  YM <- as.yearmon(year)
  temp23 <- as.data.frame(round(as.matrix(y2),4))
  ind <- matrix(NA,nrow(temp23),2)
  ind[,2] <- match(temp23$year,year)
  ind[,1] <- match(temp23$long*10^7+temp23$lat,coord@coords[,1]*10^7+coord@coords[,2])
  stsdf2 <- STSDF(sp=coord,time=YM,data=y2,index=ind)
  vario.stsdf2 <- variogramST(level~1,data=stsdf2,tlags=0,width=10,na.omit=T)


# Semi-variogram model fitting

coordinates(y2) <- ~ long + lat
  proj4string(y2) <- CRS('+proj=longlat +datum=WGS84')
  vario2 <- variogram(level~1,data=y2,width=10)
  ind <- match(as.numeric(rownames(vario.stsdf2))-1,rownames(vario2))
  ind2 <- rep(F,nrow(vario2))
  ind2[ind] <- T
  vario2 <- vario2[ind2,]
  if(sum(vario.stsdf2$dist==0)>0) vario2[,1:3] <- vario.stsdf2[ind[-1],1:3] else vario2[,1:3] <- vario.stsdf2[ind,1:3]
  
  vario.fit2 <- fit.variogram(vario2,vgm(model="Exp",nugget=NA),fit.kappa=T,fit.method=7)

#Bias corrected
vario.fit2$psill <- c(0.146, 0.739)
vario.fit2$range <- c(0, 299.099)

# Kriging
year.start <- range(temp2$year)[1]
year.end <- range(temp2$year)[2]
year.all <- year.start:year.end
nyear <- year.end - year.start + 1
ncase1 <- rep(NA, nyear)
for (i in 1:nyear) {
  ncase1[i] <- sum(temp2$year == year.start + i - 1)
  if (ncase1[i] > 0) {
    temp21 <- temp2[temp2$year == (year.start + i - 1), ]
  }
}
m <- 1
year2 <- year.all[which(ncase1 >= m)]

arr.pred <- arr.std <- array(NA, c(12, 14, length(year2)))
loc <- expand.grid(long = seq(97.5, 125, by = 2.5), lat = seq(18, 42.63158, by = 1.89473692308))
coordinates(loc) <- ~ long + lat
proj4string(loc) <- CRS('+proj=longlat +datum=WGS84')

#past temperature data

for(i in 1:length(year2)) {
    temp <- y2[y2$year==year2[i],]
    y.pred <- krige(level~1,temp,loc,model=vario.fit2,beta=0)
    arr.pred[,,i] <- matrix(y.pred$var1.pred,12,14)
    arr.std[,,i] <- matrix(sqrt(pmax(y.pred$var1.var,0)),12,14)
    if(i==1) sk.all <- as.data.frame(y.pred)[,1:2]
    sk.all <- cbind(sk.all,c(arr.pred[,,i]))
}
sk.all2 <- sk.all[sk.all$lat > 19 & sk.all$long > 99,]
sk.all1 <- sk.all2[-c(1,2,3,4,6,7,8,9,10,11,12,13,14,15,18,19,20,21,22,32,33,44,45,55,56,66,67,77,78,79,88,89,98,99,100,110,111,112,113,114,121,122,123,124,125,126,132,133,134,135,136,137,138,139,140,141,143 ),]
tempe_all <- sk.all1

#Quantile mapping
library(EnvStats)
library(fGarch)
library(ggplot2)

year4<- year2 - 1349
haa1<-haa[,c(year4)]
y2 <- array(0, dim = nrow(data1))

for (i in 1:nrow(data1)) {
  z.mean<-mean(as.numeric(tempe_all[i,-c(1,2)]))

  z.sd<-sd(as.numeric(tempe_all[i,-c(1,2)]))
  
  y.snorm<-snormFit(as.numeric(haa1[i,]))
  
  z<-pnorm(as.numeric(tempe_all[i,-c(1,2)]), z.mean,z.sd)

  y<-qsnorm(z, y.snorm$par[1], y.snorm$par[2], y.snorm$par[3])
  
  y2[i]<- qsnorm(pnorm(0, z.mean,z.sd), y.snorm$par[1], y.snorm$par[2], y.snorm$par[3])
  
}

dff<-data.frame(long = c(tempe_all$long) , lat = c(tempe_all$lat) , y2 = y2)

#Figure9

#jpeg("~/Downloads/value0-v1.png",width=5,height=4 , res = 300, units = "in")
    print(
      ggplot(dff,aes(long,lat)) +
      geom_point(aes(colour=y2),cex=8.8, shape=15) +
      coord_map(xlim=c(98,124.5),ylim=c(18,42.5)) + 
      scale_color_gradientn(colors = c("blue", "cyan", "green", "yellow", "red"), limits=c(-6,25),na.value="transparent", guide="colourbar") +
      #scale_colour_gradientn(colours=rev(brewer.pal(n=9,name='RdBu')),
      #                       limits=c(-6,25),na.value="transparent",
      #                       guide="colourbar") +
      borders(database="world",xlim=c(76,126),ylim=c(18,45),fill=NA,colour="grey50") +
      theme(text=element_text(size=15),legend.title=element_blank(),
        legend.position="right")
    )
#    dev.off()
