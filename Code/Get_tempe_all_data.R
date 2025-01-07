
library(readxl)
temperature<- read_excel("./temperature index value.v1.xlsx",col_type = c("skip","skip","numeric","numeric","skip","skip","skip","skip","skip","numeric","numeric","skip","skip"))
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


# Kriging on 0.5x0.5 grid points
year.start <- range(temp2$year)[1]
year.end <- range(temp2$year)[2]
year.all <- year.start:year.end

ncase1 <- rep(NA, nyear)
for (i in 1:nyear) {
  ncase1[i] <- sum(temp2$year == year.start + i - 1)
  if (ncase1[i] > 0) {
    temp21 <- temp2[temp2$year == (year.start + i - 1), ]
  }
}
m <- 1
year2 <- year.all[which(ncase1 >= m)]

arr.pred <- arr.std <- array(NA, c(53, 49, length(year2)))

loc <- expand.grid(long = seq(98.25, 124.25, by = 0.5), lat = seq(18.25, 42.25, by = 0.5))
coordinates(loc) <- ~ long + lat
proj4string(loc) <- CRS('+proj=longlat +datum=WGS84')

#predicted

library(MASS)
library(mvtnorm)

sigma2 <- vario.fit2$psill[2] - vario.fit2$psill[1] 
alpha <- vario.fit2$range[2]/80
sigma2_epsilon <- vario.fit2$psill[1]

cov_Y <- function(s1, s2) {
  (sigma2) * exp(-sqrt(sum((s1-s2)^2))/alpha)
}

calc_c_Y <- function(s0, locations) {
  n <- nrow(locations)
  c_Y <- numeric(n)
  for (i in 1:n) {
    si <- locations[i, ]
    c_Y[i] <- cov_Y(s0, si)
  }
  return(c_Y)
}

  mu_z <- -pnorm(-1.5, mean = 0, sd = sqrt(sigma2 + sigma2_epsilon))

  calc_Sigma_zz <- function(locations) {
  n <- nrow(locations)
  Sigma_zz <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      si <- locations[i, ]
      sj <- locations[j, ]
      
      
      Sigma_zz[i, j] <- cov_Y(si, sj)
      if (i == j) {
        
        Sigma_zz[i, j] <- Sigma_zz[i, j] + (sigma2_epsilon)
      }
    }
  }
  return(Sigma_zz)
}

  for(i in 1:length(year2)) {
    temp14 <- y2[y2$year==year2[i],]
    temp15 <- rep(0, length(loc@coords[,1]))
    temp16 <- rep(0, length(loc@coords[,1]))
    
    locations<- temp14@coords
    
    Sigma_zz <- calc_Sigma_zz(locations)
    invsigma_zz<- solve(Sigma_zz)
    c_Y <- matrix(0, nrow = length(loc@coords[,1]), ncol = length(temp14$level))
    Mat1<- matrix(mu_z, nrow = 1, ncol = length(temp14$level))
    
    for (j in 1:length(loc@coords[,1])) {
      
      s0 <- as.numeric(loc@coords[j,])
      c_Y[j,] <- calc_c_Y(s0, locations)
      temp15[j]<- t(as.matrix(c_Y[j,])) %*% as.matrix(invsigma_zz) %*% (as.matrix(temp14$level) - t(Mat1) )
      temp16[j] <- cov_Y(s0,s0) - t(c_Y[j,]) %*% as.matrix(invsigma_zz) %*% (c_Y[j,])
    }
    arr.pred[,,i] <- matrix(temp15,53,49)
    arr.std[,,i] <- matrix(temp16,53,49)
    
    if(i==1) sk.all <- as.data.frame(loc@coords)
    sk.all <- cbind(sk.all,c(arr.pred[,,i]))
    
  }
  
#Write csv
  names(sk.all)[-c(1,2)] <- year2
#The point of Hong Kong, Shanghai, Beijing
  write_excel_csv(sk.all[c(456,1425,2316),],"tempe_all.csv")
  
