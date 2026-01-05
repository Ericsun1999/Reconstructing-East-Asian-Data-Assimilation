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
  
