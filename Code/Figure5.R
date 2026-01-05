# Kriging on 0.5x0.5 grid points

year.start <- range(temp2$year)[1]
year.end <- range(temp2$year)[2]
year.all <- year.start:year.end

year2 <- year.all[c(98,484)]
arr.pred <- arr.std <- array(NA, c(53, 49, length(year2)))
loc <- expand.grid(long = seq(98.25, 124.25, by = 0.5), lat = seq(18.25, 42.25, by = 0.5))
coordinates(loc) <- ~ long + lat
proj4string(loc) <- CRS('+proj=longlat +datum=WGS84')
y2 <- var.fit2$y2

library(RColorBrewer)
library(ggplot2)

#Figure5(a1), 5(b1)

  for(i in 1:length(year2)) {
    temp.dat <- as.data.frame(y2[y2$year==year2[i],])
    #jpeg(paste(getwd(),"/data-temp-",i,".jpg",sep=""),width=5,height=4, res=300, units = "in")
     print(
     ggplot(temp.dat,aes(long,lat)) +
      geom_point(aes(colour=level),cex=2) +
      ggtitle(paste0("Year ",year2[i])) + 
      coord_map(xlim=c(98,124.5),ylim=c(18,42.5)) + 
      scale_colour_gradientn(colours=rev(brewer.pal(n=9,name='RdBu')),limits=c(-2,2),
      na.value="transparent",guide="colourbar") +
      borders(database="world",xlim=c(76,132),ylim=c(19,52),fill=NA,colour="grey50") +
      theme(text=element_text(size=15),legend.title=element_blank(),
        legend.position=c(1.15, 0.7),legend.key.height=unit(2.4,"cm"))
     )
     #dev.off()
  }



#predicted

library(MASS)
library(mvtnorm)


# Define parameters
sigma2 <- vario.fit2$psill[2] - vario.fit2$psill[1]  
alpha <- vario.fit2$range[2]/100  
sigma2_epsilon <- vario.fit2$psill[1]  

# Define cov 
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

#Plot Figure5(b), 6(b)

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
    temp.dat1 <- cbind(as.data.frame(loc@coords)[,1:2],mu=c(arr.pred[,,i]))
    # jpeg(paste(getwd(),"/krige-temp-",i,".jpg",sep=""),width=5,height=4, res=300, units = "in")
    print(
      ggplot(temp.dat1,aes(long,lat)) +
      geom_point(aes(colour=mu),cex=9) +
      ggtitle(paste0("Year ",year2[i])) + 
      coord_map(xlim=c(98,124.5),ylim=c(18,42.5)) + 
      scale_colour_gradientn(colours=rev(brewer.pal(n=9,name='RdBu')),limits=c(-2,2),
      na.value="transparent",guide="colourbar") +
      borders(database="world",xlim=c(76,132),ylim=c(19,52),fill=NA,colour="grey50") +
      theme(text=element_text(size=15),legend.title=element_blank(),
        legend.position="right",legend.key.height=unit(2.2,"cm"),legend.spacing.y = unit(-3,"cm"))
    )
    #dev.off()

  }


