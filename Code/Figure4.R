#Load data

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
  
  print(vario.fit2)

#Create Var(Y+\epsilon) v.s. Var(Z)

  library(mgcv)
  library(ggplot2)
  set.seed(10)

  n=500 #binwidth

  vary <- array(0,dim = n)
  varz <- array(0,dim = n)


for (i in 1:n) {
  k=2*i/n
  var=k
  y<-rnorm(10000, 0, sqrt(var))
  z<-array(0,dim = 10000)

  for (j in 1:10000) {
    z[j]<- round(y[j], digits = 0)
    if (z[j] > 1) z[j]=1
    if (z[j] < -2) z[j]=-2
  }
  vary[i] <- var(y)
  varz[i] <- var(z)
}


df <- data.frame(vary = vary, varz = varz)

#Figure4(a), calibration functions for model parameters

#jpeg("~/Downloads/varf1.png",width=6,height=6, res = 300, units = "in")
      print(
      ggplot(df, aes(x=vary, y=varz)) +
  geom_point(aes(x=0.885, y=0.761), colour="red", cex=3) +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=25)) +
  xlab(expression(sigma[Y]^2 + sigma[epsilon]^2)) +
  ylab(expression(f[1]~"("~sigma[Y]^2 + sigma[epsilon]^2~")")) +
  ylim(0, 1.02) +
  xlim(0, 1.02) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  geom_abline(intercept = 0.761, slope = 0, linetype = "longdash", color = "red")+
  geom_smooth(method = 'gam', formula = y ~ s(x), se = F, size=0.7,color="firebrick", linetype = "solid")
      )
#  dev.off()


#True Var(\epsilon) v.s. Round Var(\epsilon)
library(mvtnorm)
set.seed(10)

n=200 #binwidth

sigmaepsilon <- array(0,dim = n)
sigmaepsilon1 <- array(0,dim = n)


for (i in 1:n) {
  k=i*0.885/n
  var=0.885-k
  
  y<-rnorm(10000, sd = sqrt(var))
  e<-rmvnorm(10000, sigma = matrix(c(0.885-var,0,0,0.885-var), ncol=2))
  
  z1<-array(0,dim = 10000)
  z2<-array(0,dim = 10000)

  for (j in 1:10000) {
    z1[j]<- round(y[j]+e[j,1], digits = 0)
    if (z1[j] > 1) z1[j]=1
    if (z1[j] < -2) z1[j]=-2
    z2[j]<- round(y[j]+e[j,2], digits = 0)
    if (z2[j] > 1) z2[j]=1
    if (z2[j] < -2) z2[j]=-2
  }
  
  sigmaepsilon[i] <- var(e[,1]-e[,2])/2 
  sigmaepsilon1[i] <- var(z1-z2)/2
  
}

df1 <- data.frame(sigmaepsilon = sigmaepsilon, sigmaepsilon1 = sigmaepsilon1)

#Figure4(b), calibration functions for model parameters

#jpeg("~/Downloads/varf2.png",width=6,height=6, res = 300, units = "in")
      print(
      ggplot(df1, aes(x=sigmaepsilon, y=sigmaepsilon1)) +
  geom_point(aes(x=0.146, y=0.196), colour="red", cex=3) +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=25)) +
  xlab(expression(sigma[epsilon]^2)) +
  ylab(expression(f[2]~"("~sigma[epsilon]^2~")")) +
  ylim(0, 0.9) +
  xlim(0, 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  geom_abline(intercept = 0.196, slope = 0, linetype = "longdash", color = "red")+
  geom_smooth(method = 'gam', formula = y ~ s(x), se = F, size=0.7,color="firebrick", linetype = "solid")
      )
#  dev.off() 

#Calculate \alpha

#sigma2y = 0.885-0.146 = 0.739
#sigma2epsilon = 0.146

set.seed(320)

n=30 #binwidth

alpha2 <- array(0,dim = n+1)
alphastar <- array(0,dim = n+1)

# Generate the grid points
a<-array(0,dim = c(2,2))
cz<- array(0,dim = c(241,121))

var=0.739

for (aaa in 1:121) {
    aa = 95 + 5*aaa  #Different alpha setting
    
    for (hh in 0:240) {  
    a<- array(c(0,0,0,hh*5),dim = c(2,2))
    
    sigma22 <- matrix(0,nrow = 2,ncol = 2) 
    for (k in 1:2) {
      for (j in 1:2) {
        sigma22[k,j]=var*exp(-sqrt(t(a[,k]-a[,j]) %*% (a[,k]-a[,j]))/aa)
        if (k == j){
          sigma22[k,j] = sigma22[k,j] + 0.146
        }
      }
    }
    
    y<-rmvnorm(1000, sigma = sigma22)
    
    z1<-array(0, dim = c(1000,2))
  
    for (k in 1:1000) {
      for (j in 1:2) {
        z1[k,j]<- round(y[k,j], digits = 0)
        if (z1[k,j] > 1) z1[k,j]=1
        if (z1[k,j] < -2) z1[k,j]=-2
      }
    }
    
    cz[hh+1,aaa] <- cov(z1[,1], z1[,2])
    }
}
  
np1 = 0
kk=0

for (i in 0:n) {
  k = 100 + i*10
  var = 0.739
  
  alpha1 = k
  
  sum2 = 100000000
  
  kk=kk+1
  
  for (aaa in kk:121) {
    sum1 = 0
    aa = 95 + 5*aaa
    for (hh in 0:240) {
      sum1 = sum1 + (cz[hh+1, aaa] - var*exp(-(hh*5)/alpha1))^2
    }
    if (sum2 >= sum1){
      sum2 = sum1
      np = aa
    }
  }
  
  kk= ((np)-95)/5
  
  if (kk > 121) {kk=121}
  
  
  np1 = np
  
  
  alphastar[i+1] = np
  alpha2[i+1] = alpha1
}

df2 <- data.frame(alphastar = alphastar, alpha2 = alpha2)

#Figure4(c), calibration functions for model parameters
#jpeg("~/Downloads/varf3.png",width=5,height=6, res = 300, units = "in")
      print(
      ggplot(df2, aes(x= alpha2, y= alphastar)) +
  geom_point(aes(x=299.099, y=426.66), colour="red", cex=3) +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=27)) +
  xlab(expression(alpha)) +
  ylab(expression(f[3]~"("~alpha~")")) +
  ylim(90, 540) +
  xlim(90, 450) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  geom_abline(intercept = 426.66, slope = 0, linetype = "longdash", color = "red")+
  geom_smooth(method = 'gam', formula = y ~ s(x), se = F, size=0.7,color="firebrick", linetype = "solid")
      )
 # dev.off() 
      
