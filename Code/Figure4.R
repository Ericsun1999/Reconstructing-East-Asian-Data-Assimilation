#Load data

  library(readxl)
  temperature<- read_excel("./temperature index value.v1.xlsx",col_type = c("skip","skip","numeric","numeric","skip","skip","skip","skip","skip","numeric","numeric","skip","skip"))
  colnames(temperature)<-c("level","year","long","lat")

  library(dplyr)
  
  temp2 <-temperature%>% 
    group_by(year,long,lat) 

# Compute empirical variograms based on all the data at various years
# Semi-variogram model fitting

  library(sp)
  library(spacetime)
  library(zoo)
  library(gstat)
  library(mgcv)
  library(ggplot2)
  library(mvtnorm)
  library(MASS)


par_est_initial <- function(temp2 = temp2){
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
  vario.stsdf3 <- subset(vario.stsdf2, dist > 1)
  coordinates(y2) <- ~ long + lat
  proj4string(y2) <- CRS('+proj=longlat +datum=WGS84')
  vario2 <- variogram(level~1,data=y2,width=10)
  ind <- match(as.numeric(rownames(vario.stsdf3))-1,rownames(vario2))
  ind2 <- rep(F,nrow(vario2))
  ind2[ind] <- T
  vario2 <- vario2[ind2,]
  if(sum(vario.stsdf3$dist==0)>0) vario2[,1:3] <- vario.stsdf3[ind[-1],1:3] else vario2[,1:3] <- vario.stsdf3[ind,1:3]
  vario.fit2 <- fit.variogram(vario2, vgm(model="Exp",nugget=NA), fit.kappa=T,fit.method=2)
  
  list(psill = vario.fit2$psill, range = vario.fit2$range[2], y2 =y2)
}

var.fit2 <- par_est_initial(temp2)


#Create Var(Y+\epsilon) v.s. Var(Z)
#True Var(\epsilon) v.s. Round Var(\epsilon)
#Calculate \alpha

library(mgcv)
library(ggplot2)
library(mvtnorm)

# Replace y_target with the y-value you are interested in
findgam <- function(y,z, y_target){
  df <- data.frame(vary = y, varz = z)
  # Fit the GAM model
  gam_fit <- gam(varz ~ s(vary), data = df)
  # Generate a range of x-values for prediction
  x_values <- seq(min(df$vary), max(df$vary), length.out = 1000)
  # Predict the y-values for the generated x-values using the fitted GAM model
  predicted_y <- predict(gam_fit, newdata = data.frame(vary = x_values))
  # Find the index of the x-value that gives the closest y-value to y_target
  closest_index <- which.min(abs(predicted_y - y_target))
  # Retrieve the corresponding x-value
  x_closest <- x_values[closest_index]
  # Output the result
  print(paste("The x-value corresponding to y =", y_target, "is:", x_closest))
  x_closest
}

par_estimation_calibrate<- function(psill1, psill, range, n1=500, n2=200, n3=40, plot = F){
  #Var(Y+\epsilon) v.s. Var(Z)
  vary <- array(0,dim = n1)
  varz <- array(0,dim = n1)

  for (i in 1:n1) {
    k=2*i/n1
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
  
  x_closest <- findgam(vary, varz, psill[1]+psill[2])
  df <- data.frame(vary = vary, varz = varz)
  
  if (plot == T){
    #jpeg("~/Downloads/varf1.png",width=6,height=6, res = 300, units = "in")
      print(
      ggplot(df, aes(x=vary, y=varz)) +
  geom_point(aes(x=x_closest, y= psill[1]+psill[2]), colour="red", cex=3) +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=25)) +
  xlab(expression(sigma[Y]^2 + sigma[epsilon]^2)) +
  ylab(expression(f[1]~"("~sigma[Y]^2 + sigma[epsilon]^2~")")) +
  ylim(0, 1.42) +
  xlim(0, 1.42) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  geom_abline(intercept = psill[1]+psill[2], slope = 0, linetype = "longdash", color = "red")+
  geom_smooth(method = 'gam', formula = y ~ s(x), se = F, size=0.7,color="firebrick", linetype = "solid")
      )
  #dev.off()
  }
  
  #Var(\epsilon) v.s. Round Var(\epsilon)
  sigmaepsilon <- array(0,dim = n2)
  sigmaepsilon1 <- array(0,dim = n2)

  for (i in 1:n2) {
    k=i*x_closest/n2
    var=x_closest-k
    y<-rnorm(10000, sd = sqrt(var))
    e<-rmvnorm(10000, sigma = matrix(c(x_closest-var,0,0,x_closest-var), ncol=2))
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
  
  x_closest1 <- findgam(sigmaepsilon, sigmaepsilon1, psill1)
  df1 <- data.frame(sigmaepsilon = sigmaepsilon, sigmaepsilon1 = sigmaepsilon1)
  
  if (plot == T){
  #jpeg("~/Downloads/varf2.png",width=6,height=6, res = 300, units = "in")
      print(
      ggplot(df1, aes(x=sigmaepsilon, y=sigmaepsilon1)) +
  geom_point(aes(x=x_closest1, y=psill1), colour="red", cex=3) +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=25)) +
  xlab(expression(sigma[epsilon]^2)) +
  ylab(expression(f[2]~"("~sigma[epsilon]^2~")")) +
  ylim(0, 0.9) +
  xlim(0, 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  geom_abline(intercept = psill1, slope = 0, linetype = "longdash", color = "red")+
  geom_smooth(method = 'gam', formula = y ~ s(x), se = F, size=0.7,color="firebrick", linetype = "solid")
      )
  #dev.off() 
  }

  alpha2 <- array(0,dim = n3+1)
  alphastar <- array(0,dim = n3+1)

  # Generate the grid points
  a<-array(0,dim = c(2,2))
  cz<- array(0,dim = c(241,151))

  var = x_closest - x_closest1

  for (aaa in 1:151) {
    aa = 95 + 5*aaa
    for (hh in 0:240) {
    a<- array(c(0,0,0,hh*5),dim = c(2,2))
    sigma22 <- matrix(0,nrow = 2,ncol = 2) #10*10
    for (k in 1:2) {
      for (j in 1:2) {
        sigma22[k,j]=var*exp(-sqrt(t(a[,k]-a[,j]) %*% (a[,k]-a[,j]))/aa)
        if (k == j){
          sigma22[k,j] = sigma22[k,j] + x_closest1
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

  for (i in 0:n3) {
    #k = 100 + i*10
    alpha1 = 100 + i*10
    sum2 = 100000000
    kk=kk+1
  for (aaa in kk:151) {
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
    if (kk > 150) kk=150
    np1 = np
    alphastar[i+1] = np
    alpha2[i+1] = alpha1
  }
  x_closest2 <- findgam(alpha2, alphastar, range)
  df2 <- data.frame(alphastar = alphastar, alpha2 = alpha2)
  
  if (plot == T) {
  #jpeg("~/Downloads/varf3.png",width=5,height=6, res = 300, units = "in")
      print(
      ggplot(df2, aes(x= alpha2, y= alphastar)) +
  geom_point(aes(x=x_closest2, y=range), colour="red", cex=3) +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=27)) +
  xlab(expression(alpha)) +
  ylab(expression(f[3]~"("~alpha~")")) +
  ylim(90, 810) +
  xlim(90, 600) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  geom_abline(intercept = range, slope = 0, linetype = "longdash", color = "red")+
  geom_smooth(method = 'gam', formula = y ~ s(x), se = F, size=0.7,color="firebrick", linetype = "solid")
      )
  #dev.off() 
  }
  
  list(psill1 = x_closest1, psill2 = x_closest - x_closest1, range = x_closest2)
  }
  
#Figure4, Calibration functions for model parameters
set.seed(10)
vario.fit2 <-par_estimation_calibrate(var.fit2$psill[1], var.fit2$psill, var.fit2$range, 500, 200, 40, plot = T)

