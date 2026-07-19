#Figure6abcd 7abcd

#Kriged REACHES
tempe_all_data <- read.csv("./tempe_all_v3.csv", header=FALSE)
year3<- as.integer(tempe_all_data[1,-c(1,2)]) 
tempe_all<-tempe_all_data[c(2:4),]

#Data1 Hong Kong
#Data2 Shanghai
#Data3 Beijing

Data1 <- read.csv("./d1.csv", row.names=1)
Data2 <- read.csv("./d2.csv", row.names=1)
Data3 <- read.csv("./d3.csv", row.names=1)

#Beijing
DData<-Data3 #Plot Shanghai or hong kong just change it to Data2 or Data1
loc <- 3  #Shanghai or hong kong just change 3 to 2 or 1

tempe_use<-tempe_all[loc,] # Shanghai or hong kong just change 3 to 2 or 1
nu1 <- (read.csv("./tempe_all_std.csv"))
nu <- nu1[loc, -c(1,2)] 

haa<- (DData[,-c(1,2)]) - 273.15
haave <- colMeans(haa)
haave2 <- as.numeric(t(haa))

#QM

library(ggplot2)
library(np)

# ---- Input you need to prepare ----
# haave     : vector of observed LME temperatures (historical sample)
# yhat  : vector of kriging predictions Ŷ_t for each year
# nu    : vector of kriging RMSE (ν_t), one per year
  
yhat <- as.numeric(as.matrix(tempe_use[,-c(1,2)]))
nu <- as.numeric(as.matrix(nu))

# ---- Estimate F_X ----
Fx_hat <- function(q, fx_bw) {
  fx_ob <- npudist(bws = fx_bw, edat = data.frame(x = q))
  # Return estimated CDF values for given q using npudist
  fitted(fx_ob) 
}

# ---- Estimate F_Y ----
FY_hat <- function(y, yhat, nu) {
  # Ensure numeric vectors
  yhat <- as.numeric(yhat)
  nu   <- as.numeric(nu)
  stopifnot(length(yhat) == length(nu))
  # Avoid division by zero
  nu_safe <- pmax(nu, 1e-8)
  # Vectorized computation
  sapply(y, function(yy) {
    mean(pnorm((yy - yhat) / nu_safe))
  })
}

Fx_inv <- function(u, qx_min, qx_max, fx_bw) {
  # u must be between 0 and 1
  u <- pmin(pmax(u, 0 + 1e-8), 1 - 1e-8)
  sapply(u, function(uu) {
    f <- function(q) Fx_hat(q, fx_bw) - uu
    uniroot(f, interval = c(qx_min, qx_max), tol = 1e-6)$root
  })
}

qmapping <- function(yhat = yhat, std = nu, x = haave2, ymax = 1, ymin = -2, xhat = tempe_use[,-c(1,2)]){
  # ---- Estimate F_X using npudist with bandwidth chosen by LS-CV ----
  fx_bw <- npudistbw(dat = c(haave2))       # select bindwidth h via LS-CV
  
  # ---- Inverse of F_X (quantile function) using uniroot ----
  qx_min <- min(x) - 5 * sd(x)
  qx_max <- max(x) + 5 * sd(x)
  
  y_seq <- seq(ymin, ymax, length.out = 150)
  FY_values <- FY_hat(y_seq, yhat, std)
  FX_values <- Fx_inv(FY_values, qx_min, qx_max, fx_bw)
  
  ycorrected <- Fx_inv(FY_hat(as.numeric(xhat), yhat, std), qx_min, qx_max, fx_bw)
  
  list(y_seq = y_seq, FX_values = FX_values, ycorrected = ycorrected)
}

qm <- qmapping(yhat, nu, x = haave2, 1, -2, tempe_use[,-c(1,2)])

ycorrected <- qm$ycorrected


# Figure6abcd

df111<-data.frame(reach=as.numeric(tempe_use[,-c(1,2)]), 
                  year =  c(year3))

#Figure7(a)
 #jpeg("~/Downloads/reachtime.png",width=6,height=4 , res = 300, units = "in")
      print(
      ggplot(data = df111, aes(x = year, y = reach)) +
  geom_line(color = "darkorchid1") +
  theme(text=element_text(size=19),legend.position="right",legend.key.height=unit(1.5,"cm"),plot.title = element_text(hjust = 0.5)) +
  xlab("year") +
  ylab("level") 
      )
   #dev.off()

#Figure7(c)   
  #jpeg("~/Downloads/Reachtemp.png",width=6,height=3.5 , res = 300, units = "in")
      print(
  ggplot(df111, aes(x = reach)) +
  geom_histogram(aes(y = after_stat(density)), 
                 breaks = seq(-2.5, 1.5, by = 0.2),
                 binwidth = 0.1,  
                 fill = "darkorchid1", 
                 color = "white") +
  geom_density(color = "black", size = 0.4) +
  theme(text=element_text(size=18),legend.position="right",legend.key.height=unit(1.5,"cm"),plot.title = element_text(hjust = 0.5)) +
  xlab("temperature") +
  ylab("density")
  )
 #dev.off()

  df<-data.frame(z=qm$y_seq, y1 = qm$FX_values)

#Figure7 (b)
  #jpeg("~/Downloads/fyinverse.png",width=6,height=4 , res = 300, units = "in")
      print(
      ggplot(data = df,aes(x=z,y=y1)) +
  geom_smooth(method = 'loess', formula = 'y ~ x', se = F) +
  xlab("index")+
  ylab("temperature")+
  #ggtitle("Transformation function") +
  theme(text=element_text(size=19),legend.position="right",legend.key.height=unit(1.5,"cm"),plot.title = element_text(hjust = 0.5)) 
      )
    #dev.off()   
      
  dsnorm.fit2 <- data.frame(temper=c(haave2))

#Figure7 (d)
#jpeg("~/Downloads/snormpdfhist.png",width=6,height=3.5 , res = 300, units = "in")
      print(
ggplot(data = dsnorm.fit2, aes(x=temper)) +
  geom_histogram(aes(y = after_stat(density)),breaks = seq(9, 14, by = 0.1),binwidth = 0.1,colour="white",fill="deepskyblue") +
  geom_density(color = "black", size = 0.4) +
  #ggtitle("Skew-normal pdf on LME") +
  theme(text=element_text(size=18),legend.position="right",legend.key.height=unit(1.5,"cm"),plot.title = element_text(hjust = 0.5)) +
  xlim(9,14)+
  xlab("temperature")
  )
# dev.off()

df7<-data.frame(YEAR=c(year3),y= ycorrected)
df7<-cbind(df7,t(as.vector("REACH")))
colnames(df7)<-c("YEAR","temperature","type")
dsnorm.fit1 <- data.frame(temper=c(haave)[(year3-1350)], year = c(year3))

#Figure8(bcd)

#jpeg("~/Downloads/time reach lme corrected.png",width=6,height=3, res=300, units = "in")
      print(
      ggplot() +
  geom_line(data = df7[c(1:524),],aes(x=YEAR,y=temperature),color="palevioletred1", alpha = 0.75 ) +
  geom_smooth(data = df7[c(1:524),],aes(x=YEAR,y=temperature),method = 'loess', formula = 'y ~ x', se = F, size=1.1,color="firebrick", linetype = "dashed", span = 0.75) +
  geom_line(data = dsnorm.fit1[1:524,], aes(x = year, y = temper), color = "skyblue1", alpha = 0.75) +
  geom_smooth(data = dsnorm.fit1[1:524,], aes(x = year, y = temper), method = 'loess', formula = 'y ~ x', se = F, size=1.1,color="deepskyblue", linetype = "dashed") +
  xlab("year")+
  ylab("temperature")+
  theme(text=element_text(size=12),legend.position="bottom",legend.key.height=unit(-0.5,"cm"),plot.title = element_text(hjust = 0.5),  panel.background = )
      )
#    dev.off()  
    
df8<-data.frame(year = c(year3), REACHES = ycorrected, std = nu)

#Figure8(a)
#jpeg("~/Downloads/bPredicted.png",width=6,height=3, res=300, units = "in")
      print(
      ggplot(data = df8, aes(x = year, y = REACHES)) +
      geom_line(color = "red") +
      theme(text=element_text(size=11),legend.position="right",
            legend.key.height=unit(1.5,"cm"),
            plot.title = element_text(hjust = 0.5)) +
      geom_ribbon(aes(ymin = REACHES - std,
                      ymax = REACHES + std),
                  alpha = 0.5, fill = "grey3") +
      xlab("year") +
      ylab("temperature") 
      )
#    dev.off()
