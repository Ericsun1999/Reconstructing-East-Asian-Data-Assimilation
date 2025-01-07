#This code should be executed after running Figure7-9.R

library(penalized)
library(tidyverse)

lam1 = 0.03
lam2 = 10
lam3 = 10

mt<-array(0.5,543)

#fused lasso

mu<-haave

rt<-array(0,600)

tte<-matrix(0,13,7)

for (i in 1:600) {
  for (j in 1:13) {
      rt[i] = rt[i] + (haa[j,i]-mu[i])^2
  }
  rt[i]=rt[i]/12
}


for (z in 1:100) {
  #set fused lasso y,x
  yh<-matrix(haa[,(20:562)],nrow = 13*543)
  xh<-matrix(0,nrow = 13*543, ncol = 543)
  
  
    for (i in 1:543) {
    for (j in 1:13) {
      t=13*(i-1)
      xh[t+j,i]=(haa[j,i+20]-mu[i+20])/sqrt(rt[i+20])
    }
  }

  for (i in 1:543) {
    for (j in 1:13) {
      yh[13*(i-1)+j,1]=(yh[13*(i-1)+j,1]-mu[i+19])/sqrt(rt[i+20])
   }
  }
  
  dat1 = data.frame(cbind(yh,xh))
  
  fit3 = penalized(response = dat1[,1], penalized = dat1[,2:544], unpenalized = ~0, lambda1 = lam1, lambda2 = lam2 ,fusedl = T, positive = T, data = dat1, model = c("linear"), maxiter = 300)
  
  fit4 = lm(yh~xh-1)
  
  mt[1:543] = fit3@penalized
  
  yv<-matrix(0,nrow = 13*544, ncol = 1)
  xv<-matrix(0,nrow = 13*544, ncol = 544)
  
  for (i in 1:13) {
    yv[i,1] = (haa[i,19])/sqrt(rt[19])
    xv[i,1] = 1/sqrt(rt[19])
  }

  for (k in 1:543) {
    for (i in 1:13) {
      yv[13*(k)+i,1] = (haa[i,k+20]-mt[k]*haa[i,k+19])/sqrt(rt[k+20])
      xv[13*(k)+i,k] = -mt[k]/sqrt(rt[k+20])
      xv[13*(k)+i,k+1] = 1/sqrt(rt[k+20])      
    }
  }
  dat2 = data.frame(cbind(yv,xv))
  
  fit2 = penalized(response = dat2[,1], penalized = dat2[,2:545], unpenalized = ~0, lambda1 = 0, lambda2 = lam3 ,fusedl = T , data = dat2, model = c("linear"), maxiter = 300)
  
  fit5 = lm(yv~xv-1)
  
  mu[20:563] = fit2@penalized
  
}

dfc<- data.frame(year = c(1368:1910), penalized = mt[1:543], unpenalized = coefficients(fit4))

dfc<-tidyr::gather(dfc,key = "coefficient", value = "value", -year)

dfc[1:543,2]<-"FLasso"
dfc[544:1086,2]<-"ML"


df2c<- data.frame(year = c(1369:1911), penalized = mu[20:562], lme = haave[19:561])

df2c<-tidyr::gather(df2c,key = "coefficient", value = "value", -year)

df2c[1:544,2]<-"FLasso"

write_excel_csv(dfc,"mtB.csv")
write_excel_csv(df2c,"muB.csv")
