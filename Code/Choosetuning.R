#This code should be executed after running Figure7-9.R

library(penalized)
library(tidyverse)

mt<-array(0.5,543)

#fused lasso

mu<-haave
rt<-array(0,600)
tte<-array(0, c(13, 7, 7, 7))


lam1<-3*10^seq(-3, 3)
lam2<-10^seq(-3, 3)
lam3<-10^seq(-3, 3)


for (ii in 1:7) {
  for (jj in 1:7) {
    for (kk in 1:7) {
      llambda1 <- lam1[ii]
      llambda2 <- lam2[jj]
      llambda3 <- lam3[kk]
      
      for (i in 1:600) {
        for (j in 1:13) {
          rt[i] = rt[i] + (haa[j,i]-mu[i])^2
        }
        rt[i]=rt[i]/12
      }
      
      for (z in 1:13) {
        #z=13
        haaa=haa[-z,]
        for (iter in 1:100) {
          #set fused lasso y,x
        yh<-matrix(haaa[,(20:562)],nrow = 12*543)
        xh<-matrix(0,nrow = 12*543, ncol = 543)
  
        for (i in 1:543) {
          for (j in 1:12) {
            t=12*(i-1)
            xh[t+j,i]=(haaa[j,i+20]-mu[i+20])/sqrt(rt[i+20])
          }
        }

        for (i in 1:543) {
          for (j in 1:12) {
            yh[12*(i-1)+j,1]=(yh[12*(i-1)+j,1]-mu[i+19])/sqrt(rt[i+20])
          }
        }
  
        dat1 = data.frame(cbind(yh,xh))
  
        fit3 = penalized(response = dat1[,1], penalized = dat1[,2:544],
                         unpenalized = ~0, lambda1 = llambda1, 
                         lambda2 = llambda2 ,fusedl = T, positive = T, 
                         data = dat1, model = c("linear"), 
                         maxiter = 300)
  
        mt[1:543] = fit3@penalized
  
        yv<-matrix(0,nrow = 12*544, ncol = 1)
        xv<-matrix(0,nrow = 12*544, ncol = 544)
  
        for (i in 1:12) {
          yv[i,1] = (haaa[i,19])/sqrt(rt[19])
          xv[i,1] = 1/sqrt(rt[19])
        }

        for (k in 1:543) {
        for (i in 1:12) {
        yv[12*(k)+i,1]=(haaa[i,k+20]-mt[k]*haaa[i,k+19])/sqrt(rt[k+20])
        xv[12*(k)+i,k] = -mt[k]/sqrt(rt[k+20])
        xv[12*(k)+i,k+1] = 1/sqrt(rt[k+20])      
        }
        }
        dat2 = data.frame(cbind(yv,xv))
  
        fit2 = penalized(response = dat2[,1], penalized = dat2[,2:545],
                         unpenalized = ~0, lambda1 = 0, 
                         lambda2 = llambda3 ,fusedl = T , data = dat2, 
                         model = c("linear"), maxiter = 300)
  
        mu[20:563] = fit2@penalized
        for (i in 1:600) {
          for (j in 1:12) {
            rt[i] = rt[i] + (haaa[j,i]-mu[i])^2
          }
           rt[i]=rt[i]/11
        }

        }
        
        sum=0
        for (i in 1:543) {
          sum = sum + (((haa[z,i+20]-mu[i+20])-mt[i]*(haa[z,i+19]-mu[i+19]))^2)/(2*rt[i+1])
        }
        tte[z,ii,jj,kk]=sum
      }
      
      
    }
  }
}

# Calculate the mean along the first dimension
mean_array <- apply(tte, c(2, 3, 4), mean)

# Find the minimum value of the resulting mean array
min_value <- min(mean_array)

# Find the indices where the minimum value occurs
min_indices <- which(mean_array == min_value, arr.ind = TRUE)

min_indices

#lam1 = 1e-2
#lam2 = 1e1
#lam3 = 1e1



