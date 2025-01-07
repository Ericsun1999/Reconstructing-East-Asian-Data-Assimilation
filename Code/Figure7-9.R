#Load data

#Kriged REACHES
tempe_all_data <- read.csv("./tempe_all.csv", header=FALSE)
year3<- as.integer(tempe_all_data[1,-c(1,2)]) 
tempe_all<-tempe_all_data[c(2:4),]

Data1 <- read.csv("~/Downloads/DA1/a1.csv", row.names=1)
Data2 <- read.csv("~/Downloads/DA1/a2.csv", row.names=1)
Data3 <- read.csv("~/Downloads/DA1/a3.csv", row.names=1)

#Data1 Hong Kong
#Data2 Shanghai
#Data3 Beijing

#Clean data

#Beijing as an example
tempe_use<-tempe_all[3,] # Shanghai or hong kong just change 3 to 2 or 1
Data<-Data3 #Shanghai or hong kong just change it to Data2 or Data1

temp1=array(0,dim = 12)

haa<-matrix(0, nrow = 13, ncol = 600)

for (k in 1:13) {
  for (i in 1:600) {
    for (j in 1:12) {
      temp1[j]=Data[k,-10+12*i+j]-273
    } 
  haa[k,i]=mean(temp1)
  }
}

haave<-array(0,dim=600)

for (i in 1:600) {
  haave[i]=colMeans(haa)[i]
}

#QM (The figures 7 & 8 here adjust to plot Beijing, but the code still can run for Shanghai & Hong Kong)
library(EnvStats)
library(fGarch)
library(ggplot2)

z.mean<-mean(as.numeric(tempe_use[,-c(1,2)]))

z.sd<-sd(as.numeric(tempe_use[,-c(1,2)]))

y.snorm<-snormFit(c(haave))

dsnorm.fit <- data.frame(temper=c(haave))

#Figure8 (b)
#jpeg("~/Downloads/snormpdfhist.png",width=6,height=3.5 , res = 300, units = "in")
      print(
ggplot(data = dsnorm.fit, aes(x=temper)) +
  geom_histogram(aes(y=..density..),breaks = seq(10, 12.5, by = 0.1),binwidth = 0.1,colour="white",fill="deepskyblue") +
  stat_function(fun = dsnorm, args = list(mean = y.snorm$par[1], sd = y.snorm$par[2], xi=y.snorm$par[3])) +
  #ggtitle("Skew-normal pdf on LME") +
  theme(text=element_text(size=18),legend.position="right",legend.key.height=unit(1.5,"cm"),plot.title = element_text(hjust = 0.5)) +
  xlim(10,12.5)+
  xlab("temperature")
  )
# dev.off()


df111<-data.frame(reach=as.numeric(tempe_use[,-c(1,2)]), 
                  year =  c(year3))

#Figure8 (a)
#jpeg("~/Downloads/Reachtemp.png",width=6,height=3.5 , res = 300, units = "in")
      print(
ggplot(data = df111, aes(x=reach)) +
  geom_histogram(aes(y=..density..),breaks = seq(-2.5, 1.5, by = 0.2),binwidth = 0.1,colour="white",fill="darkorchid1") +
  stat_function(fun = dnorm, args = list(mean = z.mean, sd = z.sd)) +
  #ggtitle("Histogram on the REACHES") +
  theme(text=element_text(size=18),legend.position="right",legend.key.height=unit(1.5,"cm"),plot.title = element_text(hjust = 0.5)) +
  xlab("temperature")
  )
# dev.off()

#Figure7 (a)
 #jpeg("~/Downloads/reachtime.png",width=6,height=4 , res = 300, units = "in")
      print(
      ggplot(data = df111, aes(x = year, y = reach)) +
  geom_line(color = "darkorchid1") +
  theme(text=element_text(size=19),legend.position="right",legend.key.height=unit(1.5,"cm"),plot.title = element_text(hjust = 0.5)) +
  xlab("year") +
  ylab("level") 
      )
#   dev.off()


z<-pnorm(as.numeric(tempe_use[,-c(1,2)]), z.mean,z.sd)

y<-qsnorm(z,y.snorm$par[1],y.snorm$par[2],y.snorm$par[3])
y1<-qsnorm(z,y1.snorm$par[1],y1.snorm$par[2],y1.snorm$par[3])

df<-data.frame(z=as.numeric(tempe_use[,-c(1,2)]),y=y)

#Figure7 (b)
#jpeg("~/Downloads/fyinverse.png",width=6,height=4 , res = 300, units = "in")
      print(
      ggplot(data = df,aes(x=z,y=y)) +
  geom_smooth(method = 'loess', formula = 'y ~ x', se = F) +
  xlab("index")+
  ylab("temperature")+
  #ggtitle("Transformation function") +
  theme(text=element_text(size=19),legend.position="right",legend.key.height=unit(1.5,"cm"),plot.title = element_text(hjust = 0.5)) 
      )
#    dev.off()


df7<-data.frame(YEAR=c(year3),y=y)
df7<-cbind(df7,t(as.vector("REACH")))
colnames(df7)<-c("YEAR","temperature","type")
dsnorm.fit1 <- data.frame(temper=c(haave)[(year3-1350)], year = c(year3))

#Figure9 (Edit line 130-132 for Shanghai and Hong Kong)
#jpeg("~/Downloads/time reach lme63.png",width=6,height=3, res=300, units = "in")
      print(
      ggplot() +
  #geom_point(cex=1.2,color="palevioletred1") +
  geom_line(data = df7[c(1:524),],aes(x=YEAR,y=temperature),color="palevioletred1", alpha = 0.75 ) +
  geom_smooth(data = df7[c(1:524),],aes(x=YEAR,y=temperature),method = 'loess', formula = 'y ~ x', se = F, size=1.1,color="firebrick", linetype = "dashed", span = 0.75) +
  geom_smooth(data = dsnorm.fit1[1:524,], aes(x = year, y = temper), method = 'loess', formula = 'y ~ x', se = F, size=1.1,color="deepskyblue", linetype = "dashed") +
  xlab("year")+
  ylab("temperature")+
  ylim(10.4,12.4)+ #Beijing
  #ylim(16.45, 18.45)+ #Shanghai
  #ylim(23.05,24.45)+ #Hong Kong
  theme(text=element_text(size=12),legend.position="bottom",legend.key.height=unit(-0.5,"cm"),plot.title = element_text(hjust = 0.5),  panel.background = )
      )
#    dev.off()  
    

