#This code should be executed after running Figures7-9.R

fz <-  dnorm(z.mean, z.mean, z.sd)
fw <-  dsnorm(y.snorm$par[1], mean = y.snorm$par[1],
              y.snorm$par[2],y.snorm$par[3])

fwfz<- (fz/fw)^2

#For Shanghai and Hong kong change mtB, muB, rtB to mtS, muS, rtS & mtH, muH, rtH
  dfc <- read.csv("./mtB.csv")
  df2c <- read.csv("./muB.csv")
  df3c <- read.csv("./rtB.csv")

  sigma2<-array(0,524)

#For Shanghai and Hong kong change Bstd to BstdS & BstdH
  sigma2<- read.csv("./Bstd.csv")
 
  sigma3 = sigma2^2 *fwfz

a<-array(0,523)
b<-array(0,525)
c<-array(0,525)
d<-array(0,525)
kt<-array(0,525)
jt<-array(0,525)
e<-array(0,525)

a[1]=df2c$value[19]
b[1]=df2c$value[20]
c[1]=df3c$rt[19]

mumu <- df2c$value[year3-1367]
mtmt <- dfc$value[year3-1368]
rtrt <- df3c$rt[year3-1367]

for (i in 1:522) {
  b[i+1] = mumu[i+1] - mtmt[i]*mumu[i] + mtmt[i]*a[i]
  d[i] = rtrt[i+1] + mtmt[i]*mtmt[i]*c[i]
  kt[i+1] = d[i]/(d[i] + sigma3[i+1,1])
  a[i+1] = b[i] + kt[i+1]*(y[i] - b[i])
  c[i+1] = (1 - kt[i+1])*d[i]
}

pred = a

for (i in 1:522) {
  jt[523-i] = c[523-i]*mtmt[523-i]/d[523-i]
  pred[523-i] = a[523-i] + jt[523-i]*(pred[524-i] - b[524-i])
  e[523-i] = c[523-i] + jt[523-i]*jt[523-i]*(e[524-i]- d[523-i])
}
    
year4<- year3[1:523]

df87<-data.frame(year = c(year4), predicted = pred, sigma = c[1:523], REACHES = y[1:523], LME = haave[year4-1350])

#Figure13, Predicted temperatures in Beijing from 1368 to 1911
#jpeg("~/Downloads/Predicted.png",width=6,height=3, res=300, units = "in")
      print(
      ggplot(data = df87, aes(x = year, y = predicted)) +
      geom_line(color = "red") +
      theme(text=element_text(size=11),legend.position="right",
            legend.key.height=unit(1.5,"cm"),
            plot.title = element_text(hjust = 0.5)) +
      geom_ribbon(aes(ymin = predicted - sqrt(sigma),
                      ymax = predicted + sqrt(sigma)),
                  alpha = 0.5, fill = "grey3") +
      xlab("year") +
      ylim(10.4,12.4) +
      ylab("temperature") 
      )
#    dev.off()
    
