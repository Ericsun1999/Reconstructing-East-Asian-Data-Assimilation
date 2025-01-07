
#Data1 Hong Kong
#Data2 Shanghai
#Data3 Beijing

Data1 <- read.csv("./a1.csv", row.names=1)
Data2 <- read.csv("./a2.csv", row.names=1)
Data3 <- read.csv("./a3.csv", row.names=1)


#Beijing
Data<-Data3 #Plot Shanghai or hong kong just change it to Data2 or Data1

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

haave1<-array(haa,dim=600*13)


#ggplot LME 

library(ggplot2)

dfp<-data.frame(
  data1=haa[1,],
  data2=haa[2,],
  data3=haa[3,],
  data4=haa[4,],
  data5=haa[5,],
  data6=haa[6,],
  data7=haa[7,],
  data8=haa[8,],
  data9=haa[9,],
  data10=haa[10,],
  data11=haa[11,],
  data12=haa[12,],
  data13=haa[13,],
  dataave=haave,
  Year=c(1350:1949)
)

#Yearly temperature time series data from 13 LME simulations, Figure3(a) 
#(This show Bejing, for Shanghai and Hong Kong, edit line 80-82)

#jpeg("~/Downloads/lmetime63.png",width=6,height=3, res=300, units = "in")
      print(
      ggplot(data=dfp[19:562,], aes(x=Year)) +
  geom_line(aes(y=data1,col="data1"))+
  geom_line(aes(y=data2,col="data2"))+
  geom_line(aes(y=data3,col="data3"))+
  geom_line(aes(y=data4,col="data4"))+
  geom_line(aes(y=data5,col="data5"))+
  geom_line(aes(y=data6,col="data6"))+
  geom_line(aes(y=data7,col="data7"))+
  geom_line(aes(y=data8,col="data8"))+
  geom_line(aes(y=data9,col="data9"))+
  geom_line(aes(y=data10,col="data10"))+
  geom_line(aes(y=data11,col="data11"))+
  geom_line(aes(y=data12,col="data12"))+
  geom_line(aes(y=data13,col="data13"))+
  geom_line(aes(y=dataave,col="average"))+
  scale_colour_manual("", values = c( "average"="deepskyblue")) +
  ylab("temperature")+
  xlab("year")+
  xlim(1368,1911)+
  ylim(9,13.5)+
# ylim(14,20) + #Shanghai
# ylim(21.5,25) + #Hong Kong
  theme(text=element_text(size=12), plot.title = element_text(hjust = 0.5), legend.title = element_text(),legend.position = "none")
      )
 #  dev.off()
