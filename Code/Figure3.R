here::i_am("Code/Figure3.R")

#Data1 Hong Kong
#Data2 Shanghai
#Data3 Beijing

Data1 <- read.csv("./d1.csv", row.names=1)
Data2 <- read.csv("./d2.csv", row.names=1)
Data3 <- read.csv("./d3.csv", row.names=1)

library(ggplot2)

#Beijing
Data <- Data3 #Plot Shanghai or hong kong just change it to Data2 or Data1

haa<- (DData[,-c(1,2)]) - 273
mean_row <- colMeans(haa, na.rm = TRUE)
dfp <- data.frame(t(haa))
dfp$dataave <- as.numeric(mean_row) 
dfp$Year <- as.numeric(1350:1949)
colnames(dfp) <- c(
  paste0("data", 1:13),
  "dataave",
  "Year"
)
    
#Figure3, Yearly temperature time series data from 13 LME simulations

#jpeg("./Figure3.png",width=6,height=3, res=300, units = "in")
      print(
      ggplot(data=dfp, aes(x=Year)) +
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
#  ylim(9,13.5)+ #Beijing
# ylim(13,18) + #Shanghai
# ylim(20,24) + #Hong Kong
  theme(text=element_text(size=12), plot.title = element_text(hjust = 0.5), legend.title = element_text(),legend.position = "none")
      )
#   dev.off()
