#Load data

library(readxl)
temperature<- read_excel("./temperature index value.v1.xlsx",col_type = c("skip","skip","numeric","numeric","skip","skip","skip","skip","skip","numeric","numeric","skip","skip"))
colnames(temperature)<-c("level","year","long","lat")

library(ggplot2)

#Annual counts of temperature records in the REACHES dataset from 1368 to 1911
ggplot(temperature, aes(x=year)) +
    geom_histogram(binwidth = 1,color ="red",fill="red") +
    theme(text = element_text(size=15))

#Frequency of each temperature category in the REACHES dataset
ggplot(temperature,aes(x=level)) +
    geom_histogram(bins=9) +
    theme(text = element_text(size=15))

jpeg("./anntemp.png",width=4,height=2.5 , res = 300, units = "in")
    print(
       ggplot(temperature, aes(x=year)) +
    geom_histogram(binwidth = 1,color ="red",fill="red") +
    theme(text = element_text(size=15))
    )
    dev.off()
  
jpeg("./temphisto.png",width=4,height=2.5 , res = 300, units = "in")
    print(
       ggplot(temperature,aes(x=level)) +
    geom_histogram(bins=9) +
    theme(text = element_text(size=15))
    )
    dev.off()
    

