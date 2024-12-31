#Load data

library(readxl)
temperature<- read_excel("./temperature index value.v1.xlsx",col_type = c("skip","skip","numeric","numeric","skip","skip","skip","skip","skip","numeric","numeric","skip","skip"))
colnames(temperature)<-c("level","year","long","lat")

# Take mode for events at duplicated locations and rearrange the data by year

library(dplyr)
sum(duplicated(temperature[,c(2,3,4)])) #check for duplicate
  
temp2 <-temperature%>% 
   group_by(year,long,lat) %>%  
   summarise(level = mean(level))


# A histogram of temperature levels
  library(ggplot2)
  
  ggplot(temperature,aes(x=level)) +
    geom_histogram(bins=9) +
    theme(text = element_text(size=15))
  
  ggplot(temperature, aes(x=year)) +
    geom_histogram(binwidth = 1,color ="red",fill="red") +
    theme(text = element_text(size=15))
  
  jpeg("./temphisto.png",width=4,height=2.5 , res = 300, units = "in")
    print(
       ggplot(temperature,aes(x=level)) +
    geom_histogram(bins=9) +
    theme(text = element_text(size=15))
    )
    dev.off()
    
  jpeg("./anntemp.png",width=4,height=2.5 , res = 300, units = "in")
    print(
       ggplot(temperature, aes(x=year)) +
    geom_histogram(binwidth = 1,color ="red",fill="red") +
    theme(text = element_text(size=15))
    )
    dev.off()
