here::i_am("Code/Supplementary/FigureS1.R")

# A map of temperature event levels for different periods, S1 (a)-(e) 

library(readxl)
library(dplyr)
library(RColorBrewer)

temperature<- read_excel(here::here("Data", "temperature index value.v1.xlsx"),col_type = c("skip","skip","numeric","numeric","skip","skip","skip","skip","skip","numeric","numeric","skip","skip"))
colnames(temperature)<-c("level","year","long","lat")

temp2 <-temperature%>% 
    group_by(year,long,lat) 

tempa<- temp2[which(temp2$year < 1501),]
tempb<- temp2[which(temp2$year > 1500 & temp2$year < 1601),]
tempc<- temp2[which(temp2$year > 1600 & temp2$year < 1701),]
tempd<- temp2[which(temp2$year > 1700 & temp2$year < 1801),]
tempe<- temp2[which(temp2$year > 1800),]


plot_originREACHES_map <- function(data, text_size = 15, legend_height_cm = 2.5, legend_spacing_cm = 30) {
  ggplot(data, aes(long, lat)) +
    borders(database = "world", xlim = c(95, 126), ylim = c(19, 45), fill = NA, colour = "grey30") +
    geom_point(aes(colour = level), cex = 1) +
    coord_map(xlim = c(98, 124.5), ylim = c(19, 42.5)) +
    scale_colour_gradientn(
      colours = rev(brewer.pal(n = 9, name = 'RdBu')),
      limits = c(-2, 2),
      na.value = "transparent",
      guide = "colourbar"
    ) +
    theme(
      text = element_text(size = text_size),
      legend.position = c(1.12,0.61),
      legend.title = element_blank(),
      legend.key.height = unit(legend_height_cm, "cm"),
      #legend.key.spacing.y = unit(legend_spacing_cm, "cm")
    )
}
      
  # Figure S1 (a)-(e)

    #jpeg("~/Downloads/recordstem1368-1500.png",width=5,height=4, res=300, units = "in")
      print(
      plot_originREACHES_map(tempa)
      )
    #dev.off() 
      
  #jpeg("~/Downloads/recordstem1501-1600.png",width=5,height=4, res=300, units = "in")
      print(
      plot_originREACHES_map(tempb)
      )
    #dev.off() 
      
  #jpeg("~/Downloads/recordstem1601-1700.png",width=5,height=4, res=300, units = "in")
      print(
      plot_originREACHES_map(tempc)
      )
    #dev.off() 

  #jpeg("~/Downloads/recordstem1701-1800.png",width=5,height=4, res=300, units = "in")
      print(
      plot_originREACHES_map(tempd)
      )
    #dev.off()       
    
  #jpeg("~/Downloads/recordstem1801-1911.png",width=5,height=4, res=300, units = "in")
      print(
      plot_originREACHES_map(tempe)
      )
    #dev.off() 
```
