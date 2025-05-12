#Load data

library(readxl)
temperature<- read_excel("./temperature index value.v1.xlsx",col_type = c("skip","skip","numeric","numeric","skip","skip","skip","skip","skip","numeric","numeric","skip","skip"))
colnames(temperature)<-c("level","year","long","lat")

library(ggplot2)

#Figure 1(a), Annual counts of temperature records in the REACHES dataset from 1368 to 1911

#jpeg("./anntemp.png",width=4,height=2.5 , res = 300, units = "in")
    print(
       ggplot(temperature, aes(x=year)) +
    geom_histogram(binwidth = 1,color ="red",fill="red") +
    theme(text = element_text(size=15))
    )
#    dev.off()

#Figure 1(b), Frequency of each temperature category in the REACHES dataset

#jpeg("./temphisto.png",width=4,height=2.5 , res = 300, units = "in")
    print(
       ggplot(temperature,aes(x=level)) +
    geom_histogram(bins=9) +
    theme(text = element_text(size=15))
    )
 #   dev.off()

# Take mode for events at duplicated locations and rearrange the data by year

library(dplyr)

temp2 <-temperature%>% 
    group_by(year,long,lat) %>%  
    summarise(level = mean(level))

plot_originREACHES_map <- function(data, text_size = 17, legend_height_cm = 3.5, legend_spacing_cm = -3) {
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
      legend.position = "right",
      legend.title = element_blank(),
      legend.key.height = unit(legend_height_cm, "cm"),
      legend.spacing.y = unit(legend_spacing_cm, "cm")
    )
}

#Figure2, A map of temperature event levels
  
#  jpeg("~/Downloads/recordstem1.png",width=6,height=5.5, res=300, units = "in")
      print(
      plot_originREACHES_map(temp2)
      )
#    dev.off() 

