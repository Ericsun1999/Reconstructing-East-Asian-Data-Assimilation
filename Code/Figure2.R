here::i_am("Code/Figure2.R")

#Load data

library(readxl)
temperature<- read_excel("./temperature index value.v1.xlsx",col_type = c("skip","skip","numeric","numeric","skip","skip","skip","skip","skip","numeric","numeric","skip","skip"))
colnames(temperature)<-c("level","year","long","lat")

library(ggplot2)
library(dplyr)
  
#Figure2b, Empirical frequencies of the temperature categories

#jpeg("./Figure2(b).png",width=4,height=2.5 , res = 300, units = "in")
    print(
       ggplot(temperature,aes(x=level)) +
    geom_histogram(bins=9) +
    theme(text = element_text(size=15))
    )
 #   dev.off()

year_count <- temperature %>%
  count(year)
  
#Figure2a, Annual counts of temperature records in the REACHES dataset

#jpeg("./Figure2(a).png",width=4,height=2.5 , res = 300, units = "in")
    print(
       ggplot(year_count, aes(x = year, y = n)) +
  geom_col(fill = "red", color = "red") +
  geom_smooth(method = "loess", span = 0.3, color = "black", linewidth = 1, se = FALSE) +
  theme(text = element_text(size = 15)) +
  labs(y = "count")
    )
#    dev.off()



# Take mode for events at duplicated locations and rearrange the data by year

library(dplyr)

temp2 <-temperature%>% 
    group_by(year,long,lat) 

library(RColorBrewer)

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

#Figure2c, Spatial distribution of temperature records across East Asia

  #jpeg("./Figure2(c).png",width=5,height=4, res=300, units = "in")
      print(
      plot_originREACHES_map(temp2)
      )
    #dev.off() 

