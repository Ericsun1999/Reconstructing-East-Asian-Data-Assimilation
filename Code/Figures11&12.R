dfc <- read.csv("./mtB.csv")
df2c <- read.csv("./muB.csv")

#Figure11, Time series plot of the estimated parameter M_t from years 1368 to 1911
 #jpeg("~/Downloads/mt.png",width=6,height=4, res=300, units = "in")
      print(
      ggplot(data = dfc,aes(x=year, y=value)) +
  geom_line(aes(color = coefficient, linetype = coefficient)) + 
  scale_color_manual(values = c( "darkred", "steelblue")) +
  xlab("year")+
  ylab(expression(M))+
  #ggtitle(bquote(Time~series~plot~of~M[t])) +
  theme(text=element_text(size=12),plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_linetype_manual(values=c("solid","dotted"))
      )
    #dev.off() 

#Figure12, Time series plot of the estimated parameter µ_t from years 1368 to 1911
 #jpeg("~/Downloads/mu1.png",width=6,height=4, res=300, units = "in")
      print(
      ggplot(data = df2c,aes(x=year, y=value)) +
  geom_line(aes(color = coefficient, linetype = coefficient)) + 
  scale_color_manual(values = c("darkred", "steelblue")) +
  xlab("year")+
  ylab(expression(mu))+
  #ggtitle(expression(Time~series~plot~of~mu[t]~and~LME)) +
  theme(text=element_text(size=12),plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_linetype_manual(values=c("solid", "dotted"))
      )
 #   dev.off() 
