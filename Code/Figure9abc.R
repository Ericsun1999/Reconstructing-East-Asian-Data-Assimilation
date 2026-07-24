here::i_am("Code/Figure9abc.R")

#For Shanghai, replace the last character with "S" 
#For Hong Kong, replace the last character with "H"

dfc <- read.csv("./mtB.csv")
df2c <- read.csv("./muB.csv")
df3c <- read.csv("./rtB.csv")

#Figure9a, Time series plot of the estimated parameter M_t from years 1368 to 1911
 #jpeg("Figure8a.png",width=6,height=4, res=300, units = "in")
      print(
      ggplot(data = dfc,aes(x=year, y=value)) +
  geom_line(aes(color = coefficient, linetype = coefficient)) + 
  scale_color_manual(values = c(  "steelblue","darkred")) +
  xlab("year")+
  ylab(expression(M))+
  theme(text=element_text(size=12),plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title = element_blank()) +
  scale_linetype_manual(values=c("dotted","solid"))
      )
    #dev.off() 

#Figure9b, Time series plot of the estimated parameter µ_t from years 1368 to 1911
 #jpeg("~/Downloads/mu1.png",width=6,height=4, res=300, units = "in")
      print(
      ggplot(data = df2c,aes(x=year, y=value)) +
  geom_line(aes(color = coefficient, linetype = coefficient)) + 
  scale_color_manual(values = c(  "steelblue","darkred")) +
  xlab("year")+
  ylab(expression(mu))+
  theme(text=element_text(size=12),plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title = element_blank()) +
  scale_linetype_manual(values=c("dotted","solid"))
      )
 #   dev.off() 

 #Figure9c, Time series plot of the estimated parameter r_t^2 from years 1368 to 1911      
  # jpeg("~/Downloads/rt1.png",width=6,height=3, res=300, units = "in")
      print(
      ggplot(data = df3c,aes(x=year, y=value)) +
  geom_line(aes(color = coefficient, linetype = coefficient)) + 
  scale_color_manual(values = c(  "steelblue","darkred")) +
  xlab("year")+
  ylab(expression(r^2))+
  theme(text=element_text(size=12),plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title = element_blank()) +
  scale_linetype_manual(values=c("dotted","solid"))
      )
  #  dev.off() 
