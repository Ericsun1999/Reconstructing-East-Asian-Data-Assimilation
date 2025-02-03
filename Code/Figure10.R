#This code should be executed after running Figures7-9.R

df7<-data.frame(YEAR=c(year3),y=y)
df7<-cbind(df7,t(as.vector("REACH")))
colnames(df7)<-c("YEAR","temperature","type")
dsnorm.fit1 <- data.frame(temper=c(haave)[(year3-1350)], year = c(year3))

#Figure10, Time series plots of Celsius-scaled REACHES data from 1368 to 1911
# (Edit line 19-21 for Shanghai and Hong Kong)
#jpeg("~/Downloads/time reach lme63.png",width=6,height=3, res=300, units = "in")
      print(
      ggplot() +
  #geom_point(cex=1.2,color="palevioletred1") +
  geom_line(data = df7[c(1:524),],aes(x=YEAR,y=temperature),color="palevioletred1", alpha = 0.75 ) +
  geom_smooth(data = df7[c(1:524),],aes(x=YEAR,y=temperature),method = 'loess', formula = 'y ~ x', se = F, size=1.1,color="firebrick", linetype = "dashed", span = 0.75) +
  geom_smooth(data = dsnorm.fit1[1:524,], aes(x = year, y = temper), method = 'loess', formula = 'y ~ x', se = F, size=1.1,color="deepskyblue", linetype = "dashed") +
  xlab("year")+
  ylab("temperature")+
  ylim(10.4,12.4)+ #Beijing
  #ylim(16.45, 18.45)+ #Shanghai
  #ylim(23.05,24.45)+ #Hong Kong
  theme(text=element_text(size=12),legend.position="bottom",legend.key.height=unit(-0.5,"cm"),plot.title = element_text(hjust = 0.5),  panel.background = )
      )
#    dev.off()  
