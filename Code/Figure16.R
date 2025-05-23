#This code should be executed after running Figures13.R

#Get GHCN temperature data, Figure16

library(ggplot2)
library(readxl)

ghcnn<-read_excel("~/Downloads/GHCNv4.xlsx",col_type = c("skip","numeric","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip"))

#Hong Kong
ghcn<-ghcnn[839:895, ]
ghcnave<-array(0,dim = nrow(ghcn))

for (i in 1:nrow(ghcn)) {
  ghcnave[i]<-rowMeans(ghcn[i,-1],na.rm = T)/100
}

na_per_row <- apply(ghcn, 1, function(x) sum(is.na(x)))
aa<-which(na_per_row < 1)
    
    df92<-data.frame(year=ghcn[aa,1], temperature=ghcnave[aa])
    colnames(df92)<- c("year", "temperature")
    df93<-inner_join(df87,df92, by="year")
    m=nrow(df93)
    df94<- df93[1:m,]


fitt1<-lm(df94$temperature~df94$predicted)
fitt2<-lm(df94$temperature~df94$REACHES)
fitt3<-lm(df94$temperature~df94$LME)

#Figure16, Scatter plots between GHCN temperatures and estimates in Hong Kong from 1853 to 1911
#jpeg("~/Downloads/scpregh.png",width=6,height=6, res = 300, units = "in")
      print(
      ggplot(df94, aes(x=predicted, y= temperature))+
      geom_point(cex=2.2)+
      geom_abline(intercept = fitt1$coefficients[1], 
                  slope = fitt1$coefficients[2],color="blue") +
      ylab("GHCN")+
      xlab("Assimilated") +
      theme(text=element_text(size=23))+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                  color = "red") +
      ylim(22, 24.5) +
      xlim(22, 24.5) + 
      annotate("text", x = 22.93, y = 24.35, label = "Cor = 0.32", 
               hjust = 1, vjust = 0, color = "black", size = 9) 
      )
#  dev.off() 
  
    
#jpeg("~/Downloads/screagh.png",width=6,height=6, res = 300, units = "in")
    print(
      ggplot(df94, aes(x=REACHES, y= temperature))+
      geom_point(color = "black",cex=2.2)+
      geom_abline(intercept = fitt2$coefficients[1], 
                  slope = fitt2$coefficients[2],color="blue") +
      ylab("GHCN")+
      theme(text=element_text(size=23))+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                  color = "red") +
      coord_cartesian(ylim = c(22, 24.5), xlim = c(22, 24.5)) +
      annotate("text", x = 22.93, y = 24.35, label = "Cor = 0.17", 
               hjust = 1, vjust = 0, color = "black", size = 9)
      )
    
#    dev.off()
    
#jpeg("~/Downloads/sclmgh.png",width=6,height=6, res = 300, units = "in")
    print(
      ggplot(df94, aes(x=LME, y= temperature))+
      geom_point(color = "black",cex=2.2)+
      geom_abline(intercept = fitt3$coefficients[1] , 
                  slope = fitt3$coefficients[2],color="blue") +
      ylab("GHCN")+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                  color = "red") +
      theme(text=element_text(size=23)) + 
      coord_cartesian(ylim = c(22, 24.5), xlim = c(22, 24.5)) +
      annotate("text", x = 22.93, y = 24.35, label = "Cor = 0.28", 
               hjust = 1, vjust = 0, color = "black", size = 9)
      )
#    dev.off()

