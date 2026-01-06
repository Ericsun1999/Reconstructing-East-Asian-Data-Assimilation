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

    
    df87 <- read.csv("~/Downloads/DA3/tempHv5.csv")
    df92<-data.frame(year=ghcn[aa,1], temperature=ghcnave[aa])
    colnames(df92)<- c("year", "temperature")
    df94<-inner_join(df87,df92, by="year")


fitt1<-lm(df94$temperature~df94$predicted)
fitt2<-lm(df94$temperature~df94$REACHES)
fitt3<-lm(df94$temperature~df94$LME)

  
#jpeg("~/Downloads/scpregh.png",width=6,height=6, res = 300, units = "in")
      print(
      ggplot()+
      geom_point(data = df94, aes(x=predicted, y= temperature),cex=2.2)+
      geom_abline(intercept = fitt1$coefficients[1], 
                  slope = fitt1$coefficients[2],color="blue") +
      ylab("GHCN")+
      xlab("Assimilated") +
      theme(text=element_text(size=23))+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                  color = "red") +
      ylim(21.5, 24.5) +
      xlim(21.5, 24.5) + 
      annotate("text", x = 22.93, y = 24.35, label = "Cor = 0.38", 
               hjust = 1, vjust = 0, color = "black", size = 9) 
      )
#  dev.off() 
  
    
#jpeg("~/Downloads/screagh.png",width=6,height=6, res = 300, units = "in")
    print(
      ggplot()+
      geom_point(data = df94, aes(x=REACHES, y= temperature),cex=2.2)+
      geom_abline(intercept = fitt2$coefficients[1], 
                  slope = fitt2$coefficients[2],color="blue") +
      ylab("GHCN")+
      theme(text=element_text(size=23))+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                  color = "red") +
      ylim(21.5, 24.5) +
      xlim(21.5, 24.5) + 
      annotate("text", x = 22.93, y = 24.35, label = "Cor = 0.15", 
               hjust = 1, vjust = 0, color = "black", size = 9)
      )
    
  #  dev.off()
    
#jpeg("~/Downloads/sclmgh.png",width=6,height=6, res = 300, units = "in")
    print(
      ggplot()+
      geom_point(data = df94, aes(x=LME, y= temperature), color = "black",cex=2.2)+
      geom_abline(intercept = fitt3$coefficients[1] , 
                  slope = fitt3$coefficients[2],color="blue") +
      ylab("GHCN")+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                  color = "red") +
      theme(text=element_text(size=23)) + 
      ylim(21.5, 24.5) +
      xlim(21.5, 24.5) + 
      annotate("text", x = 22.93, y = 24.35, label = "Cor = 0.32", 
               hjust = 1, vjust = 0, color = "black", size = 9)
      )
#    dev.off()
