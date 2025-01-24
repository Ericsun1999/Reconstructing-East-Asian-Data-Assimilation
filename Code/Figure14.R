#This code should be executed after running Figures13.R

#Get GHCN temperature data
library(ggplot2)
library(readxl)

ghcnn<-read_excel("./GHCNv4.xlsx",col_type = c("skip","numeric","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip"))

#Beijing
ghcn<-ghcnn[1:161,]
ghcn<-ghcn[-c(103,104),]

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

#Valid Beijing, Figure14
                    
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
      ylim(10.65, 12.65) +
      xlim(10.65, 12.65) +
      annotate("text", x = 12.65, y = 10.65, label = "Cor = 0.45", 
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
      coord_cartesian(ylim = c(10.65, 12.65), xlim = c(10.65,12.65)) +
      annotate("text", x = 12.65, y = 10.65, label = "Cor = 0.50", 
               hjust = 1, vjust = 0, color = "black", size = 9) 
      )
    
#    dev.off()
    
#jpeg("~/Downloads/sclmgh.png",width=6,height=6, res = 300, units = "in")
    print(
      ggplot(df94, aes(x=LME, y= temperature))+
      geom_point(color = "black",cex=2.2)+
      geom_abline(intercept = fitt3$coefficients[1], 
                  slope = fitt3$coefficients[2],color="blue") +
      ylab("GHCN")+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                  color = "red") +
      theme(text=element_text(size=23)) +
      coord_cartesian(ylim = c(10.65, 12.65), xlim = c(10.65,12.65)) +
      annotate("text", x = 12.65, y = 10.65, label = "Cor = 0.33", 
               hjust = 1, vjust = 0, color = "black", size = 9)
      )
#    dev.off()
