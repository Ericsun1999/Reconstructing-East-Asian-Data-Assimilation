here::i_am("Code/Supplementary/FigureS2.R")

library(ggplot2)
library(readxl)

ghcnn<-read_excel("~/Downloads/GHCNv4.xlsx",col_type = c("skip","numeric","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip"))

#Beijing
ghcn<-ghcnn[1:161,]
ghcn<-ghcn[-c(103,104),]
df87 <- read.csv("~/Downloads/DA3/tempBv5.csv")

#Shanghai
#ghcn<-ghcnn[443:503, ]
#df87 <- read.csv("~/Downloads/DA3/tempSv5.csv")

#Hong Kong
#ghcn<-ghcnn[839:895, ]
#df87 <- read.csv("~/Downloads/DA3/tempHv5.csv")


ghcnave<-array(0,dim = nrow(ghcn))
for (i in 1:nrow(ghcn)) {
  ghcnave[i]<-rowMeans(ghcn[i,-1],na.rm = T)/100
}
na_per_row <- apply(ghcn, 1, function(x) sum(is.na(x)))
aa<-which(na_per_row < 1)

  beta <- df87$beta
  df87$predicted1 <- (df87$predicted*beta + df87$alpha)
  df92<-data.frame(year=ghcn[aa,1], temperature=ghcnave[aa])
  colnames(df92)<- c("year", "temperature")
  
  df96<- df87
  
  # Beijing year 1837
  # Shanghai year 1847
  # Hong Kong year 1854
  df96_sub <- df96[which(df96$year > 1837), ]
  
  
  df96_sub$se <- sqrt(df96_sub$sigmasmooth2)   
  df96_sub$lo <- df96_sub$predicted - 1.96 * df96_sub$se
  df96_sub$hi <- df96_sub$predicted + 1.96 * df96_sub$se
  
# FigureS2 (a)-(c)
                    
#  jpeg("~/Downloads/time reach lme assimilated ghcn3.png",width=6,height=3, res=300, units = "in")
      print(
      ggplot() +
  geom_line(data = df96_sub,aes(x=year,y=predicted, color = "Assimilated",
                                         linetype = "Assimilated"),
            linewidth = 0.7, alpha = 0.75 ) +
  geom_line(data = df96_sub,aes(x=year,y=predicted1, color = "REACHES", 
                                         linetype = "REACHES"),
            linewidth = 0.7, alpha = 0.75 ) +
  geom_line(data = df96_sub,aes(x=year,y=LME,color="LME",
                                         linetype = "LME"),
            linewidth = 0.7, alpha = 0.75 ) +
  geom_point(data = df92[c(1:47),],aes(x=year,y=temperature,color="GHCN"), cex = 0.7 ) +
  geom_ribbon(
    data = df96_sub,
    aes(x = year, ymin = lo, ymax = hi),
    fill = "grey30", alpha = 0.25
  ) +
  xlab("year")+
  ylab("temperature")+
  scale_color_manual(values = c("Assimilated" = "black", "REACHES" = "firebrick", "LME" = "deepskyblue", "GHCN" = "green3")) +
   scale_linetype_manual(values = c(
    "Assimilated" = "solid",
    "REACHES" = "dashed",
    "LME" = "dotted"
  )) +
  theme(text=element_text(size=12),legend.position="bottom",legend.key.height=unit(-0.5,"cm"),plot.title = element_text(hjust = 0.5)) +
  guides(
    linetype = "none",
    color = guide_legend(
      title = NULL,
      override.aes = list(
        linetype = c("solid", "blank", "dotted", "dashed")
      )
    )
  )
      )
#    dev.off() 
