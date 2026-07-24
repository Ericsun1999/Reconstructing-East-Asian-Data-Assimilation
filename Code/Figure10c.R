here::i_am("Code/Figure10c.R")

library(ggplot2)
library(readxl)
library(tidyverse)

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

df94_H <- df94 %>%
  rename(GHCN = temperature)


df_long_H <- df94_H %>%
  select(year, GHCN, predicted, REACHES, LME) %>%
  pivot_longer(
    cols = c(predicted, LME, REACHES),
    names_to = "Method",
    values_to = "Estimate"
  ) %>%
  mutate(
    Method = recode(
      Method,
      predicted = "Assimilated",
      LME = "LME",
      REACHES = "REACHES"
    ),
    Method = factor(
      Method,
      levels = c("Assimilated", "LME", "REACHES")
    )
  )


cor_lab_H <- df_long_H %>%
  group_by(Method) %>%
  summarise(
    cor_val = cor(GHCN, Estimate, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      Method,
      ": Cor = ",
      round(cor_val, 2)
    ),
    x = 21.58,
    y = c(24.35, 24.10, 23.85)
  )

cor_lab_H


p_scatter_H  <- make_scatter_plot(df_long_H,  cor_lab_H,  c(21.5, 24.5), c(21.5, 24.5))


#jpeg("~/Downloads/Figure10c.png",width=6,height=6, res = 300, units = "in")
print(p_scatter_H)
#dev.off
