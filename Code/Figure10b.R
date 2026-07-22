library(ggplot2)
library(readxl)
library(tidyverse)

ghcnn<-read_excel("~/Downloads/GHCNv4.xlsx",col_type = c("skip","numeric","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip","numeric","skip","skip"))

#Shanghai
ghcn<-ghcnn[443:503, ]



ghcnave<-array(0,dim = nrow(ghcn))
for (i in 1:nrow(ghcn)) {
  ghcnave[i]<-rowMeans(ghcn[i,-1],na.rm = T)/100
}
na_per_row <- apply(ghcn, 1, function(x) sum(is.na(x)))
aa<-which(na_per_row < 1)
    
    df87 <- read.csv("~/Downloads/DA3/tempSv5.csv")
    df92 <- data.frame(year=ghcn[aa,1], temperature=ghcnave[aa])
    colnames(df92)<- c("year", "temperature")
    df94<-inner_join(df87,df92, by="year")

df94_S <- df94 %>%
  rename(GHCN = temperature)


df_long_S <- df94_S %>%
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

cor_lab_S <- df_long_S %>%
  group_by(Method) %>%
  summarise(
    cor_val = cor(GHCN, Estimate, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(Method, ": Cor = ", round(cor_val, 2)),
    x = 15.05,
    y = c(17.88, 17.63, 17.38)
  )

p_scatter_S  <- make_scatter_plot(df_long_S,  cor_lab_S,  c(15, 18), c(15, 18))


#jpeg("~/Downloads/Figure10(b).png",width=6,height=6, res = 300, units = "in")
print(p_scatter_S)
#dev.off
