here::i_am("Code/Figure10a.R")

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
     
library(tidyverse)

df87 <- read.csv("~/Downloads/DA3/tempBv5.csv")
df92 <- data.frame(year = ghcn[aa,1], temperature = ghcnave[aa])
colnames(df92) <- c("year", "GHCN")

df94 <- inner_join(df87, df92, by = "year")


df_long <- df94 %>%
  select(year, GHCN, predicted, REACHES, LME) %>%
  pivot_longer(
    cols = c(predicted, REACHES, LME),
    names_to = "Method",
    values_to = "Estimate"
  ) %>%
  mutate(
    Method = recode(Method,
                    predicted = "Assimilated",
                    REACHES  = "REACHES",
                    LME      = "LME")
  )


cor_lab <- df_long %>%
  group_by(Method) %>%
  summarise(cor_val = cor(GHCN, Estimate, use = "complete.obs")) %>%
  mutate(
    label = paste0(Method, ": Cor = ", round(cor_val, 2)),
    x = 10.05,
    y = c(12.9, 12.65, 12.4)
  )

base_text_size   <- 23   
corr_text_size   <- 7    
legend_text_size <- 18   

make_scatter_plot <- function(df_long, cor_lab, xlim_range, ylim_range) {
  ggplot(df_long,
         aes(x = GHCN, y = Estimate, color = Method, shape = Method)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_text(
      data = cor_lab,
      aes(x = x, y = y, label = label, color = Method),
      inherit.aes = FALSE,
      show.legend = FALSE,
      hjust = 0,
      size = corr_text_size
    ) +
    scale_color_manual(
      values = c(
        "Assimilated" = "#F8766D",
        "LME" = "#00BA38",
        "REACHES" = "#619CFF"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Assimilated" = 16,
        "LME" = 17,
        "REACHES" = 15
      )
    ) +
    coord_equal(xlim = xlim_range, ylim = ylim_range) +
    labs(x = "GHCN", y = "Estimate") +
    theme_bw() +
    theme(
      text = element_text(size = base_text_size),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = legend_text_size)
    )
}

p_scatter_BJ <- make_scatter_plot(df_long, cor_lab, c(10, 13), c(10, 13))

#jpeg("~/Downloads/Figure10(a).png",width=6,height=6, res = 300, units = "in")
print(p_scatter_BJ)
#dev.off
