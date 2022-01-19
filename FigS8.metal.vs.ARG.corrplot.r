library(data.table)
library(dplyr)


ARG_df.l <- fread("_data/nonRegARG.abund.csv")
totalARG_df <- ARG_df.l %>% group_by(sampleID ) %>% summarise(totalARG= sum(DepthPG))


#tailings ----------
phychem <- fread("_data/physico-chemical.tailings.csv")

tailings_df <- merge(totalARG_df, phychem, by.x="sampleID", by.y = "V1") %>% 
  select_if(is.numeric) %>% 
  select(-longitude, -latitude, -MAT, -MAP, -TP, -pH, -EC, -TN, -TC, -multiMetalIndex_total, -multiMetalIndex_avail)

cor_df  <- cor(tailings_df) %>% as.data.frame() %>% select(totalARG)

tailings_df[,-1] <- log(tailings_df[,-1])
cor_df_log <- cor(tailings_df) %>% as.data.frame() %>% select(totalARG)

library(ggplot2)
p_tailings <- 
  ggplot(data = tailings_df, aes(x=Zn, y=totalARG)) +
  geom_point(color="#CD534CFF", size=2) +
  geom_point(shape=21,size=2) +
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.grid = element_blank())

cor.test(tailings_df$totalARG, tailings_df$Zn)
# r=0.397,  p-value = 1.316e-05


p_tailings

# AMD sediment --------------
library(readxl)
excel_sheets("_data/physico_chemical.amdsediments.xlsx")
phychem <- read_excel("_data/physico_chemical.amdsediments.xlsx",sheet = "metals.log")

load("_data/1.GeoInformaton.RData")
phychem$sample <-
  sapply(phychem$sample, function(x) AMDsedi_geo_df$sampleID[which(AMDsedi_geo_df$sample == x)])
AMDsedi_df <- merge(totalARG_df, phychem, by.x = "sampleID", by.y = "sample")  %>%
  select_if(is.numeric) 
cor_df_log <- cor(AMDsedi_df) %>% as.data.frame() %>% select(totalARG)


library(ggplot2)
p_AMDsedi <- 
  ggplot(data = AMDsedi_df, aes(x=avail_Mn, y=totalARG)) +
  geom_point(color="#EFC000FF", size=2) +
  geom_point(shape=21,size=2) +
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.grid = element_blank())

cor.test(AMDsedi_df$totalARG, AMDsedi_df$avail_Mn)
# r=0.237,  p-value = 0.025

library(ggpubr)
ggsave(ggarrange(p_tailings, p_AMDsedi, nrow = 1, ncol = 2),
       device = "pdf", filename = "FigS6.scatterplot.pdf", width = 6, height = 2.5)
