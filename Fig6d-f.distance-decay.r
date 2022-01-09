load("_data/1.GeoInformaton.RData")


library(dplyr)
library(tibble)
library(geosphere)
library(vegan)

# equations  -------------------------------------
lm_rp_plain <- function(df, y_column, x_column){
  y = y_column
  x = x_column
  m <- lm(y ~ x, df);
  corrTest=cor.test(y, x);
  a = round(unname(coef(m)[1]), digits = 2)
  b = round(unname(coef(m)[2]), digits = 2)
  R = round(corrTest$estimate,digits = 3)
  pvalue = ifelse(corrTest$p.value == 0, 
                  "< 2.2e-16", 
                  paste("= ",formatC(corrTest$p.value, format = "e", digits = 3),sep = ""))
  
  return(paste("y = ",a," + ",b,"x, r = ",R,", p ", pvalue, sep = ""))
}



# public AMD data -------------------------------------

# data 
load("_data/nonRegARG.count.publicAMD.PcoaSourceData.RData")
all(groups_pubAMD$sampleID %in% publicAMD_info_df$sample_no)

geo.pubAMD_bySite <- groups_pubAMD %>% 
  mutate(longitude = sapply(sampleID, function(x)publicAMD_info_df$longitude[publicAMD_info_df$sample_no == x])) %>%
  mutate(latitude = sapply(sampleID, function(x)publicAMD_info_df$latitude[publicAMD_info_df$sample_no == x])) %>%
  group_by(location_Eng) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)


arg_df_bySite <- arg_df %>%
  rownames_to_column("sampleID") %>%
  mutate(site = sapply(sampleID, function(x)publicAMD_info_df$location_Eng[publicAMD_info_df$sample_no==x])) %>%
  group_by(site) %>%
  summarise_if(is.numeric, sum, na.rm=T) %>%
  column_to_rownames("site")

# distance
rownames(arg_df_bySite) == geo.pubAMD_bySite$location_Eng

geo_dist <- as.vector(distm(geo.pubAMD_bySite[,c('longitude','latitude')], geo.pubAMD_bySite[,c('longitude','latitude')], fun=distVincentyEllipsoid))

bray_dis <- as.matrix(vegdist(arg_df_bySite, method = 'bray', binary = T))
bray_simil <- as.vector(1- bray_dis)

plotDat <- cbind.data.frame(geo_dist, bray_simil)

# Plotting

library(ggplot2)

disDecay_pubAMD <- ggplot(data=plotDat) +
  geom_point(aes(x=geo_dist, y=bray_simil), shape=21, fill="#6A6599FF",size=2) +
  geom_smooth(aes(x=geo_dist, y=bray_simil), method = "lm") +
  theme_bw() + theme(panel.grid = element_blank())
disDecay_pubAMD
  
eqt_lm <- lm_rp_plain(df=plotDat, plotDat$geo_dist/1000, plotDat$bray_simil)
eqt_lm
"y = 11021 + -8947.58x, r = -0.315, p = 6.890e-06"

# Tailings -------------------------------------

# data 
ARG_df.l <- fread("_data/nonRegARG.abund.csv")

arg_df_bySite <- ARG_df.l %>% filter(grepl("\\.T", sampleID)) %>%
  dcast(sampleID ~ ARG, value.var = "DepthPG") %>% 
  mutate(site = sapply(sampleID, function(x) soil_geo_df$location_Eng[soil_geo_df$sampleID == x])) %>%
  group_by(site) %>%
  summarise_if(is.numeric, mean, na.rm=T) %>%
  column_to_rownames("site")
geo_bySite <- soil_geo_df %>% 
  filter(sampleID %in% ARG_df.l$sampleID) %>%
  group_by(location_Eng) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

# distance
rownames(arg_df_bySite) == geo_bySite$location_Eng
geo_dist <- as.vector(distm(geo_bySite[,c('longitude','latitude')], geo_bySite[,c('longitude','latitude')], fun=distVincentyEllipsoid))

bray_dis <- as.matrix(vegdist(arg_df_bySite, method = 'bray', binary = F))
bray_simil <- as.vector(1- bray_dis)

# Plotting
plotDat <- cbind.data.frame(geo_dist, bray_simil)

disDecay_tailings <- ggplot(data=plotDat) +
  geom_point(aes(x=geo_dist, y=bray_simil), shape=21, fill= "#CD534CFF",size=2) +
  geom_smooth(aes(x=geo_dist, y=bray_simil), method = "lm") +
  theme_bw() + theme(panel.grid = element_blank())
disDecay_tailings

eqt_lm <- lm_rp_plain(df=plotDat, plotDat$geo_dist/1000, plotDat$bray_simil)
eqt_lm
"y = 2374.95 + -1824.71x, r = -0.341, p = 1.268e-42"

# AMD sediments -------------------------------------

# data 
ARG_df.l <- fread("_data/nonRegARG.abund.csv")

arg_df_bySite <- ARG_df.l %>% filter(grepl("^LD", sampleID)) %>%
  dcast(sampleID ~ ARG, value.var = "DepthPG") %>% 
  mutate(site = sapply(sampleID, function(x) AMDsedi_geo_df$location_Eng[AMDsedi_geo_df$sampleID == x])) %>%
  group_by(site) %>%
  summarise_if(is.numeric, mean, na.rm=T) %>%
  column_to_rownames("site")
geo_bySite <- AMDsedi_geo_df %>% 
  filter(sampleID %in% ARG_df.l$sampleID) %>%
  group_by(location_Eng) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

# distance
rownames(arg_df_bySite) == geo_bySite$location_Eng
geo_dist <- as.vector(distm(geo_bySite[,c('longitudes','latitudes')], geo_bySite[,c('longitudes','latitudes')], fun=distVincentyEllipsoid))

bray_dis <- as.matrix(vegdist(arg_df_bySite, method = 'bray', binary = F))
bray_simil <- as.vector(1- bray_dis)

# Plotting
plotDat <- cbind.data.frame(geo_dist, bray_simil)

disDecay_AMDsedi <- ggplot(data=plotDat) +
  geom_point(aes(x=geo_dist, y=bray_simil), shape=21, fill= "#EFC000FF",size=2) +
  geom_smooth(aes(x=geo_dist, y=bray_simil), method = "lm") +
  theme_bw() + theme(panel.grid = element_blank())
disDecay_AMDsedi

eqt_lm <- lm_rp_plain(df=plotDat, plotDat$geo_dist/1000, plotDat$bray_simil)
eqt_lm
"y = 1512.74 + -1284.84x, r = -0.373, p = 1.114e-14"


 
# integrate three figures 

library(ggpubr)

ggsave(ggarrange(disDecay_pubAMD, disDecay_tailings, disDecay_AMDsedi, nrow = 1, ncol = 3),
       device = "pdf", width = 12,height = 3,filename = "Fig5d-f.distanceDecay.pdf")


