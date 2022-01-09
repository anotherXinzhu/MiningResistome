library(data.table)
library(dplyr)
library(reshape2)
library(grid) #function grobTree
library(tibble)
library(ggsci)

# Figure 5a. pcoa of global public AMD samples ---------------------------------------
load("_data/nonRegARG.count.publicAMD.PcoaSourceData.RData")
# arg_df are the counts of ARG orfs in each sample
# arg_df_subset are the ARG counts incuding only ARGs detected in European samples where reads are not available


library(ggplot2)
library(vegan)
dat_arg <- arg_df
groups_df <- groups_pubAMD


arg.dist <- vegdist(dat_arg, method = "bray", binary = T)

pcoa <- cmdscale(arg.dist, k = (nrow(dat_arg) - 1), eig = TRUE)
pcoa_eig <- pcoa$eig
pcoa_exp <- pcoa$eig/sum(pcoa$eig) 


#ggplot2 作图

#添加分组信息
site <- data.frame(pcoa$point)[1:2] 
site$name <- rownames(site)
groups_df$sampleID == site$name #确定groups的排列与数据排序一致


#前 2 轴解释量
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')

site$group <- groups_df$position #color by position
arg.ano <- with(dat_arg, anosim(arg.dist, groups_df[,"position"]))
#summary(arg.ano)
arg.ano$statistic # R
arg.ano$signif # p

#ggplot2 作图
grob1 <- grobTree(textGrob(paste("Anosim.R=",round(arg.ano$statistic, 4),
                                 ", pvalue=",arg.ano$signif,
                                 sep = ""), 
                           x=0.05,  y=0.1, hjust=0, gp=gpar(col="black", fontsize=13)))
numColors <- length(unique(groups_df$position))
Pcoa_pubAMD <- ggplot(data = site, aes(X1, X2)) +
  geom_point(aes(color = group),size=2.8) +
  geom_point(size=2.8,shape=21) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.2, show.legend = TRUE) +    #添加置信椭圆，注意不是聚类
  scale_color_manual(values = pal_d3("category20")(20)[1:numColors]) +
  scale_fill_manual(values = pal_d3("category20")(20)[1:numColors]) +
  theme(#panel.grid.major = element_line(color = 'gray', size = 0.2), 
        panel.grid.major = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        plot.title = element_text(hjust = 0.5), legend.position = 'right') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  #geom_text(data = ARGs_top, aes(label = name), color = 'black', size = 3) +    #添 top20 丰度物种标签
  # labs(x = pcoa1, y = pcoa2, title = paste('PCoA (with top ',N,' abundant ARG),\n',distMeth,", Binary: ",ifBinary,sep = "")) #+
  labs(x = pcoa1, y = pcoa2,sep = "")+
  annotation_custom(grob1)

Pcoa_pubAMD



# Fig5b Tailings ----------------------------------------
ARG_df.l <- fread("_data/nonRegARG.abund.csv") %>% dcast( sampleID ~ ARG, value.var = "DepthPG") 
load("_data/1.GeoInformaton.RData")

arg_df_Tailing <- ARG_df.l %>% filter(!grepl("^LD",sampleID)) %>% filter(!(grepl("pubAMD",sampleID))) %>% column_to_rownames("sampleID")

soil_geo_df$position <- sapply(soil_geo_df$location_Eng,
                               function(x) all_geoBySite_df$position[which(all_geoBySite_df$location_Eng==x)])
groups_Tailing <- soil_geo_df[match(rownames(arg_df_Tailing),soil_geo_df$sampleID),] %>% select(sampleID, location_Eng,position)
groups_Tailing$sampleID == rownames(arg_df_Tailing)



dat_arg <- arg_df_Tailing
groups_df <- groups_Tailing %>% mutate(NS = sub("([NS])[EWM]", "\\1", position))


arg.dist <- vegdist(dat_arg, method = "bray", binary = TRUE)

pcoa <- cmdscale(arg.dist, k = (nrow(dat_arg) - 1), eig = TRUE)
pcoa_eig <- pcoa$eig
pcoa_exp <- pcoa$eig/sum(pcoa$eig) 


#ggplot2 作图

#添加分组信息
site <- data.frame(pcoa$point)[1:2] 
site$name <- rownames(site)
groups_df$sampleID == site$name #确定groups的排列与数据排序一致


#前 2 轴解释量
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')


g="position"
site$group <- groups_df[,g] #point color by position

#为了画椭圆的：
site$circle_grp <- groups_df[,"NS"] # circle color by NS

arg.ano <- with(dat_arg, anosim(arg.dist, groups_df[,g]))
grob1 <- grobTree(textGrob(paste("Anosim.R=",round(arg.ano$statistic, 4),
                                 ", pvalue=",arg.ano$signif,
                                 sep = ""), 
                           x=0.05,  y=0.1, hjust=0, gp=gpar(col="black", fontsize=13)))

Pcoa_Tailings <- ggplot(data = site, aes(X1, X2)) +
  geom_point(aes(color = group),size=2.8) +
  geom_point(size=2.8,shape=21) +
  stat_ellipse(aes(fill = circle_grp), geom = 'polygon', level = 0.95, alpha = 0.2, show.legend = TRUE) +    #添加置信椭圆，注意不是聚类
  scale_color_manual(values = c('lightskyblue2', "steelblue3",'orange1',"red1",'orangered3')) +
  scale_fill_manual(values = c('blue', 'red1')) +
  theme(#panel.grid.major = element_line(color = 'gray', size = 0.2), 
        panel.grid.major = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        plot.title = element_text(hjust = 0.5), legend.position = 'right') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  #geom_text(data = ARGs_top, aes(label = name), color = 'black', size = 3) +    #添 top20 丰度物种标签
  # labs(x = pcoa1, y = pcoa2, title = paste('PCoA (with top ',N,' abundant ARG),\n',distMeth,", Binary: ",ifBinary,sep = "")) #+
  labs(x = pcoa1, y = pcoa2,sep = "") +
  #geom_text(aes(label = name), color = 'black', size = 3) +  #看下几个特殊点是什么样品
  annotation_custom(grob1)

Pcoa_Tailings

# Fig5c. AMDsediments -------------------------------------
arg_df_AMDSedi <- ARG_df.l %>% filter(grepl("^LD",sampleID)) %>% column_to_rownames("sampleID")
# create grouping table ####
AMDsedi_geo_df$position <- sapply(AMDsedi_geo_df$location_Eng,
                                  function(x) all_geoBySite_df$position[which(all_geoBySite_df$location_Eng==x)])

groups_AMDSedi <- AMDsedi_geo_df[match(rownames(arg_df_AMDSedi),AMDsedi_geo_df$sampleID),] %>% select(sampleID, location_Eng, sample,position,longitudes,latitudes)
groups_AMDSedi$sampleID == rownames(arg_df_AMDSedi)


dat_arg <- arg_df_AMDSedi
groups_df <- groups_AMDSedi 
rownames(dat_arg) == groups_df$sampleID #确认顺序一致


arg.dist <- vegdist(dat_arg, method = "bray", binary = F)
pcoa <- cmdscale(arg.dist, k = (nrow(dat_arg) - 1), eig = TRUE)
pcoa_eig <- pcoa$eig
pcoa_exp <- pcoa$eig/sum(pcoa$eig) 

#添加分组信息
site <- data.frame(pcoa$point)[1:2] 
site$name <- rownames(site)
#sampleInfo_df <-  sampleInfo_df[match(site$name, sampleInfo_df$sampleID),] #为了排列颜色，画图里面用
groups_df$sampleID == site$name #确定groups的排列与数据排序一致

#前 2 轴解释量
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')


site$group <- groups_df$position #color by position
arg.ano <- with(dat_arg, anosim(arg.dist, groups_df[,"position"]))
#summary(arg.ano)
arg.ano$statistic # R
arg.ano$signif # p

#ggplot2 作图
grob1 <- grobTree(textGrob(paste("Anosim.R=",round(arg.ano$statistic, 4),
                                 ", pvalue=",arg.ano$signif,
                                 sep = ""), 
                           x=0.05,  y=0.1, hjust=0, gp=gpar(col="black", fontsize=13)))

Pcoa_AMDsedi <- ggplot(data = site, aes(X1, X2)) +
  geom_point(aes(color = group),size=2.8) +
  geom_point(size=2.8,shape=21) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.2, show.legend = TRUE) +    #添加置信椭圆，注意不是聚类
  scale_color_manual(values = c('gold2', 'orangered3', 'mediumseagreen',"steelblue3")) +
  scale_fill_manual(values = c('gold2', 'orangered3', 'mediumseagreen',"steelblue3")) +
  theme(#panel.grid.major = element_line(color = 'gray', size = 0.2), 
        panel.grid.major = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        plot.title = element_text(hjust = 0.5), legend.position = 'right') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  labs(x = pcoa1, y = pcoa2) +
  annotation_custom(grob1) 

Pcoa_AMDsedi

library(ggpubr)
ggsave(ggarrange(Pcoa_pubAMD, Pcoa_Tailings, Pcoa_AMDsedi, 
                 nrow = 1,ncol = 3),
       device = "pdf", width = 12,height = 2.7,filename = "Fig5a-c.integratedPCoA.pdf")



