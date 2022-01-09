library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer) #自己根据所需要的颜色数量来制备调色板


leastDistanceDf_ARG_MetalRG <- fread("_data/nonRegARG-MRG_coOccurrence.txt", )
load("_data/otherARGtypes.RData")
load("_data/otherMetalRGtypes.RData")


dat <- leastDistanceDf_ARG_MetalRG %>% filter(leastDistance < 100000)
dat$core_label_1[which(dat$core_label_1 %in% otherARGtypes)] <- "other"
dat$neighbor_label_2[which(dat$neighbor_label_2 %in% otherMetals)] <- "other_metals"


x_rank_df_pubAMD <- c("sulfonamide","multidrug","MLS","bacitracin","other","phenicol","aminoglycoside","tetracycline",
                      "beta−lactam","glycopeptide","mupirocin")
x_rank_df_Tailing <- c("sulfonamide","other","glycopeptide","bacitracin","multidrug","phenicol","tetracycline",
                       "MLS","aminoglycoside","beta−lactam","mupirocin")
x_rank_df_AMDSedi <- c("sulfonamide","phenicol","MLS","other","tetracycline","multidrug","bacitracin","aminoglycoside",
                       "beta−lactam","glycopeptide","mupirocin")


dat_pairs <- dat %>% 
  group_by(sampleCollection, core_label_1, neighbor_label_2) %>%
  summarise(n=n()) %>% as.data.frame() %>% 
  merge((dat %>% group_by(sampleCollection,core_label_1) %>% summarise(sum=n())),by=c("sampleCollection","core_label_1")) %>%
  mutate(freq_inDrug=n/sum)

MRGtype_rank <- c("Arsenic_resistance", "Chromium_resistance", "Cobalt_resistance","Copper_resistance",
                  "Iron_resistance","Lead_resistance","Mercury_resistance","Multi-metal_resistance",
                  "Nickel_resistance","Zinc_resistance","other_metal_resistance") #从_colors.RData拷过来的

dat_pairs$neighbor_label_2[which(dat_pairs$neighbor_label_2== "other_metals") ] <- "other_metal_resistance"
dat_pairs$sampleCollection <- factor(dat_pairs$sampleCollection,levels = c("pubAMD","Tailing","AMDSedi"))
dat_pairs$neighbor_label_2 <- factor(dat_pairs$neighbor_label_2,
                                     levels = MRGtype_rank)

# colors 
library(scales)
library(ggsci)
tmp = pal_npg("nrc")(10); tmp[8] = "#EFC000FF"; tmp[4] = "steelblue3"
tmp = replace(tmp, c(1, 8), tmp[c(8,1)]);
MetalType_colorsDf2 <- cbind.data.frame(
  Colors=c(tmp, "gray"),
  DrugTypes=c("Arsenic_resistance", "Chromium_resistance", "Cobalt_resistance","Copper_resistance",
              "Iron_resistance","Lead_resistance","Mercury_resistance","Multi-metal_resistance",
              "Nickel_resistance","Zinc_resistance","other_metal_resistance"),
  stringsAsFactors=F
) 

Colors = sapply(MRGtype_rank, function(x)MetalType_colorsDf2$Colors[which(MetalType_colorsDf2$DrugTypes == x)])


for(dc in levels(dat_pairs$sampleCollection)){
  #  dc = levels(dat_pairs$sampleCollection)[1]
  x_rank_df <- eval(parse(text = paste("x_rank_df_", dc, sep = "")))
  
  dat_pairs_tmp <- dat_pairs %>% filter(sampleCollection == dc) %>% filter(core_label_1 != "unclassified")
  dat_pairs_tmp$core_label_1 <- factor(dat_pairs_tmp$core_label_1,
                                       levels = x_rank_df)
  P <- ggplot(dat_pairs_tmp,aes(x=core_label_1, y=freq_inDrug, fill=neighbor_label_2)) +
    geom_col() + 
    scale_fill_manual(values = Colors) +
    theme_bw() +
    theme (axis.text.x = element_text (angle = 90, hjust = 1 ) ) + xlab("") +ylab("Nearest MRG composition") +
    #guides(fill = guide_legend(nrow = 1)) +   
    #theme(legend.position="none") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  P
  assign(paste("P_",dc,sep = ""),P, envir = .GlobalEnv)
}

ggsave(ggarrange(P_pubAMD+theme(legend.position="none") ,
                 P_Tailing+theme(legend.position="none") ,
                 P_AMDSedi+theme(legend.position="none") ,
                 nrow = 1,ncol = 3),
       filename = "Fig4d-f.MRGcompP_drugs_ncol3.pdf", width = 9,height = 4)



# plot in overall -------


dat_pairs_overall <- dat %>% group_by(sampleCollection, neighbor_label_2) %>% summarise(n=n()) %>% as.data.frame() %>% 
  merge((dat %>% group_by(sampleCollection) %>% summarise(sum=n())),by=c("sampleCollection")) %>%
  mutate(freq=n/sum)

dat_pairs_overall$sampleCollection <- factor(dat_pairs_overall$sampleCollection,levels = c("pubAMD","AMDSedi","Tailing"))

dat_pairs_overall$neighbor_label_2[which(dat_pairs_overall$neighbor_label_2== "other_metals") ] <- "other_metal_resistance"
dat_pairs_overall$neighbor_label_2 <- factor(dat_pairs_overall$neighbor_label_2,
                                             levels = MRGtype_rank)


MRGcompP_overall <- ggplot(dat_pairs_overall,aes(x=1, y=freq,fill =neighbor_label_2)) +
  geom_col() + 
  scale_fill_manual(values = Colors) +
  #guides(fill = guide_legend(nrow = 1)) +   
  theme(legend.position="none") +
  facet_wrap(vars(sampleCollection),ncol = 3,nrow = 1) +
  theme ( axis.text.x = element_blank (), axis.text.y = element_blank() ) + xlab("") +ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

MRGcompP_overall

ggsave(MRGcompP_overall+theme(legend.position="none"),
       filename = "Fig4d-f.Arg-MrgComposition-overall.pdf",device = "pdf", width = 3,height = 3.5)
