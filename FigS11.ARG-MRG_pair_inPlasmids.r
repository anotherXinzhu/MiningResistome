setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")

library(dplyr)
library(data.table)

load("_data/otherARGtypes.RData")
load("_data/otherMetalRGtypes.RData")

leastDistanceDf.ArgMrg_plasmid <- 
  fread("_data/nonRegARG-MRG_pairs_inPlasmids.csv", data.table = F) %>%
  filter(grepl("plasmid", label)) %>%
  filter(leastDistance < 100000) %>% filter(same_ORF == "") 

dat <- leastDistanceDf.ArgMrg_plasmid

dat$core_label_1[which(dat$core_label_1 %in% otherARGtypes)] <- "other"
dat$neighbor_label_2[which(dat$neighbor_label_2 %in% otherMetals)] <- "other_metal_resistance"

# plot in drugs -------
dat_pairs <- dat %>%
  group_by(sampleCollection, core_label_1, neighbor_label_2) %>%
  summarise(n=n()) %>% as.data.frame() %>% 
  merge((dat %>% group_by(sampleCollection,core_label_1) %>% summarise(sum=n())),
        by=c("sampleCollection","core_label_1")) %>%
  mutate(freq_inDrug=n/sum)

MRGtype_rank <- c("Arsenic_resistance", "Chromium_resistance", "Cobalt_resistance","Copper_resistance",
                  "Iron_resistance","Lead_resistance","Mercury_resistance","Multi-metal_resistance",
                  "Nickel_resistance","Zinc_resistance","other_metal_resistance") 

dat_pairs$neighbor_label_2[which(dat_pairs$neighbor_label_2== "other_metals") ] <- "other_metal_resistance"
dat_pairs$sampleCollection <- factor(dat_pairs$sampleCollection,levels = c("pubAMD","Tailing","AMDSedi"))
dat_pairs$neighbor_label_2 <- factor(dat_pairs$neighbor_label_2,
                                     levels = MRGtype_rank)

MetalType_colorsDf2 <- cbind.data.frame(
  Colors=c("#EFC000FF","#4DBBD5FF","#00A087FF","steelblue3","#F39B7FFF","#8491B4FF","#91D1C2FF","#E64B35FF" ,
           "#7E6148FF","#B09C85FF","gray"),
  DrugTypes=c("Arsenic_resistance", "Chromium_resistance", "Cobalt_resistance","Copper_resistance",
              "Iron_resistance","Lead_resistance","Mercury_resistance","Multi-metal_resistance",
              "Nickel_resistance","Zinc_resistance","other_metal_resistance"),
  stringsAsFactors=F
)
Colors = sapply(MRGtype_rank, function(x)MetalType_colorsDf2$Colors[which(MetalType_colorsDf2$DrugTypes == x)])

x_rank_df_pubAMD <- c("sulfonamide","multidrug","MLS","bacitracin","other","phenicol","aminoglycoside","tetracycline",
                      "beta−lactam","glycopeptide","mupirocin")
x_rank_df_Tailing <- c("sulfonamide","other","glycopeptide","bacitracin","multidrug","phenicol","tetracycline",
                       "MLS","aminoglycoside","beta−lactam","mupirocin")
x_rank_df_AMDSedi <- c("sulfonamide","phenicol","MLS","other","tetracycline","multidrug","bacitracin","aminoglycoside",
                       "beta−lactam","glycopeptide","mupirocin")



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

P_AMDSedi

ggsave(ggarrange(P_pubAMD+theme(legend.position="none") ,
                 P_Tailing+theme(legend.position="none") ,
                 P_AMDSedi+theme(legend.position="none") ,
                 nrow = 1,ncol = 3),
       filename = "FigS10.MRGcompP_drugs_plasmid.pdf", width = 9,height = 4)
