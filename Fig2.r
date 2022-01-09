setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")
library(dplyr)
library(ggplot2)
library(xlsx)
library(data.table)


ARG_df.l <- fread("_data/nonRegARG_sewage.freshwaterSedim.mine.txt")


# total ARG abundance -------------------------------------

total.ARGDepth_df <- 
  ARG_df.l %>% 
  group_by(sampleID) %>% 
  summarise(totalARG.DepthPG = sum(DepthPG),
            numARG = n(),
            sampleType = unique(sampleType)) 

total.ARGDepth_df$sampleType <- factor(total.ARGDepth_df$sampleType, 
                                       levels = c("freshwater sediments","mine","sewage")) 


library(ggpubr)

spType_colorDf2 <- cbind.data.frame(spType=c("freshwater sediments","mine","sewage"),
                                    Color=c("#2CA02CFF", "#1F77B4FF", "#FF7F0EFF"),
                                    stringsAsFactors=F)
Colors <- sapply(levels(total.ARGDepth_df$sampleType),
                 function(x){spType_colorDf2$Color[spType_colorDf2$spType == x]}) 

  
my_comparisons <- list(c("freshwater sediments","mine"),c("mine","sewage"))
p1<- ggboxplot(data = total.ARGDepth_df, 
               x="sampleType", y="totalARG.DepthPG", fill = "sampleType", alpha=0.5, width = 0.6, size = 1, 
               color="lightslategray",outlier.shape = NA) +
  geom_jitter(aes(x=sampleType, y=totalARG.DepthPG, color=sampleType), width = 0.15, size=2, alpha=0.8) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  scale_fill_manual(values = Colors) +
  scale_color_manual(values = Colors) +
  ylab("Total ARG (coverage /Gb)") + xlab("") +
  theme_bw() + scale_x_discrete(labels=c("Freshwater\nsediments","Mine","Sewage")) +
  theme(legend.position = "none", panel.grid = element_blank())

# mean(total.ARGDepth_df$totalARG.DepthPG[total.ARGDepth_df$sampleType == "mine"])
# mean(total.ARGDepth_df$totalARG.DepthPG[total.ARGDepth_df$sampleType == "freshwater sediments"])

p1


# ARG number ------------------------------------
p2 <- ggboxplot(data = total.ARGDepth_df, 
          x="sampleType", y="numARG", fill = "sampleType", alpha=0.5, width = 0.6, size = 1,
          color="lightslategray", outlier.shape = NA) +
  geom_jitter(aes(x=sampleType, y=numARG, color=sampleType), width = 0.15, size=2, alpha=0.8) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  scale_fill_manual(values = Colors) +
  scale_color_manual(values = Colors) +
  ylab("Number of ARG") + xlab("") +
  theme_bw() + scale_x_discrete(labels=c("Freshwater\nsediments","Mine","Sewage")) +
  theme(legend.position = "none", panel.grid = element_blank())

p2

# ARG composition ------------------------------
# calculate other ARG typs
tmp <- 
  ARG_df.l %>% 
  group_by(sampleID, drug_type) %>%
  summarise(DepthPG = sum(DepthPG), sampleType = unique(sampleType))



relAbund.df <- tmp %>% 
  group_by(sampleType,drug_type ) %>% summarise(avgDepth = mean(DepthPG)) %>% 
  mutate(freq = avgDepth/sum(avgDepth)) %>%
  reshape2::dcast(drug_type~sampleType) %>% 
  arrange(desc(mine))


otherTypes <- relAbund.df$drug_type[11:nrow(relAbund.df) ]
tmp$drug_type[tmp$drug_type %in% otherTypes] <- "other"


# plot by sample with significance 
plotDat_3bars <- tmp %>% filter(drug_type != "unclassified") %>%
  group_by(sampleType,drug_type ) %>% summarise(avgDepth = mean(DepthPG)) %>% 
  mutate(freq = avgDepth/sum(avgDepth)) %>%
  as.data.frame()

plotDat_3bars %>% group_by(sampleType) %>% summarise(s = sum(freq))


library(ggsci)
drugType_colorsDf <- cbind.data.frame(
  Colors=c(pal_npg("nrc")(10), "gray","gold3","orchid"),
  DrugTypes=c("aminoglycoside", "bacitracin", "beta-lactam","glycopeptide",
              "MLS","tetracycline","phenicol","multidrug",
              "mupirocin","sulfonamide","other","fosmidomycin","rifamycin"), # 后来添加了fosmidomycin 和 rifamycin是SARG
  stringsAsFactors=F
)
Colors = sapply(unique(plotDat_3bars$drug_type)[order(unique(plotDat_3bars$drug_type) )],
                function(x)drugType_colorsDf$Colors[drugType_colorsDf$DrugTypes == x] )
drugType_levels <- c(unique(plotDat_3bars$drug_type)[unique(plotDat_3bars$drug_type)!="other"],"other")
plotDat_3bars$drug_type <- factor(plotDat_3bars$drug_type, levels=rev(drugType_levels))

p3<-ggplot(plotDat_3bars) +
  geom_col(aes(x=sampleType ,y =freq, fill=drug_type), width = 0.6,color="white") +
  scale_fill_manual(values = Colors)+
  theme_bw() + theme(#axis.text.x = element_text(angle = 90),
    panel.grid = element_blank()) 
p3

# calculate significance of means comparison of different ARG types ---------------------
ARGtype_df.l <- tmp %>% filter(drug_type != "unclassified") %>%
  group_by(sampleID) %>% mutate(freq = DepthPG/sum(DepthPG)) %>%
  arrange(sampleType)  %>%
  data.frame()

# calculate significance between freshwater sediments and mine ####
ARGtype_df.l %>%
  filter(sampleType != "sewage") %>%
  group_by(drug_type) %>%
  do(tt = t.test(freq~sampleType, data=.)) %>%
  summarise(drug_type, ttest.p = tt$p.value) %>%
  mutate(p.sig = cut(ttest.p, breaks=c(-Inf, 0.001, 0.01, 0.05,Inf), labels=c("***","**","*",""))) 

# calculate significance between sewage and mine ####
ARGtype_df.l %>%
  filter(sampleType != "freshwater sediments") %>%
  group_by(drug_type) %>%
  do(tt = t.test(freq~sampleType, data=.)) %>%
  summarise(drug_type, ttest.p = tt$p.value) %>%
  mutate(p.sig = cut(ttest.p, breaks=c(-Inf, 0.001, 0.01, 0.05,Inf), labels=c("***","**","*",""))) 

#ggsave(ggarrange(p1,p2,p3, widths = c(3,3,4), nrow = 1), device = "pdf",
#       filename = "Fig2.mine.vs.sewage.vs.sediment_nonRegGenes.pdf", width = 11, height = 4.5)

ggsave(ggarrange(p1,p2, nrow = 1), device = "pdf",
       filename = "Fig2.a-b.pdf", width = 5, height = 3)


p3 + coord_flip() 
ggsave( p3 + coord_flip() ,
       device = "pdf",
       filename = "Fig2.c.pdf", width = 8, height = 2.8)
