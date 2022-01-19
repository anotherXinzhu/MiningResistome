rm(list = ls())

setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")
load("_data/plasmid.sourceData.nonRegARGs.RData")
library(dplyr)
library(ggplot2)
library(ggsci)
library(scales)


# combine data into one df -------------------------
combined_dat <- dplyr::bind_rows(cbind.data.frame(AMDSedi_df, DataCollection="AMDSedi",stringsAsFactors=F),
                                 cbind.data.frame(Tailing_df, DataCollection="National",stringsAsFactors=F),
                                 cbind.data.frame(pubAMD_df, DataCollection="pubAMD",stringsAsFactors=F))
combined_dat <- combined_dat %>% 
  filter(plasmid != "unclassified")

# calculate number of ARG/kb on plasmid, chrmosome and unclassified contigs ------
PlaChr_sumLength_df <- combined_dat %>% dplyr::select(sample_name,contig_name,contig_length,plasmid) %>% unique() %>%
  dplyr::group_by(sample_name,plasmid) %>% 
  dplyr::summarise(sumLength = sum(contig_length)) %>% as.data.frame()

numARG_summary_df <- combined_dat %>% dplyr::group_by(sample_name,plasmid) %>% 
  dplyr::summarise(numARG=n(),DataCollection=unique(DataCollection)) %>% as.data.frame() %>%
  base::merge(PlaChr_sumLength_df, by=c("sample_name","plasmid")) %>%
  dplyr::mutate(numARG_perKB = numARG/(sumLength/1000)) 

numARG_summary_df$DataCollection <- factor(numARG_summary_df$DataCollection, 
                                           levels = c( "pubAMD", "National", "AMDSedi"))
p <- ggplot(numARG_summary_df,aes(x=DataCollection,y=numARG_perKB))+
  geom_boxplot(aes(color = plasmid, fill = plasmid), width = 0.6, size = 0.4, position = position_dodge(0.8),outlier.shape = NA) +
  scale_color_jco() +   scale_fill_jco(alpha = 0.7) + coord_cartesian(ylim=c(0, 0.6)) + 
  geom_jitter(aes( color = plasmid), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8))+
  xlab("") + theme_bw() + theme(panel.grid = element_blank())

p
ggsave(p, filename = "FigS11.numARGperKB_plasmid.pdf", width = 6, height = 4, device = "pdf")

#  calculate significance
numARG_summary_df %>%
  group_by(DataCollection) %>%
  do(w = wilcox.test(numARG_perKB ~ plasmid, data = .)) %>%
  summarise(DataCollection, wilcox.p = w$p.value)


numARG_summary_df %>%
  group_by(plasmid) %>%
  summarise(numARG = mean(numARG_perKB))
