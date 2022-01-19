library(data.table)
library(dplyr)

rm(list = ls())
load("_data/SPAdes_vs_MEGAhit.RData")
# plot contig stats -----------------------------
library(ggpubr)
contigStats <- bind_rows(MEGAhit_contigStats, SPAdes_contigStats)


FigS2a <- ggpaired(contigStats %>% filter(variable %in% c("numContig.1000","avgLen.contig1000","N50")),
                   x = "assemMethod", y = "value",
                   color = "assemMethod", line.color = "gray", line.size = 0.4,
                   palette = "jco")+
  stat_compare_means(paired = TRUE) +
  facet_wrap( vars(variable), scales = "free") +
  theme_bw() + theme(panel.grid = element_blank())
FigS2a

# plot ARG overlapped proportion ------------------------------------
# randomly select combinations of n samples 

Results <- NULL
for(n in c(1:length(sps)) ){
  # n is the number of mixed samples in a SPAdes-MEGAhit comparison\
  # n=2
  
  n_combns <- min(30, choose(length(sps),n))
  i_selected <- sample(1:choose(length(sps),n), n_combns)
  
  for(i_combn in i_selected){
    #  i_combn = i_selected[1]
    
    n_sps <- combn(sps, n)[,i_combn]
    
    SPAdes_ARGsubtypes <- unique(spades_argDf$ARG[which(spades_argDf$sample %in% n_sps)])
    MEGAhit_ARGsubtypes <- unique(megahit_argDf$ARG[which(megahit_argDf$sample %in% n_sps)])
    SPAdes_ARGtypes <- unique(spades_argDf$ARGtype[which(spades_argDf$sample %in% n_sps)])
    MEGAhit_ARGtypes <- unique(megahit_argDf$ARGtype[which(megahit_argDf$sample %in% n_sps)])
    
    # calculate total ARG subtypes and types detected 
    totalARGs <- unique(c(SPAdes_ARGsubtypes, MEGAhit_ARGsubtypes))
    totalARGtypes <- unique(c(SPAdes_ARGtypes, MEGAhit_ARGtypes))
    
    # calculate common ARG subtypes and types detected
    commonARGs <- intersect(SPAdes_ARGsubtypes, MEGAhit_ARGsubtypes)
    commonARGtypes <- intersect(SPAdes_ARGtypes, MEGAhit_ARGtypes)
    
    res_c  <- c("num_samples" = n, 
                "num_totalARGs" = length(totalARGs),
                "num_totalARGtypes" = length(totalARGtypes),
                "num_overlapped.ARGs" = length(commonARGs), 
                "num_overlapped.ARGtypes" = length(commonARGtypes))
    
    Results <- bind_rows(Results, res_c)
  }# loop through random combinations
}# loop through n, the number of mixed samples in a SPAdes-MEGAhit comparison


Results <- Results %>% mutate(perc.overlappedARG = num_overlapped.ARGs/num_totalARGs,
                              perc.overlappedARGtype = num_overlapped.ARGtypes/num_totalARGtypes)


# plot overlapped ARG subtypes  -------
library(ggplot2)
library(ggpubr)

FigS2b <- ggboxplot(data = Results %>% filter(num_samples %in% c(1,2,4,6,8)), 
                    x = "num_samples", y="perc.overlappedARG", fill = "num_samples",
                    outlier.shape = NA, add = "jitter" ) +
  theme_bw() +
  scale_fill_brewer(palette = "GnBu")

FigS2b

# plot Venn for overlapped ARG subtypes by 10 samples ------------

megahit_ARGsubtypes <- unique(megahit_argDf$ARG)
spades_ARGsubtypes <- unique(spades_argDf$ARG)

plotDat <- merge(cbind.data.frame(megahit_ARGsubtypes, megahit = T, stringsAsFactors=F),
                 cbind.data.frame(spades_ARGsubtypes, spades = T, stringsAsFactors=F), 
                 by.x = "megahit_ARGsubtypes", by.y = "spades_ARGsubtypes", all = T)
plotDat[is.na(plotDat)] <- F

library(eulerr)
?euler
fit <- euler(plotDat[,-1])
FigS2c <- plot(fit, col="gray", fills = c("#0073C2FF", "#EFC000FF"), alpha = 0.6, quantities = TRUE)



pdf("FigS2.MEGAhit_vs_SPAdes.pdf",width = 10,height = 7)
ggarrange(FigS2a + theme(legend.position = "none"), 
          ggarrange(FigS2b + theme(legend.position = "none"), FigS2c, nrow = 1, widths=c(0.34, 0.66)),
          ncol = 1)
dev.off()
