setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")

library(dplyr)
library(ggplot2)

# Figure S?(a-b). plot histogram for ARG carrying contig lengths in bins ------------
load("_data/ARG-carryingContigLengths_bins.RData")


# plot contig lengths for bins recovered from Tailings  #####                                                                                                                          

p1 <- ggplot(data=Tailings.contigLengths, aes(x=lengths.kb)) + 
  geom_histogram(breaks=seq(2, 1243, by=1), 
                 col="darkgray", 
                 fill="#CD534CFF", 
                 alpha = .2) + 
  #xlim(c(0,520)) + 
  scale_x_log10() +
  #ylim(c(0,30)) +
  labs(title="Histogram for ARG-carrying contig lengths", x="Contig length (kb)", y="Count") +
  theme_bw()+theme(panel.grid = element_blank())

# calculate percentages of each length range
sum(Tailings.contigLengths$lengths.kb < 10)/nrow(Tailings.contigLengths) # 9.56%
sum(Tailings.contigLengths$lengths.kb >= 10 & Tailings.contigLengths$lengths.kb < 100)/nrow(Tailings.contigLengths) # 65.69%
sum(Tailings.contigLengths$lengths.kb >= 100 )/nrow(Tailings.contigLengths) # 24.75%


# plot contig lengths for bins recovered from AMD sediment   ##### 

library(ggplot2)
p2 <- ggplot(data=AMD.contigLengths, aes(x=lengths.kb)) + 
  geom_histogram(breaks=seq(2, 1578, by=1), 
                 col="darkgray", 
                 fill="#EFC000FF", 
                 alpha = .2) + 
  #xlim(c(0,520)) + 
  scale_x_log10() +
  #ylim(c(0,30)) +
  labs(title="Histogram for ARG-carrying contig lengths", x="Contig length (kb)", y="Count") +
  theme_bw()+theme(panel.grid = element_blank())

# calculate percentages of each length range
sum(AMD.contigLengths$lengths.kb < 10)/nrow(AMD.contigLengths) # 8.03%
sum(AMD.contigLengths$lengths.kb >= 10 & AMD.contigLengths$lengths.kb < 100)/nrow(AMD.contigLengths) # 61.43%
sum(AMD.contigLengths$lengths.kb >= 100 )/nrow(AMD.contigLengths) # 30.54%



# Figure S?(c). within bin eveness MAGs vs pseudo MAGs ----------------
library(xlsx)

# 20 MAGs mapping results --------------
MAG_covDf <- read.xlsx("_data/20Bin_mapping.xlsx",sheetIndex = 1)

MAGs_df <- MAG_covDf %>% 
  group_by(bin) %>% 
  summarise(n =  n(),
            sd =  sd(as.numeric(coverage)),
            avg = mean(as.numeric(coverage))) %>%
  mutate(sd.perc = sd/avg,
         se = sd/sqrt(n)) %>%
  arrange(desc(sd))


# 20 artificial MAGs mapping results --------------
PseudoBin.cov_df.all <- read.xlsx("_data/20pseudoBin_mapping.xlsx",sheetIndex = 1)

pseudo.MAGs_df <- PseudoBin.cov_df.all %>% 
  group_by(bin) %>% 
  summarise(n =  n(),
            sd =  sd(coverage),
            avg = mean(coverage)) %>%
  mutate(sd.perc = sd/avg,
         se = sd/sqrt(n)) %>%
  arrange(desc(sd)) 


# plotting --------------
library(ggplot2)
plotDat <- bind_rows(cbind.data.frame(MAGs_df, bin.type="MAG"),
                     cbind.data.frame(pseudo.MAGs_df, bin.type = "Pseudo.MAG"))

p3<-ggplot(plotDat) +
  geom_boxplot(aes(x=bin.type, y=sd.perc, fill=bin.type), outlier.shape = NA, alpha=0.2) +
  geom_jitter(aes(x=bin.type, y=sd.perc, color=bin.type), width = 0.2, size=2 ) +
  ylab("sd of coverage/avg of coverage") +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.title.x = element_blank(),
                     legend.position = "none")
p3

# p value 
wilcox.test(sd.perc ~ bin.type, data = plotDat)

library(ggpubr)

ggsave(ggarrange(p1+geom_vline(xintercept = 10, linetype="dashed") + geom_vline(xintercept = 100, linetype="dashed"), 
                 p2+geom_vline(xintercept = 10, linetype="dashed") + geom_vline(xintercept = 100, linetype="dashed"),
                 p3, widths = c(4,4,2), ncol = 3), 
       filename = 'FigS.binReliability.pdf', 
       device = "pdf", width = 10, height = 3)
