setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/22.VirulenceFactor")

load("4.organized_Vir_54Pathog.30nonPatho.Bins.RData")
load("3.organized_tox_53Pathog.30nonPatho.Bins.RData")



library(dplyr)

plotDat.VF <- 
  Vir.taxon_df %>% 
  filter(Virulence_confidence_level != "-") %>%
  group_by(bin, Virulence_confidence_level) %>%
  summarise(numORF = n(), pathogenicity = unique(pathogenicity), species = unique(corrected_species)) %>%
  mutate(bin.abb = sub("\\.filtered","",sub("\\.contigs\\.fa","", bin))) %>%
  mutate(Label = paste(bin.abb, "|  ", species, sep = "")) %>% 
  mutate(VF.type = Virulence_confidence_level) %>%
  select(-Virulence_confidence_level)


plotDat.T <- 
  Tox.taxon_df %>% 
  filter(Toxin_confidence_level != "-") %>%
  group_by(bin, Toxin_confidence_level) %>%
  summarise(numORF = n(), pathogenicity = unique(pathogenicity), species = unique(corrected_species)) %>%
  mutate(bin.abb = sub("\\.filtered","",sub("\\.contigs\\.fa","", bin))) %>%
  mutate(Label = paste(bin.abb, "|  ", species, sep = ""))  %>% 
  mutate(VF.type = Toxin_confidence_level) %>%
  select(-Toxin_confidence_level)

plotDat <- bind_rows(plotDat.VF %>% filter(grepl("^[12]", VF.type, perl = T)), plotDat.T)


orderDf <- plotDat %>% select(Label, species, pathogenicity) %>% unique() %>% 
  arrange(species) %>% arrange(pathogenicity)
plotDat$Label <- factor(plotDat$Label, levels = orderDf$Label)


# data for average value segment ------------------
avg_df <- plotDat %>% 
  group_by(pathogenicity, VF.type) %>%
  summarise(avg.numORF = mean(numORF), sd.numORF = sd(numORF)) 

test <- sapply(avg_df$pathogenicity, 
               function(x){
                 if(x =="") return(c(1,30)) else return(c(31,84))
               }) %>% t() %>% as.data.frame()

colnames(test)<- c("X1","X2")

avg_df <- cbind(avg_df, test)
avg_df <- as.data.frame(avg_df)
sapply(avg_df, class)

# plot -------------------
# plot individual ####
library(ggplot2)
p1<- ggplot() +
  geom_point(data =plotDat ,
             aes(x=Label, y=numORF, group=VF.type, color=VF.type), 
             size=2) +
  geom_segment(data = avg_df, aes(x = X1, y = avg.numORF, xend = X2, yend = avg.numORF, color = VF.type), 
               # linetype ="dashed",
               size=1.5) +
  # scale_color_manual(values = c( "#0073C2FF", "#EFC000FF", "#CD534CFF","#868686FF")) +
  scale_color_manual(values = c("red3","steelblue4","pink2","lightsteelblue3")) +
  # scale_x_discrete(labels = sapply(strsplit(levels(plotDat$Label),"|", fixed = T),"[[",2))
  # facet_grid(.~pathogenicity, scales = "free_x") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) 

p1


# plot average and significance ####
library(ggpubr)
p2<-ggplot(avg_df) +
  geom_bar(aes(x=pathogenicity, y=avg.numORF, fill=VF.type,), stat="identity", width = 0.7) +
  geom_pointrange( aes(x=pathogenicity, y=avg.numORF, ymin=avg.numORF-sd.numORF, ymax=avg.numORF+sd.numORF), 
                   colour="darkgray", alpha=1, size=0.7) +
  facet_grid(VF.type~., scales = "free") +
  scale_fill_manual(values = c("red3","steelblue4","pink2","lightsteelblue3")) +
  theme_bw()+ theme( panel.grid = element_blank()) 


p2
# calculate significance 
plotDat %>%
  group_by(VF.type ) %>%
  do(w = wilcox.test(numORF ~ pathogenicity, data=.)) %>%
  summarise(VF.type, wilcox.p = w$p.value)

ggsave(ggarrange(p1,  widths = c(9,1)),
       filename = "FigS.VF_n_toxinGene.pdf", device = "pdf", width = 15, height = 8)
ggsave(p2+theme(legend.position = "none"),
       filename = "FigS.VF_n_toxinGene_avg.pdf",device = "pdf", width = 1.8, height = 5.5)
