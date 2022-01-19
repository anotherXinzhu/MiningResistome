rm(list = ls())


library(data.table)
plotDat <- fread("_data/pathogen_vs_nonpathogen.txt", data.table = F)

plotDat$variable <- factor(plotDat$variable, levels =  c("numORF_asARG", "numORF_multidrug", "numORF_MRG","numORF_multimetal"))
library(ggpubr)

p<-ggplot(plotDat) + 
  geom_violin(aes(x=pathogenicity, y=value, fill=pathogenicity), alpha=0.6, trim = F) +
  scale_fill_manual(values = c("darkslategray3", "goldenrod3")) +
  facet_grid(.~variable) +
  geom_boxplot(aes(x=pathogenicity, y=value, fill=pathogenicity),width=0.05) +
  theme_bw() + theme(panel.background = element_blank(),
                     panel.grid = element_blank()) 


p
ggsave(p, filename = "FigS19.violin_pathogen.other.pdf", device = "pdf", width = 12, height = 3)


p<-ggviolin(plotDat.l, x = "pathogenicity", y = "value",
            color = "pathogenicity", palette = "jco",
            add = "mean_sd", #add.params = list(width  = "0.3"),
            facet.by = "variable", short.panel.labs = FALSE)+
  stat_compare_means( label = "p.format")
p
