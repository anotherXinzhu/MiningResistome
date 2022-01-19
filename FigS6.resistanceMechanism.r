library(dplyr)
library(data.table)


ARG_df.l <- fread("_data/nonRegARG.abund.csv") %>%
  mutate(dataColl = sapply(sampleID,
                           function(x){
                             if(grepl("^LD", x, perl = T)) "AMDSedi" else if(grepl("\\.T",x,perl=T)) "Tailings" else "pubAMD"
                           }))

plotDat_pubAMD <- ARG_df.l %>% filter(dataColl == "pubAMD") %>%
  group_by(sampleID, res_mechanism) %>% summarise(DepthPG=sum(DepthPG)) %>% 
  group_by(res_mechanism) %>% summarise(DepthPG = mean(DepthPG)) %>% mutate(Freq= DepthPG/sum(DepthPG))




plotDat_Tailing <- ARG_df.l %>% filter(dataColl == "Tailings") %>% 
  group_by(sampleID, res_mechanism) %>% summarise(DepthPG=sum(DepthPG)) %>% 
  group_by(res_mechanism) %>% summarise(DepthPG = mean(DepthPG)) %>% mutate(Freq= DepthPG/sum(DepthPG))




plotDat_AMDSedi <- ARG_df.l %>% filter(dataColl == "AMDSedi") %>%
  group_by(sampleID, res_mechanism) %>% summarise(DepthPG=sum(DepthPG)) %>% 
  group_by(res_mechanism) %>% summarise(DepthPG = mean(DepthPG)) %>% mutate(Freq= DepthPG/sum(DepthPG))



pdf("FigS4.resMech/pubADM.pdf")
pie(plotDat_pubAMD$Freq,labels = plotDat_pubAMD$res_mechanism,
    col = c( "#FED439FF", "#709AE1FF", "#FD7446FF", "#D5E4A2FF", "#197EC0FF","#46732EFF", "#D2AF81FF"))
dev.off()

pdf("FigS4.resMech/Tailing.pdf")
pie(plotDat_Tailing$Freq,labels = plotDat_Tailing$res_mechanism,
    col = c( "#FED439FF", "#709AE1FF", "#FD7446FF", "#D5E4A2FF", "#197EC0FF","#46732EFF", "#D2AF81FF"))
dev.off()



pdf("FigS4.resMech/AMDSedi.pdf")
pie(plotDat_AMDSedi$Freq,labels = plotDat_AMDSedi$res_mechanism,
    col = c( "#FED439FF", "#709AE1FF", "#FD7446FF", "#D5E4A2FF", "#197EC0FF","#46732EFF", "#D2AF81FF"))
dev.off()
