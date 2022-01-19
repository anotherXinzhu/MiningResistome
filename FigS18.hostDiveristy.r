setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")

library(data.table)

dat.original <- fread("_data/nonRegARG.bins.txt", data.table = F)
load("_data/otherARGtypes.RData")


# plotting for total numbers (#Bins, #Phylum, #Family, #Genus)-------------------------------
library(ggsci)
library(scales)
show_col(pal_jama("default")(7));pal_jama("default")(7)
show_col(pal_jama("default", alpha = 0.7)(7));pal_jama("default", alpha = 0.7)(7)
myPallette6 <- c("#0073C299", "#EFC00099",  "#CD534C99", "#79AF97B2","#6A6599B2","#DF8F44B2")


for(ecosystem in c("AMDSedi","Tailing")){
  
  # ecosystem = "Tailing"
  
  dat <- dat.original %>% filter(grepl(ecosystem, binFileName))
  
  plotDat = merge(merge(dat %>% select(drugType, binFileName) %>% unique() %>% group_by(drugType) %>% summarise(numBins=n()),
                        dat %>% select(drugType, Phylum) %>% unique() %>% group_by(drugType) %>% summarise(numPhylum=n()),
                        by = "drugType"),
                  merge(dat %>% select(drugType, Family) %>% unique() %>% group_by(drugType) %>% summarise(numFamily=n()),
                        dat %>% select(drugType, Genus) %>% unique() %>% group_by(drugType) %>% summarise(numGenus=n()),
                        by = "drugType"),
                  by="drugType") %>% reshape2::melt(id.var="drugType")
  
  x_rankDf <- dat %>% select(drugType, binFileName) %>% unique() %>% group_by(drugType) %>% summarise(numBins=n()) %>% arrange(desc(numBins))
  plotDat$drugType <- factor(plotDat$drugType, levels = x_rankDf$drugType)
  x_label <- sapply(x_rankDf$drugType, 
                    function(x) if(x %in% otherARGtypes) paste(x, "(*)", sep = "") else x)
  
  p<-ggplot(plotDat, aes(x=drugType, y=value)) +
    geom_col(aes(fill=variable), position = position_dodge(),width = 0.6) +
    scale_fill_manual(values = myPallette6) +
    facet_wrap(vars(variable),ncol = 1,nrow = 4, scales = "free_y")+
    theme_bw() + theme(panel.grid = element_blank()) + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position = "none") +
    scale_x_discrete(labels= x_label)
  
  pdfName <- paste("FigS16.HostDiversiyPlot_",ecosystem,".pdf",sep = "")
  
  ggsave(p, filename = pdfName, device = "pdf", width = 7, height = 10)
  
}
