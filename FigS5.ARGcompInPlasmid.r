library(data.table)

dat <- fread("_data/nonRegARG_onPlasmid.txt")

dat <- dat %>% 
  mutate(dataColl =
           sapply(sample_name,function(x){
             if(grepl("^LD",x,perl = T)) "AMDSedi" else if(grepl("tailing",x)) "tailings" else "pubAMD" 
           }))

load("_data/otherARGtypes.RData")
dat$drugType[dat$drugType %in% otherARGtypes] <- "other"

plotDat <- dat %>% filter(drugType != 'unclassified')  %>%
  group_by(dataColl, drugType) %>% summarise(count= n()) %>% mutate(freq=count/sum(count))

plotDat$dataColl <- factor(plotDat$dataColl, levels = c("pubAMD","tailings","AMDSedi"))

library(ggsci)
drugType_colorsDf <- cbind.data.frame(
  Colors=c(pal_npg("nrc")(10), "gray","gold3","orchid"),
  DrugTypes=c("aminoglycoside", "bacitracin", "beta-lactam","glycopeptide",
              "MLS","tetracycline","phenicol","multidrug",
              "mupirocin","sulfonamide","other","fosmidomycin","rifamycin"), 
  stringsAsFactors=F
)

plotDat$drugType <- factor(plotDat$drugType,
                           c(unique(plotDat$drugType)[unique(plotDat$drugType)!="other"],"other") )

Colors <- sapply(levels(plotDat$drugType), 
                 function(x) drugType_colorsDf$Colors[drugType_colorsDf$DrugTypes == x] )

library(ggplot2)

P_overall <- ggplot(plotDat) +
  geom_col(aes(x=dataColl, y=freq, fill=drugType)) + 
  scale_fill_manual(values = Colors) +
  theme_bw() + theme(panel.grid = element_blank())+ 
  theme ( axis.text.x = element_text ( angle = 90, hjust = 1 )) 

P_overall

ggsave(P_overall, filename = "FigS3.ARGcomposition.plasmids.pdf", device = "pdf")
