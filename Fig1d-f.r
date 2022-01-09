library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer) 
library(eulerr)


# read data ---------------------
load("_data/1.GeoInformaton.RData") # geographic information
allARG_df.l <- fread("_data/nonRegARG.abund.csv")

allARGtype_df.l <- 
  allARG_df.l %>%
  group_by(sampleID, drug_type) %>%
  summarise(DepthPG=sum(DepthPG), sampleTypes=unique(sampleTypes))


# bar plot for ARG drug type composition -----------------------------
Drug_depth_df <- 
  allARGtype_df.l %>% 
  dplyr::select(drug_type,sampleID,DepthPG,sampleTypes)  %>%
  filter(drug_type != "unclassified") 


soil_geo_df$sampleID <- as.character(soil_geo_df$sampleID)
soil_geo_df$location_Eng <- as.character(soil_geo_df$location_Eng)


Drug_depth_df$location <-
  sapply(Drug_depth_df$sampleID, 
         function(x) {
           if(grepl("^LD",x,perl = T)){
             AMDsedi_geo_df$location_Eng[which(AMDsedi_geo_df$sampleID == x)]
           }else if(grepl("^pubAMD",x,perl = T)){
             sp=strsplit(x,"_",fixed = T)[[1]][2]
             publicAMD_info_df$location_Eng[which(publicAMD_info_df$V2 == sp)]
           }else{
             soil_geo_df$location_Eng[which(soil_geo_df$sampleID==x)]
           }})


totalARG_depth_df <- 
  Drug_depth_df %>% 
  dplyr::group_by(sampleID) %>% 
  dplyr::summarise(DepthPG=sum(DepthPG)) %>%
  as.data.frame()

Drug_depth_df$totalARG_depth <- sapply(Drug_depth_df$sampleID,
                                       function(x) totalARG_depth_df$DepthPG[which(totalARG_depth_df$sampleID == x)]) 
Drug_depth_df$Drug_percent <- Drug_depth_df$DepthPG/Drug_depth_df$totalARG_depth



#ARG types beyond top 10 are other ARG types 
ARG.relAbund_df <-
  Drug_depth_df %>% 
  dplyr::group_by(drug_type) %>% dplyr::summarise(meanPercent = mean(Drug_percent)) %>%
  as.data.frame() %>% 
  arrange(desc(meanPercent))
otherARGtypes <- ARG.relAbund_df$drug_type[11:nrow(ARG.relAbund_df)]
save(ARG.relAbund_df, otherARGtypes, file= "_data/otherARGtypes.RData")


# replace other ARG types with others ------------------
Drug_depth_df$drug_type[Drug_depth_df$drug_type %in% otherARGtypes] <- "other"
Drug_depth_df <- Drug_depth_df %>% 
  group_by(sampleID,drug_type) %>%
  summarise(DepthPG=sum(DepthPG), sampleTypes=unique(sampleTypes), 
            location=unique(location), Drug_percent=sum(Drug_percent))

plotDat <- Drug_depth_df %>% 
  group_by(location, sampleTypes, drug_type) %>% 
  dplyr::summarise(Drug_percent=mean(Drug_percent)) %>% as.data.frame() 
plotDat$DataCollection <- sapply(plotDat$location,
                                 function(x)all_geoBySite_df$DataCollection[which(all_geoBySite_df$location_Eng==x)[1]])
plotDat$latitude <- sapply(plotDat$location,
                           function(x)all_geoBySite_df$latitude[which(all_geoBySite_df$location_Eng==x)[1]])


# colors ---------------------
library(scales)
library(ggsci)

drugType_colorsDf <- cbind.data.frame(
  Colors=c(pal_npg("nrc")(10), "gray","gold3","orchid"),
  DrugTypes=c("aminoglycoside", "bacitracin", "beta-lactam","glycopeptide",
              "MLS","tetracycline","phenicol","multidrug",
              "mupirocin","sulfonamide","other","fosmidomycin","rifamycin"), 
  stringsAsFactors=F
)

drugType_rankABC <- unique(plotDat$drug_type)[order(unique(plotDat$drug_type))]
drugType_fctLevel <- c(drugType_rankABC[drugType_rankABC != "other"], "other") #put other at the end
drugType_pallette <-sapply(drugType_fctLevel, 
                           function(x) drugType_colorsDf$Colors[which(drugType_colorsDf$DrugTypes == x)])


plotDat_pubAMD <- plotDat %>% 
  dplyr::filter(DataCollection=="pubAMD") %>%
  dplyr::filter(!is.na(latitude)) %>%
  mutate(location_spType = paste(location, sampleTypes,sep = "&")) %>%
  arrange(latitude)# %>% arrange(sampleTypes)


# integrate data of European samples with no reads 
total_Argdf_European <- fread("_data/nonRegARGcount.EuropeanNoReads.csv")

ArgMapping_tmp <- allARG_df.l %>% select(ARG,drug_type) %>% unique()
total_Argdf_European$drug_type <- sapply(total_Argdf_European$ARG,
                                         function(x) ArgMapping_tmp$drug_type[which(ArgMapping_tmp$ARG == x)])

total_Argdf_European <- total_Argdf_European %>% filter(drug_type != "unclassified") 
total_Argdf_European$drug_type[which(total_Argdf_European$drug_type %in% otherARGtypes)] <- "other"
plotDat_pubAMD2 <- merge(total_Argdf_European %>% group_by(location, drug_type) %>% summarise(numORF = sum(num)),
                         total_Argdf_European %>% group_by(location) %>% summarise(totalNumORF = sum(num)),
                         by='location') %>%
  mutate(Drug_percent = numORF/totalNumORF) %>%
  select(location, drug_type,Drug_percent)

plotDat_pubAMD2$sampleTypes = 
  sapply(plotDat_pubAMD2$location,
         function(x) unique(publicAMD_info_df$ecosystem[which(publicAMD_info_df$location_Eng == x)]))  

plotDat_pubAMD2$latitude <-
  sapply(plotDat_pubAMD2$location,
         function(x) unique(publicAMD_info_df$latitude[which(publicAMD_info_df$location_Eng == x)]))

plotDat_pubAMD2$DataCollection <- rep("pubAMD", nrow(plotDat_pubAMD2))
plotDat_pubAMD2$location_spType <- paste(plotDat_pubAMD2$location,plotDat_pubAMD2$sampleTypes,sep = "&")


plotDat_pubAMD <- bind_rows(plotDat_pubAMD, plotDat_pubAMD2) 
plotDat_pubAMD <- plotDat_pubAMD %>% arrange(latitude)

plotDat_pubAMD$x_fct <- factor(plotDat_pubAMD$location_spType,
                               levels = unique(plotDat_pubAMD$location_spType))
plotDat_pubAMD$drugType_fct <- factor(plotDat_pubAMD$drug_type,
                                      levels = drugType_fctLevel)

ARGCompP_pubAMD <- ggplot(plotDat_pubAMD, aes(x = x_fct, y = Drug_percent, fill = drugType_fct)) + 
  geom_col() +  
  scale_fill_manual(values = drugType_pallette) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  xlab("") + ylab("percentage of drug types") + 
  theme(legend.position = "right") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  scale_x_discrete(labels= sub("(PubSite_\\d+)\\&\\w+$","\\1",unique(plotDat_pubAMD$location_spType))) 
ARGCompP_pubAMD



plotDat_AMDSedi <- plotDat %>% dplyr::filter(DataCollection=="AMDSedi") %>% arrange(latitude) 
plotDat_AMDSedi$location <- factor(plotDat_AMDSedi$location,levels = unique(plotDat_AMDSedi$location))
plotDat_AMDSedi$drugType_fct <- factor(plotDat_AMDSedi$drug_type,
                                       levels = drugType_fctLevel)

ARGCompP_AMDSedi <- ggplot(plotDat_AMDSedi, aes(x = location, y = Drug_percent, fill = drugType_fct)) + 
  geom_col() + 
  scale_fill_manual(values = drugType_pallette) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  xlab("") + ylab("percentage of drug types") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
ARGCompP_AMDSedi



plotDat_Tail <- plotDat %>% dplyr::filter(DataCollection=="National") %>% arrange(latitude) 
plotDat_Tail$location <- factor(plotDat_Tail$location,levels = unique(plotDat_Tail$location))
plotDat_Tail$drugType_fct <- factor(plotDat_Tail$drug_type,
                                    levels = drugType_fctLevel)
ARGCompP_Tail <- ggplot(plotDat_Tail, aes(x = location, y = Drug_percent, fill = drugType_fct)) + 
  geom_col() +  
  scale_fill_manual(values = drugType_pallette) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  xlab("") + ylab("percentage of drug types") + guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(labels= sub("_Tail","",unique(plotDat_Tail$location))) 
ARGCompP_Tail




library(ggpubr)
pdf("Fig1d-f.ARGcomposition.pdf",width = 15,height = 5)
ggarrange(ARGCompP_pubAMD+theme(legend.position = "none"),
          ARGCompP_Tail+theme(legend.position = "none"),
          ARGCompP_AMDSedi,
          nrow = 1,ncol = 3,widths = c(2.1,4.7,3.0),heights = c(0.7,1,0.7))
dev.off()

