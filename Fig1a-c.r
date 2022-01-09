

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)

# load data ------------------------------------------------------------
load("_data/1.GeoInformaton.RData") # geographic information

allARG_df.l <- fread("_data/nonRegARG.abund.csv")
colnames(allARG_df.l)[which(colnames(allARG_df.l) == "DepthPG")] <- "Abundance"


# plot for public AMD data  ------------------------------------------------------------

dat_pubAMD <-  allARG_df.l %>% 
  dplyr::filter(grepl("^pubAMD", sampleID, perl = T))  %>% 
  dplyr::group_by(sampleID) %>% 
  dplyr::summarise(totalARG = sum(Abundance),numSubtype = sum(Abundance>0)) %>%
  as.data.frame() 


dat_pubAMD$shortID <- sapply(strsplit(dat_pubAMD$sampleID,"_",fixed = T),"[[",2)
dat_pubAMD$location <- sapply(dat_pubAMD$shortID, 
                              function(x) publicAMD_info_df$location[which(publicAMD_info_df$V2 == x)])
dat_pubAMD$ecosystem <- sapply(dat_pubAMD$shortID, 
                               function(x) publicAMD_info_df$ecosystem[which(publicAMD_info_df$V2 == x)])


dat_pubAMD_bySite <- dat_pubAMD %>% 
  dplyr::group_by(location, ecosystem) %>% 
  dplyr::summarise(ARG_Abund =mean(totalARG), ARG_sd=sd(totalARG), 
                   avg_numSubtype = mean(numSubtype), numSubtype_sd = sd(numSubtype)) %>%
  as.data.frame() %>% base::merge(all_geoBySite_df,by.x=c("location","ecosystem"),by.y=c("location_Eng","ecosystem")) %>%
  dplyr::arrange(latitude) %>% dplyr::filter(!is.na(latitude))


# European data without reads ####
total_Argdf_European <- fread("_data/ARGcount.EuropeanNoReads.csv")



tmp1 <- total_Argdf_European %>% 
  group_by(sample) %>% summarise(numARGsubtype=n(), location=unique(location)) %>%
  group_by(location) %>% summarise(avg_numSubtype=mean(numARGsubtype), numSubtype_sd=sd(numARGsubtype))
tmp2 <- total_Argdf_European %>% select(location,ecosystem,longitude,latitude,`latitude and longitude`) %>% unique() 
tmp3 <- total_Argdf_European %>% select(location, sample) %>% unique() %>% group_by(location) %>% summarise(n=n())

plotDat_pubAMD2  <- merge(merge(tmp1,tmp2,by="location"), tmp3, by="location")
plotDat_pubAMD2$location 
plotDat_pubAMD2$position <- c("Spain","France","Sweden") #check manually 
plotDat_pubAMD2$DataCollection <- "pubAMD"

# integrate dat ####
dat_pubAMD_bySite <- bind_rows(dat_pubAMD_bySite, plotDat_pubAMD2)

dat_pubAMD_bySite$location_ecosystem <- paste(dat_pubAMD_bySite$location,dat_pubAMD_bySite$ecosystem,sep = "|")
dat_pubAMD_bySite <- dat_pubAMD_bySite %>% arrange(latitude)
dat_pubAMD_bySite$location_ecosystem <- factor(dat_pubAMD_bySite$location_ecosystem,
                                               levels = dat_pubAMD_bySite$location_ecosystem)



# colors for ecosystem #####
spType_colorDf <- cbind.data.frame(spType=c("AMD","AMDSediment","biofilm","Tailing"),
                                   Color=c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
                                   stringsAsFactors=F)

Colors=sapply(unique(dat_pubAMD_bySite$ecosystem),
              function(x) spType_colorDf$Color[which(spType_colorDf$spType == x)])


p1 <- ggplot(dat_pubAMD_bySite) +
  geom_col(aes(x=location_ecosystem, y=ARG_Abund,fill=ecosystem),width = 0.6,color="black")  +
  scale_fill_manual(values=Colors) +
  geom_errorbar(aes(ymin = ARG_Abund - ARG_sd, ymax = ARG_Abund + ARG_sd,x = location_ecosystem),  width = 0.2, position = position_dodge(0.9)) +
  theme_bw() + theme(panel.grid.minor =element_blank(), panel.grid.major = element_blank()) +
  scale_x_discrete(labels= dat_pubAMD_bySite$location)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  ylim(0, 4100) +
  theme(legend.position="none") 
p1


p2 <- ggplot(dat_pubAMD_bySite) + 
  geom_point(aes(x=location_ecosystem, y=avg_numSubtype),size=4,shape=21,fill="white") + #如果用numDrugType就在这里改
  scale_y_continuous(limits = c(0,510), position = "right")  +
  theme_half_open(11, rel_small = 1) + 
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") 
p2


#put two graphs on one layer
#library(ggpubr)
#library(cowplot)
aligned_plots <- align_plots(p1, p2, align="hv", axis="tblr")
p_pubAMD <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p_pubAMD



# plot for AMD sediment -----------------------------------------------------------------------------------------
dat_AMDSedi <- allARG_df.l %>%  
  dplyr::filter(grepl("^LD", sampleID, perl = T))  %>% 
  dplyr::group_by(sampleID) %>% 
  dplyr::summarise(totalARG = sum(Abundance),numSubtype = sum(Abundance>0)) %>%
  as.data.frame() 


dat_AMDSedi$location_Eng <- sapply(dat_AMDSedi$sampleID, 
                                   function(x) AMDsedi_geo_df$location_Eng[which(AMDsedi_geo_df$sampleID == x)])  
dat_AMDSedi$ecosystem <- rep("AMDSediment",nrow(dat_AMDSedi))


dat_AMDSedi_bySite <- dat_AMDSedi %>% dplyr::group_by(location_Eng,ecosystem) %>% 
  dplyr::summarise(ARG_Abund =mean(totalARG), ARG_sd=sd(totalARG), 
                   avg_numSubtype = mean(numSubtype), numSubtype_sd = sd(numSubtype)) %>%
  as.data.frame() %>% base::merge(all_geoBySite_df,by=c("location_Eng","ecosystem")) %>%
  dplyr::arrange(latitude) %>% dplyr::filter(!is.na(latitude))


dat_AMDSedi_bySite$location_Eng <- factor(dat_AMDSedi_bySite$location_Eng,
                                          levels = dat_AMDSedi_bySite$location_Eng)

Colors=sapply(unique(dat_AMDSedi_bySite$ecosystem),function(x) spType_colorDf$Color[which(spType_colorDf$spType == x)])

p1 <- ggplot(dat_AMDSedi_bySite) +
  geom_col(aes(x=location_Eng, y=ARG_Abund),width = 0.6,color="black",fill=Colors)  +
  geom_errorbar(aes(ymin = ARG_Abund - ARG_sd, ymax = ARG_Abund + ARG_sd,x = location_Eng),  width = 0.2, position = position_dodge(0.9)) +
  theme_bw() + theme(panel.grid.minor =element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  ylim(0, 4100) 




p2 <- ggplot(dat_AMDSedi_bySite) +
  geom_point(aes(x=location_Eng, y=avg_numSubtype),size=4,shape=21,fill="white") +
  scale_y_continuous(limits = c(0,510), position = "right")  +
  theme_half_open(11, rel_small = 1) + 
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") 
p2

aligned_plots <- align_plots(p1, p2, align="hv", axis="tblr")
p_AMDsedi <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p_AMDsedi


# plot for Tailing ----------------------------------------------

dat_Tailing <- allARG_df.l %>% 
  filter(grepl(".T", sampleID, fixed = T))  %>% 
  dplyr::group_by(sampleID) %>% 
  dplyr::summarise(totalARG = sum(Abundance),numSubtype = sum(Abundance>0)) %>% 
  as.data.frame() 

dat_Tailing$location_Eng <- sapply(dat_Tailing$sampleID, 
                                   function(x) soil_geo_df$location_Eng[which(soil_geo_df$sampleID == x)]) 
dat_Tailing$ecosystem <- rep("Tailing",nrow(dat_Tailing))


dat_Tailing_bySite <- dat_Tailing %>% dplyr::group_by(location_Eng,ecosystem) %>% 
  dplyr::summarise(ARG_Abund =mean(totalARG), 
                   ARG_sd=if(length(totalARG)==2) (max(totalARG)-min(totalARG))/2 else sd(totalARG), 
                   avg_numSubtype = mean(numSubtype),
                   numSubtype_sd = sd(numSubtype) ) %>%
  as.data.frame() %>% base::merge(all_geoBySite_df,by=c("location_Eng","ecosystem")) %>%
  dplyr::arrange(latitude) %>% dplyr::filter(!is.na(latitude))



dat_Tailing_bySite$location_Eng <- factor(dat_Tailing_bySite$location_Eng,
                                          levels=dat_Tailing_bySite$location_Eng)

Colors=sapply(unique(dat_Tailing_bySite$ecosystem),function(x) spType_colorDf$Color[which(spType_colorDf$spType == x)])

#library(ggplot2)
p1 <- ggplot(dat_Tailing_bySite) +
  geom_col(aes(x=location_Eng, y=ARG_Abund),width = 0.6,color="black",fill=Colors)  +
  scale_y_continuous(limits = c(0,5000)) +
  geom_errorbar(aes(ymin = ARG_Abund - ARG_sd, ymax = ARG_Abund + ARG_sd,x = location_Eng),  
                width = 0.2, position = position_dodge(0.9)) +
  theme_bw() + theme(panel.grid.minor =element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  ylim(0, 4100) +
  scale_x_discrete(labels= sub("_Tail","",dat_Tailing_bySite$location_Eng)) + rremove("xlab")
p1


p2 <- ggplot(dat_Tailing_bySite) +
  geom_point(aes(x=location_Eng, y=avg_numSubtype),size=4,shape=21,fill="white") +
  scale_y_continuous( limits = c(0,510), position = "right")  +
  #scale_y_continuous(position = "right") +
  theme_half_open(11, rel_small = 1) + 
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") #+
p2

#put two graphs on one layer

aligned_plots <- align_plots(p1, p2, align="hv", axis="tblr")
p_Tailing <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p_Tailing


pdf("Fig1a-c.pdf",width = 15,height = 3)
ggarrange(p_pubAMD,p_Tailing,p_AMDsedi,nrow = 1,ncol = 3, widths = c(2.35,4.5,3.05))
dev.off()
