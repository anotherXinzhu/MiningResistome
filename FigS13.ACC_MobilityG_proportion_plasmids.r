setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")

load("_data/nonRegACCproportions_plasmids.RData")


library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

# pubAMD #####
dat_pubAMD <- numContig_df_pubAMD %>% 
  filter(contigType == "plasmid") %>%
  filter(!is.na(totalContigDepth)) %>%
  group_by(ecosystem,location_Eng) %>% 
  summarise(latitude=mean(latitude), n=n(),
            avg_ACCpercent_byNum=mean(ACCpercent_byNum,na.rm = T), sd_ACCpercent_byNum=sd(ACCpercent_byNum,na.rm = T),
            avg_ACCpercent_byDepth=mean(ACCpercent_byDepth,na.rm = T),sd_ACCpercent_byDepth=sd(ACCpercent_byDepth,na.rm = T),
            avg_ACC_MetalRG_Percent_byNum=mean(ACC_MetalRG_Percent_byNum,na.rm = T), sd_ACC_MetalRG_Percent_byNum=sd(ACC_MetalRG_Percent_byNum,na.rm = T),
            avg_ACC_MetalRG_Percent_byDepth=mean(ACC_MetalRG_percent_byDepth,na.rm = T),sd_ACC_MetalRG_Percent_byDepth=sd(ACC_MetalRG_percent_byDepth,na.rm = T),
            avg_ACC_MGEnzyme_Percent_byNum=mean(ACC_MGEnzyme_Percent_byNum,na.rm = T),sd_ACC_MGEnzyme_Percent_byNum=sd(ACC_MGEnzyme_Percent_byNum,na.rm = T),
            avg_ACC_MGEnzyme_Percent_byDepth=mean(ACC_MGEnzyme_Percent_byDepth,na.rm = T),sd_ACC_MGEnzyme_Percent_byDepth=sd(ACC_MGEnzyme_Percent_byDepth,na.rm = T)) %>%
  as.data.frame() 

dat_pubAMD_tmp1 <- dat_pubAMD %>% select(all_of(c("ecosystem", "location_Eng", "latitude", "n", colnames(dat_pubAMD)[grepl("[Pp]ercent_byDepth", colnames(dat_pubAMD),perl = T)])))
colnames(dat_pubAMD_tmp1) <- sub("^(.*[Pp]ercent)(_by\\w+)$", "\\1",colnames(dat_pubAMD_tmp1))


# European data 

dat_pubAMD_Europe <- numContig_df_pubAMD %>% filter(contigType == "plasmid") %>%
  group_by(ecosystem,location_Eng) %>% 
  summarise(latitude=mean(latitude), n=n(),
            avg_ACCpercent_byNum=mean(ACCpercent_byNum,na.rm = T), sd_ACCpercent_byNum=sd(ACCpercent_byNum,na.rm = T),
            avg_ACCpercent_byDepth=mean(ACCpercent_byDepth,na.rm = T),sd_ACCpercent_byDepth=sd(ACCpercent_byDepth,na.rm = T),
            avg_ACC_MetalRG_Percent_byNum=mean(ACC_MetalRG_Percent_byNum,na.rm = T), sd_ACC_MetalRG_Percent_byNum=sd(ACC_MetalRG_Percent_byNum,na.rm = T),
            avg_ACC_MetalRG_Percent_byDepth=mean(ACC_MetalRG_percent_byDepth,na.rm = T),sd_ACC_MetalRG_Percent_byDepth=sd(ACC_MetalRG_percent_byDepth,na.rm = T),
            avg_ACC_MGEnzyme_Percent_byNum=mean(ACC_MGEnzyme_Percent_byNum,na.rm = T),sd_ACC_MGEnzyme_Percent_byNum=sd(ACC_MGEnzyme_Percent_byNum,na.rm = T),
            avg_ACC_MGEnzyme_Percent_byDepth=mean(ACC_MGEnzyme_Percent_byDepth,na.rm = T),sd_ACC_MGEnzyme_Percent_byDepth=sd(ACC_MGEnzyme_Percent_byDepth,na.rm = T)) %>%
  as.data.frame() %>%
  filter(location_Eng %in% c("PubSite_10","PubSite_11","PubSite_14")) %>% 
  select_if(function(x){!all(is.na(x))}) # keep columns which are not all-NA

colnames(dat_pubAMD_Europe) <- sub("^(.*[Pp]ercent)(_by\\w+)$", "\\1",colnames(dat_pubAMD_Europe))

# integrate two pub_AMD dfs 

dat_pubAMD_combined <- bind_rows(dat_pubAMD_tmp1, dat_pubAMD_Europe) %>% arrange(latitude)


# arrange x
dat_pubAMD_combined$site_ecosystem <-
  paste(dat_pubAMD_combined$location_Eng,dat_pubAMD_combined$ecosystem,sep="|")
dat_pubAMD_combined <- dat_pubAMD_combined %>% filter( !is.na(latitude)) %>% arrange(site_ecosystem) %>% arrange(latitude)
dat_pubAMD_combined$site_ecosystem <- factor(dat_pubAMD_combined$site_ecosystem,
                                             levels = dat_pubAMD_combined$site_ecosystem)


# plot by variable
library(ggplot2)
spType_colorDf <- cbind.data.frame(spType=c("AMD","AMDSediment","biofilm","Tailing"),
                                   Color=c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
                                   stringsAsFactors=F)

y_vec <- c("ACC_MGEnzyme_Percent")
for(y in y_vec){
  
  # y=y_vec[2]
  i=which(y_vec==y)
  x="site_ecosystem"
  dat <- dat_pubAMD_combined %>% 
    select(ecosystem,location_Eng,all_of(c(x, colnames(dat_pubAMD_combined)[grep(y,colnames(dat_pubAMD_combined))]))) 
  
  colnames(dat)[grep("^avg",colnames(dat))] <- "Avg_proportion"
  colnames(dat)[grep("^sd",colnames(dat))] <- "sd"
  colnames(dat)[which(colnames(dat)==x)] <- "x"
  
  SpTypeABC <- unique(dat$ecosystem)[order(unique(dat$ecosystem))]
  Colors <- sapply(SpTypeABC, function(x) spType_colorDf$Color[which(spType_colorDf$spType == x)])
  
  
  p  <- ggplot(dat)+
    geom_col(aes(x=x, y=Avg_proportion,fill=ecosystem),width = 0.6,color="black")  +
    geom_errorbar(aes(ymin = Avg_proportion - sd, ymax = Avg_proportion + sd,x = x),  width = 0.2, position = position_dodge(0.9)) +
    ylim(c(-0.1,0.7)) +  # for ARG-MRG proportion
    #ylim(c(0,1)) +  # for ARG-MGE enzyme proportion
    scale_fill_manual(values = Colors)  +
    scale_x_discrete(labels=dat$location_Eng) + #改变x轴的label
    xlab("") + ylab(y)+
    theme_bw() + 
    theme ( axis.text.x = element_text (angle = 45, hjust = 1 ), legend.position = 'top') +
    theme(panel.grid.minor =element_blank(), panel.grid.major = element_blank()) 
  
  p
  assign(paste("py",i,"_pubAMD",sep = ""),p,envir = .GlobalEnv)
  
}
py1_pubAMD


# AMD sediment #####
dat_AMDSedi <- numContig_df_AMDSedi %>% filter(contigType == "plasmid") %>%
  filter(!is.na(totalContigDepth)) %>%
  group_by(location_Eng) %>% 
  summarise(latitude=mean(latitudes), n=n(),
            avg_ACCpercent_byNum=mean(ACCpercent_byNum), sd_ACCpercent_byNum=sd(ACCpercent_byNum),
            avg_ACCpercent_byDepth=mean(ACCpercent_byDepth),sd_ACCpercent_byDepth=sd(ACCpercent_byDepth),
            avg_ACC_MetalRG_Percent_byNum=mean(ACC_MetalRG_Percent_byNum), sd_ACC_MetalRG_Percent_byNum=sd(ACC_MetalRG_Percent_byNum),
            avg_ACC_MetalRG_Percent_byDepth=mean(ACC_MetalRG_percent_byDepth),sd_ACC_MetalRG_Percent_byDepth=sd(ACC_MetalRG_percent_byDepth),
            avg_ACC_MGEnzyme_Percent_byNum=mean(ACC_MGEnzyme_Percent_byNum),sd_ACC_MGEnzyme_Percent_byNum=sd(ACC_MGEnzyme_Percent_byNum),
            avg_ACC_MGEnzyme_Percent_byDepth=mean(ACC_MGEnzyme_Percent_byDepth),sd_ACC_MGEnzyme_Percent_byDepth=sd(ACC_MGEnzyme_Percent_byDepth)) %>%
  as.data.frame()

# arrange x
dat_AMDSedi <- dat_AMDSedi %>% filter( !is.na(latitude)) %>% arrange(location_Eng) %>% arrange(latitude)
dat_AMDSedi$location_Eng <- factor(dat_AMDSedi$location_Eng,levels = dat_AMDSedi$location_Eng)

# plot by variable

DataCollection = "AMDSediment"
y_vec <- c("ACC_MGEnzyme_Percent_byDepth")
for(y in y_vec){
  
  # y=y_vec[2]
  i=which(y_vec==y)
  x="location_Eng"
  dat <- dat_AMDSedi %>% select(all_of(c(x,colnames(dat_pubAMD)[grep(y,colnames(dat_pubAMD))]))) 
  
  colnames(dat)[grep("^avg",colnames(dat))] <- "Avg_proportion"
  colnames(dat)[grep("^sd",colnames(dat))] <- "sd"
  colnames(dat)[which(colnames(dat)==x)] <- "x"
  
  Colors = spType_colorDf$Color[which(spType_colorDf$spType=="AMDSediment")]
  
  p  <- ggplot(dat)+
    geom_col(aes(x=x, y=Avg_proportion),width = 0.6,color="black",fill=Colors)  +
    geom_errorbar(aes(ymin = Avg_proportion - sd, ymax = Avg_proportion + sd,x = x),  width = 0.2, position = position_dodge(0.9)) +
    xlab("") + ylab(y)+
    #ylim(c(0,0.5)) +  # for ARG-MRG proportion
    ylim(c(-0.1,0.7)) +  # for ARG-MGE enzyme proportion
    theme_bw() + 
    theme ( axis.text.x = element_text (angle = 45, hjust = 1 ), legend.position = 'top') +
    theme(panel.grid.minor =element_blank(), panel.grid.major = element_blank()) 
  
  p
  assign(paste("py",i,"_ADMSedi",sep = ""),p,envir = .GlobalEnv)
  
}
py1_ADMSedi  #  ylim(c(-0.025,0.25))


# Tailing #####
dat_Tailing <- numContig_df_Tailing %>% filter(contigType == "plasmid") %>%
  filter(!is.na(totalContigDepth)) %>%
  group_by(location_Eng) %>% 
  summarise(latitude=mean(latitude), n=n(),
            avg_ACCpercent_byNum=mean(ACCpercent_byNum), sd_ACCpercent_byNum=sd(ACCpercent_byNum),
            avg_ACCpercent_byDepth=mean(ACCpercent_byDepth),sd_ACCpercent_byDepth=sd(ACCpercent_byDepth),
            avg_ACC_MetalRG_Percent_byNum=mean(ACC_MetalRG_Percent_byNum), sd_ACC_MetalRG_Percent_byNum=sd(ACC_MetalRG_Percent_byNum),
            avg_ACC_MetalRG_Percent_byDepth=mean(ACC_MetalRG_percent_byDepth),sd_ACC_MetalRG_Percent_byDepth=sd(ACC_MetalRG_percent_byDepth),
            avg_ACC_MGEnzyme_Percent_byNum=mean(ACC_MGEnzyme_Percent_byNum),sd_ACC_MGEnzyme_Percent_byNum=sd(ACC_MGEnzyme_Percent_byNum),
            avg_ACC_MGEnzyme_Percent_byDepth=mean(ACC_MGEnzyme_Percent_byDepth),sd_ACC_MGEnzyme_Percent_byDepth=sd(ACC_MGEnzyme_Percent_byDepth)) %>%
  as.data.frame()

# arrange x
dat_Tailing <- dat_Tailing %>% filter( !is.na(latitude)) %>% arrange(location_Eng) %>% arrange(latitude)
dat_Tailing$location_Eng <- factor(dat_Tailing$location_Eng,levels = dat_Tailing$location_Eng)

# plot by variable

DataCollection = "Tailing"
y_vec <- c("ACC_MGEnzyme_Percent_byDepth")
for(y in y_vec){
  
  # y=y_vec[1]
  i=which(y_vec==y)
  x="location_Eng"
  dat <- dat_Tailing %>% select(all_of(c(x,colnames(dat_pubAMD)[grep(y,colnames(dat_pubAMD))]))) 
  
  colnames(dat)[grep("^avg",colnames(dat))] <- "Avg_proportion"
  colnames(dat)[grep("^sd",colnames(dat))] <- "sd"
  colnames(dat)[which(colnames(dat)==x)] <- "x"
  
  Colors = spType_colorDf$Color[which(spType_colorDf$spType == "Tailing")]
  
  p  <- ggplot(dat)+
    geom_col(aes(x=x, y=Avg_proportion),width = 0.6,color="black",fill=Colors)  +
    geom_errorbar(aes(ymin = Avg_proportion - sd, ymax = Avg_proportion + sd,x = x),  width = 0.2, position = position_dodge(0.9)) +
    xlab("") + ylab(y)+
    #yylim(c(0,0.5)) +  # for ARG-MRG proportion
    ylim(c(-0.1,0.7)) +  # for ARG-MGE enzyme proportion
    scale_x_discrete(labels=sub("(.*)_Tail","\\1",dat$x)) + 
    theme_bw() + 
    theme ( axis.text.x = element_text (angle = 45, hjust = 1 ), legend.position = 'top') +
    theme(panel.grid.minor =element_blank(), panel.grid.major = element_blank()) 
  
  p
  assign(paste("py",i,"_Tailing",sep = ""),p,envir = .GlobalEnv)
  
}
py1_Tailing

# plot together

library(ggpubr)
#ggarrange(py1_pubAMD,py1_Tailing,py1_ADMSedi,nrow = 3,ncol = 1,heights = c(3.7,3.3,3))
pdf("FigS12.ACC-MG_proportion_plasmids.pdf",width = 15,height = 3)
ggarrange(py1_pubAMD + theme(legend.position = "none"),py1_Tailing,py1_ADMSedi,nrow = 1,ncol = 3, widths = c(2.35,4.5,3.05))
dev.off()


