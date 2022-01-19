setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")

library(data.table)
library(dplyr)

library(ggpubr)
library(cowplot)

ARG_df.l <- fread("_data/nonRegARG.abund.csv")

ARG_df.l <- ARG_df.l %>% 
  group_by(sampleID) %>% 
  mutate(totalARG.abund = sum(DepthPG), relAbund=DepthPG/sum(DepthPG)) %>% as.data.frame() 


#calculate abundant ARGs ---------------------------------------
AvgRelAbund_df <- ARG_df.l %>%
  group_by(dataCollection, ARG) %>% summarise(avgRelAbund = mean(relAbund)) %>% as.data.frame()
 


# abundant ARGs are defined as ARGs with average relative abundance > 0.01
minRA <- 0.01
topARG_pubAMD <- (AvgRelAbund_df %>% 
                    filter(dataCollection == "pubAMD") %>% 
                    filter(avgRelAbund > minRA) %>%
                    arrange(desc(avgRelAbund)))$ARG 
topARG_Tailing <- (AvgRelAbund_df %>%
                     filter(dataCollection == "Tailings") %>% 
                     filter(avgRelAbund > minRA) %>%
                     arrange(desc(avgRelAbund)))$ARG
topARG_AMDSedi <- (AvgRelAbund_df %>%
                     filter(dataCollection == "AMDSedi") %>%
                     filter(avgRelAbund > minRA) %>%
                     arrange(desc(avgRelAbund)))$ARG

# calculate frequent ARGs ---------------------------------
numSp_df <- 
  ARG_df.l %>% 
  select(dataCollection,sampleID) %>% unique() %>% 
  group_by(dataCollection) %>% summarise(numSp=n())

FreqDat <- merge(ARG_df.l %>%
                   group_by(dataCollection,ARG) %>% summarise(presentInNumSample = sum(DepthPG>0)) %>%  
                   as.data.frame(), ##DepthPG>0.02为present
                 numSp_df,
                 by="dataCollection") %>% mutate(Freq = presentInNumSample/numSp)

minFreq = 1

FreqDat %>% group_by(dataCollection) %>% summarise(num_FreqARG = sum(Freq >= minFreq)) 
FreqARG_pubAMD <- (FreqDat %>% 
                     filter(dataCollection == "pubAMD") %>% 
                     filter(Freq >= minFreq))$ARG
FreqARG_AMDSedi <- (FreqDat %>%
                      filter(dataCollection == "AMDSedi") %>%
                      filter(Freq >= minFreq))$ARG
FreqARG_Tailing <- (FreqDat %>% 
                      filter(dataCollection == "Tailings") %>%
                      filter(Freq >= minFreq))$ARG


# FigS5.a-c. plotting bars for abundant ARGs ---------------------------------------
# colors ####
drugType_color_df <- cbind.data.frame(
  Colors=c(pal_npg("nrc")(10), "gray"),
  DrugTypes=c("aminoglycoside", "bacitracin", "beta-lactam","glycopeptide",
              "MLS","tetracycline","phenicol","multidrug",
              "mupirocin","sulfonamide","other"),
  stringsAsFactors=F
) # 

# abbreviations for ARG subtype in x axis-------------------------------
ARGsubtypes <- unique(c(topARG_AMDSedi, topARG_pubAMD, topARG_Tailing))
abb  <- unname(sapply(ARGsubtypes, function(x) {
  splited <- strsplit(x,"_",fixed = T)[[1]]
  splited[which(nchar(splited) == min(nchar(splited)))]
  
}) )

ARGsubtype_abb_df <- cbind.data.frame(ARGsubtypes,abb,stringsAsFactors=F)

for(i in c(1:length(ARGsubtypes))){
  if(ARGsubtype_abb_df$abb[i] == "to"){
    ARGsubtype_abb_df$abb[i] <- "mupA"
  }else if(ARGsubtypes[i] == "EmrB-QacA_family_major_facilitator_transporter"){
    ARGsubtype_abb_df$abb[i] <- "comD"  # reference: https://www.ncbi.nlm.nih.gov/gene/8624036
  }else if(ARGsubtypes[i] == "major_facilitator_superfamily_transporter"){
    ARGsubtype_abb_df$abb[i] <- "MFS"
  }else if(ARGsubtypes[i] == "cob(I)alamin_adenolsyltransferase"){
    ARGsubtype_abb_df$abb[i] <- "cobI_ATR"
  }else if(ARGsubtypes[i] == "chloramphenicol_and_florfenicol_resistance"){
    ARGsubtype_abb_df$abb[i] <- "radical_SAM" #geneBank ACO85172.1, match to  HMM PF04055, manually checked 
  }
}

ARGsubtype_abb_df




# point with error bar -------------------------------------------------
# pubAMD ####
dat_pubAMD <- ARG_df.l %>% filter(dataCollection == "pubAMD") %>% filter(ARG %in% topARG_pubAMD)
dat_pubAMD_sd <- 
  dat_pubAMD %>%
  group_by(ARG) %>% 
  summarise(MRA=mean(relAbund), sd=sd(relAbund), resMech=unique(res_mechanism), drugType=unique(drug_type))


dat_pubAMD$ARG <- factor(dat_pubAMD$ARG, levels=topARG_pubAMD)
p0_1 <-ggplot(dat_pubAMD) + 
  geom_point(aes( x = ARG, y = relAbund),shape=1,alpha=0.6,color="darkgray") +
  theme_bw() +
  labs(y = "relative abundance", x = "", title = "") +
  scale_y_continuous(limits = c(-0.01,0.20),oob=rescale_none) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) +
  scale_x_discrete(labels=sapply(levels(dat_pubAMD$ARG), function(x) ARGsubtype_abb_df$abb[which(ARGsubtype_abb_df$ARGsubtypes == x)] )) +
  theme(axis.text.x = element_text(angle = 90))
p0_1



# color the dots with top two most abundant drugtype and others 
dat_pubAMD_sd$drugType3 <- sapply(dat_pubAMD_sd$drugType,
                                  function(x){
                                    if(!(x %in% c("multidrug","bacitracin"))) "other" else x
                                  })
drugType3.levels <- c("multidrug","bacitracin","other")
dat_pubAMD_sd$drugType3 <- factor(dat_pubAMD_sd$drugType3,levels = drugType3.levels)

Colors= sapply(drugType3.levels,function(x) drugType_color_df$Colors[drugType_color_df$DrugTypes==x])
Colors[3] <- "darkgray"

library(scales)
dat_pubAMD_sd$ARG <- factor(dat_pubAMD_sd$ARG, levels=topARG_pubAMD)
p0_2 <- ggplot(dat_pubAMD_sd) + 
  geom_pointrange(aes(data= ,x = ARG, y = MRA,ymin=MRA-sd, ymax=MRA+sd, color=drugType),size=0.8) +
  scale_color_manual(values = Colors) + 
  scale_y_continuous(limits = c(-0.01,0.20),oob=rescale_none) +
  theme_half_open(11, rel_small = 1) + 
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  # rremove("legend") +
  rremove("y.axis") + rremove("y.text") + rremove("ylab")

p0_2


aligned_plots <- align_plots(p0_1, p0_2+theme(legend.position = "none"), align="hv", axis="tblr")
p_pubAMD <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p_pubAMD




# Tailing ####
dat_Tailing <- ARG_df.l %>%
  filter(dataCollection == "Tailings") %>% 
  filter(ARG %in% topARG_Tailing) %>%
  filter(relAbund < 0.5) # remove extreme values makes the graph look better

level.arg <- (dat_Tailing %>% group_by(ARG) %>% summarise(MRA = mean(relAbund)) %>% arrange(desc(MRA)))$ARG


dat_Tailing_sd <- 
  dat_Tailing %>% 
  group_by(ARG) %>%
  summarise(MRA=mean(relAbund), sd=sd(relAbund), resMech=unique(res_mechanism), drugType=unique(drug_type))

dat_Tailing$ARG <- factor(dat_Tailing$ARG,levels = level.arg)
p0_1 <-ggplot(dat_Tailing) + 
  geom_point(aes( x = ARG, y = relAbund),shape=1,alpha=0.6,color="darkgray") +
  theme_bw() +
  labs(y = "relative abundance", x = "", title = "") +
  scale_y_continuous(limits = c(-0.01,0.20),oob=rescale_none) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) +
  scale_x_discrete(labels=sapply(levels(dat_Tailing$ARG), function(x) ARGsubtype_abb_df$abb[which(ARGsubtype_abb_df$ARGsubtypes == x)] )) +
  theme(axis.text.x = element_text(angle = 90))

p0_1


# color the dots with top two most abundant ARG types and others
dat_Tailing_sd$drugType3 <- sapply(dat_Tailing_sd$drugType,
                                   function(x){
                                     if(!(x %in% c("multidrug","bacitracin"))) "other" else x
                                   })
drugType3.levels <- c("multidrug","bacitracin","other")
dat_Tailing_sd$drugType3 <- factor(dat_Tailing_sd$drugType3,levels = drugType3.levels)

Colors= sapply(drugType3.levels,function(x) drugType_color_df$Colors[drugType_color_df$DrugTypes==x])
Colors[3] <- "darkgray"


dat_Tailing_sd$ARG <- factor(dat_Tailing_sd$ARG, levels = level.arg)
p0_2 <- ggplot(dat_Tailing_sd) + 
  geom_pointrange(aes(data= ,x = ARG, y = MRA,ymin=MRA-sd, ymax=MRA+sd, color=drugType3),size=0.8) +
  scale_color_manual(values = Colors) + 
  scale_y_continuous(limits = c(-0.01,0.20),oob=rescale_none) +
  theme_half_open(11, rel_small = 1) + 
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  # rremove("legend") +
  rremove("y.axis") + rremove("y.text") + rremove("ylab")

p0_2




#library(ggpubr)
#library(cowplot)
aligned_plots <- align_plots(p0_1, p0_2+theme(legend.position = "none"), align="hv", axis="tblr")
p_Tailing <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p_Tailing



# AMDSedi ####
dat_AMDSedi <- ARG_df.l %>%
  filter(dataCollection == "AMDSedi") %>% 
  filter(ARG %in% topARG_AMDSedi)

dat_AMDSedi_sd <- 
  dat_AMDSedi %>% 
  group_by(ARG) %>% 
  summarise(MRA=mean(relAbund), sd=sd(relAbund), resMech=unique(res_mechanism), drugType=unique(drug_type))

dat_AMDSedi$ARG <- factor(dat_AMDSedi$ARG, levels = topARG_AMDSedi)
p0_1 <-ggplot(dat_AMDSedi) + 
  geom_point(aes( x = ARG, y = relAbund),shape=1,alpha=0.6,color="darkgray") +
  theme_bw() +
  labs(y = "relative abundance", x = "", title = "") +
  scale_y_continuous(limits = c(-0.01,0.20),oob=rescale_none) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) +
  scale_x_discrete(labels=sapply(levels(dat_AMDSedi$ARG), function(x) ARGsubtype_abb_df$abb[which(ARGsubtype_abb_df$ARGsubtypes == x)] )) +
  theme(axis.text.x = element_text(angle = 90))

p0_1

# 按照3种drug type给颜色
dat_AMDSedi_sd$drugType3 <- sapply(dat_AMDSedi_sd$drugType,
                                   function(x){
                                     if(!(x %in% c("multidrug","bacitracin"))) "other" else x
                                   })
#drugType3.levels <- c("multidrug","glycopeptide","other")
dat_AMDSedi_sd$drugType3 <- factor(dat_AMDSedi_sd$drugType3,levels = drugType3.levels)

Colors= sapply(drugType3.levels,function(x) drugType_color_df$Colors[drugType_color_df$DrugTypes==x])
Colors[3] <- "darkgray"


dat_AMDSedi_sd$ARG <- factor(dat_AMDSedi_sd$ARG, levels = topARG_AMDSedi)
p0_2 <- ggplot(dat_AMDSedi_sd) + 
  geom_pointrange(aes(data= ,x = ARG, y = MRA,ymin=MRA-sd, ymax=MRA+sd, color=drugType3),size=0.8) +
  scale_color_manual(values = Colors) + 
  scale_y_continuous(limits = c(-0.01,0.20),oob=rescale_none) +
  theme_half_open(11, rel_small = 1) + 
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  # rremove("legend") +
  rremove("y.axis") + rremove("y.text") + rremove("ylab")

p0_2




#library(ggpubr)
#library(cowplot)
aligned_plots <- align_plots(p0_1, p0_2+theme(legend.position = "none"), align="hv", axis="tblr")
p_AMDSedi <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p_AMDSedi



ggsave(ggarrange(p_pubAMD,p_Tailing,p_AMDSedi,nrow = 1,ncol = 3),
       device = "pdf", width = 12,height = 3.5,filename = "FigS5.a-c.AbundantARGsubtype.pdf")


# FigS5.d-f. Venn plots bars for abundant & frequent ARGs ---------------------------------------
ARGmapping_df <- ARG_df.l %>% select(ARG, drug_type, res_mechanism) %>% unique()

# pubAMD ####
FreqARG_pubAMD_df <- cbind.data.frame(ARG=FreqARG_pubAMD,
                                      FreqARG=T,
                                      stringsAsFactors=F)

AbundARG_pubAMD_df <- cbind.data.frame(ARG=topARG_pubAMD,
                                       topARG_pubAMD=T,
                                       stringsAsFactors=F)

pubAMD_df <- 
  merge(AbundARG_pubAMD_df, FreqARG_pubAMD_df, by="ARG",all = T) %>%
  mutate(ARGtype= sapply(ARG, function(x) ARGmapping_df$drug_type[which(ARGmapping_df$ARG == x )]))
pubAMD_df[is.na(pubAMD_df)] <- FALSE

pubAMD_df %>% filter(topARG_pubAMD & !FreqARG) %>% filter(ARGtype == "multidrug") %>% nrow()
pubAMD_df %>% filter(topARG_pubAMD & FreqARG) %>% filter(ARGtype == "multidrug") %>% nrow()
pubAMD_df %>% filter(!topARG_pubAMD & FreqARG) %>% filter(ARGtype == "multidrug") %>% nrow()


library(eulerr)
?euler
fit_pubAMD <- euler(pubAMD_df[-1], shape = "ellipse")
fit_pubAMD$original.values
plot(fit_pubAMD, fill = c("red", "steelblue4"), alpha=0.5,quantities=T,main = "pubAMD")


# Tailing ####
FreqARG_Tailing_df <- cbind.data.frame(ARG=FreqARG_Tailing,
                                       FreqARG=T,
                                       stringsAsFactors=F)

AbundARG_Tailing_df <- cbind.data.frame(ARG=topARG_Tailing,
                                        topARG_Tailing=T,
                                        stringsAsFactors=F)

Tailing_df <- 
  merge(AbundARG_Tailing_df, FreqARG_Tailing_df,  by="ARG",all = T) %>%
  mutate(ARGtype= sapply(ARG, function(x) ARGmapping_df$drug_type[which(ARGmapping_df$ARG == x )]))
Tailing_df[is.na(Tailing_df)] <- FALSE

fit_Tailing <- euler(Tailing_df[-1], shape = "ellipse")
plot(fit_Tailing, fill = c("red", "steelblue4"), alpha=0.5,quantities=T,main = "Tailing")


Tailing_df %>% filter(topARG_Tailing & !FreqARG) %>% filter(ARGtype == "multidrug") %>% nrow()
Tailing_df %>% filter(topARG_Tailing & FreqARG) %>% filter(ARGtype == "multidrug") %>% nrow()
Tailing_df %>% filter(!topARG_Tailing & FreqARG) %>% filter(ARGtype == "multidrug") %>% nrow()


# AMDSedi ####
FreqARG_AMDSedi_df <- cbind.data.frame(ARG=FreqARG_AMDSedi,
                                       FreqARG=T,
                                       stringsAsFactors=F)

AbundARG_AMDSedi_df <- cbind.data.frame(ARG=topARG_AMDSedi,
                                        topARG_AMDSedi=T,
                                        stringsAsFactors=F)

AMDSedi_df <- 
  merge(AbundARG_AMDSedi_df, FreqARG_AMDSedi_df,  by="ARG",all = T) %>%
  mutate(ARGtype= sapply(ARG, function(x) ARGmapping_df$drug_type[which(ARGmapping_df$ARG == x )]))
AMDSedi_df[is.na(AMDSedi_df)] <- FALSE

fit_AMDSedi <- euler(AMDSedi_df[-1], shape = "ellipse")
plot(fit_AMDSedi, fill = c("red", "steelblue4"), alpha=0.5,quantities=T,main = "AMDSedi")


AMDSedi_df %>% filter(topARG_AMDSedi & !FreqARG) %>% filter(ARGtype == "multidrug") %>% nrow()
AMDSedi_df %>% filter(topARG_AMDSedi & FreqARG) %>% filter(ARGtype == "multidrug") %>% nrow()
AMDSedi_df %>% filter(!topARG_AMDSedi & FreqARG) %>% filter(ARGtype == "multidrug") %>% nrow()



library(ggpubr)
ggsave(ggarrange(plot(fit_pubAMD, fill = c("orangered3", "steelblue3"), alpha=0.5,quantities=T,main = "pubAMD"),
                 plot(fit_Tailing, fill = c("orangered3", "steelblue3"), alpha=0.5,quantities=T,main = "Tailing"),
                 plot(fit_AMDSedi, fill = c("orangered3", "steelblue3"), alpha=0.5,quantities=T,main = "AMDSedi"),
                 ncol = 3, nrow = 1),
       width = 12, height = 3, device = "pdf", filename = "FigS5d-f.VennARGsubtype.pdf")



# numbers  ----------------
unique(c(topARG_AMDSedi, topARG_pubAMD, topARG_Tailing))

test <- bind_rows(
  pubAMD_df %>% filter(topARG_pubAMD) %>% select(ARG, ARGtype),
  Tailing_df  %>% filter(topARG_Tailing) %>% select(ARG, ARGtype),
  AMDSedi_df  %>% filter(topARG_AMDSedi) %>% select(ARG, ARGtype)
) %>% unique()


table(test$ARGtype)


# Table S9 -----------------------------------
ARG_inDB <- unique(c(FreqARG_AMDSedi,FreqARG_pubAMD,FreqARG_Tailing, topARG_AMDSedi, topARG_pubAMD, topARG_Tailing))
ARG_abb <- sapply(ARG_inDB, 
                  function(x) if(x %in% ARGsubtype_abb_df$ARGsubtypes) ARGsubtype_abb_df$abb[which(ARGsubtype_abb_df$ARGsubtypes == x)] else "")

DrugType <- sapply(ARG_inDB,
                   function(x) ARGmapping_df$drug_type[ARGmapping_df$ARG == x])

ResMech <- sapply(ARG_inDB,
                  function(x) ARGmapping_df$res_mechanism[ARGmapping_df$ARG == x])

AbundinPubAMD <- sapply(ARG_inDB,
                        function(x) if(x %in% topARG_pubAMD) "+" else "-")

AbundinAMDSedi <- sapply(ARG_inDB, function(x) if(x %in% topARG_AMDSedi) "+" else "-")

AbundinTailing <- sapply(ARG_inDB, function(x) if(x %in% topARG_Tailing) "+" else "-")


UbiquitPubAMD <- sapply(ARG_inDB, function(x) if(x %in% FreqARG_pubAMD) "+" else "-")

UbiquitAMDSedi <- sapply(ARG_inDB, function(x) if(x %in% FreqARG_AMDSedi) "+" else "-")

UbiquitTailing <- sapply(ARG_inDB, function(x) if(x %in% FreqARG_Tailing) "+" else "-")


TableS <- cbind.data.frame(ARG_inDB, ARG_abb, DrugType, ResMech, AbundinPubAMD, UbiquitPubAMD, AbundinTailing, UbiquitTailing, AbundinAMDSedi, UbiquitAMDSedi,
                           stringsAsFactors = F)

TableS <- TableS %>% arrange(ResMech)
write.csv(TableS, file = "TableS9.Abund.Freq.ARGsubtypes.csv", quote = F, row.names = F)
