setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")
library(dplyr)
library(ggplot2)
library(data.table)

leastDistanceDf <- fread("_data/nonRegARG-MRG_coOccurrence_inBins.csv",fill=TRUE)

plotDat <- leastDistanceDf %>% 
  mutate(binFileName=sapply(strsplit(ARG_ORF,"&",fixed = T),"[[",2)) %>%
  group_by(binFileName, core_label_1) %>% 
  summarise(sampleCollection = unique(sampleCollection),
            n=n(),
            avg_leastDistance = mean(leastDistance), 
            mediam_leastDistance = median(leastDistance)) %>%
  as.data.frame()


plotDat.l <-  plotDat %>%  
  select(-n) %>%
  reshape2::melt(id.vars=c("binFileName","core_label_1","sampleCollection"), variable.name="type_distance") %>%
  mutate(distKb = value/1000)



# plot ARG-MetalRG nearest distances based on average values  (in drug types)------------------------------------------
for(dc in unique(plotDat.l$sampleCollection)){
  # dc =  unique(plotDat_form2$sampleCollection)[3]
  
  plotDat_tmp <- plotDat.l %>% filter(sampleCollection == dc)  %>% 
    filter(type_distance == "avg_leastDistance") %>%
    filter(core_label_1 != "unclassified") 
  
  x_rank_df <- plotDat_tmp %>% 
    group_by(core_label_1) %>%
    summarise(meanDist = mean(distKb)) %>% 
    arrange(meanDist)
  assign(paste("x_rank_df_",dc,sep = ""),x_rank_df,envir = .GlobalEnv)
  plotDat_tmp$core_label_1 <- factor(plotDat_tmp$core_label_1,levels = x_rank_df$core_label_1)
  
  
  P <- plotDat_tmp %>%
    ggplot(aes(x = core_label_1, y = distKb)) +
    stat_summary(fun.y = "mean", geom = "bar", position = "identity",fill="dodgerblue3", alpha = 0.5, width = 0.5) + #,color="black"
    #stat_summary(geom = "errorbar", fun.data = mean_sdl, position = position_dodge(0.9), width = 0.2) +
    geom_point(shape=1,color="darkgray",alpha=0.6,size=1) +
    #scale_y_log10() +
    theme_bw() + theme(axis.text.x = element_text (angle = 90, hjust = 1 )) +
    xlab("") + ylab("ARG-MRG nearest distance (kb)") +
    theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) +
    coord_cartesian(ylim=c(-10, 350))
  
  P
  
  
  ?assign
  assign(x=paste("P_",dc,sep = ""), value = P, envir = .GlobalEnv)
}

P_Tailing
P_AMDSedi



# plot ARG-MetalRG nearest distances based on avg value  (overall)------------------------------------------
plotDat_tmp <- plotDat.l %>% filter(type_distance == "avg_leastDistance") %>%
  filter(core_label_1 != "unclassified") #%>% filter(distKb <= 350)
P_overall <- plotDat_tmp %>%
  ggplot(aes(x = 0, y = distKb)) +
  stat_summary(fun.y = "mean", geom = "bar", position = "identity",fill="dodgerblue3", alpha = 0.5, width = 0.5) + 
  #stat_summary(geom = "errorbar", fun.data = mean_sdl, position = position_dodge(0.9), width = 0.2) +
  geom_point(shape=1,color="darkgray",alpha=0.6,size=1) +
  #scale_y_log10() +
  theme_bw() + theme(axis.text.x = element_text (angle = 45, hjust = 1 )) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) +
  facet_wrap(vars(sampleCollection),nrow = 1,ncol = 3) + coord_cartesian(ylim=c(-10, 350)) +
  theme(axis.text.x = element_blank(), axis.ticks.x  = element_blank()) +  xlab("") + ylab("") +
  scale_x_discrete(labels="Overall")+
  coord_cartesian(ylim=c(-10, 350))

P_overall


library(ggpubr)
ggsave(ggarrange(P_Dump,P_Tailing,P_AMDSedi,nrow = 1,ncol = 3),
       filename = "FigS9.a-b.Arg-Mrg_nearestDistanceInDrugs_bins.pdf", width = 9,height = 4)
ggsave(P_overall,
       filename = "FigS9.a-b.Arg-Mrg_nearestDistanceOverall_bins.pdf", width = 3,height = 3.5)


# calculate group-group pvalue --------------------------------------------------------
library(agricolae)
Dat <- plotDat.l %>% filter(type_distance == "avg_leastDistance") %>%
  filter(core_label_1 != "unclassified") %>% # filter(distKb <= 350) %>%
  mutate(group = core_label_1) %>% select(sampleCollection, group, distKb)



dc=unique(Dat$sampleCollection)[3]

Dat_tmp <- Dat %>% filter(sampleCollection == dc)

# divide drug types into groups based on least distance to MRG
model_anov <- aov(distKb~group,data=Dat_tmp)

#  duncan.test(model,"group",alpha=0.05,group=TRUE)
com_anov = duncan.test(model_anov,"group",alpha=0.05,group=TRUE)

# model_kruskal <- kruskal.test(distKb~group,data=Dat_tmp)
# com_kruskal = duncan.test(model_kruskal,"group",alpha=0.05,group=TRUE)


group1_df <- com_anov$group

Dat_tmp$group1 <- sapply(Dat_tmp$group, 
                         function(x) group1_df$groups[which(rownames(group1_df)==x)])

in.group1_p.df <- NULL
for(grp1 in unique(Dat_tmp$group1)){
  #  grp1=unique(Dat_tmp$group1)[2]
  Dat_tmp_grp <- Dat_tmp %>% filter(group1 == grp1)
  
  if(length(unique(Dat_tmp_grp$group)) ==1 ) {Anova.p = NA; Kruskal.p = NA; next}
  
  model_tmp <- aov(distKb~group,data=Dat_tmp_grp)
  summary(model_tmp)
  Anova.p = summary(model_tmp)[[1]][["Pr(>F)"]][1]
  
  model_tmp <- kruskal.test(distKb~group,data=Dat_tmp_grp)
  #summary(model_tmp)
  model_tmp$p.value
  Kruskal.p = model_tmp$p.value
  
  
  row_tmp <- cbind.data.frame(group1=grp1,
                              anova.p = Anova.p,
                              kruskal.p = Kruskal.p,
                              stringsAsFactors=F)
  
  in.group1_p.df <- rbind.data.frame(row_tmp,
                                     in.group1_p.df,
                                     stringsAsFactors = F)
  
}
in.group1_p.df

group1_df

#manually define for each data collection
#group2_list <- list(grp1 = c("a"), grp2 = c("b"), grp3=c("bc","cde","cd"),grp4=c("de","e"))  # AMDSedi <200kb
group2_list <- list(grp1 = c("c","bc"), grp2 = c("b"),grp3=c("ab","a"), grp4=c(""))  # tailings <200kb


Dat_tmp$group2 <- sapply(Dat_tmp$group1,
                         function(x) names(which(sapply(group2_list, function(y) x %in% y)))) #根据group1找出对应的group2


in.group2_p.df <- NULL
for(grp2 in unique(Dat_tmp$group2)){
  #  grp2=unique(Dat_tmp$group2)[2]
  Dat_tmp_grp <- Dat_tmp %>% filter(group2 == grp2)
  
  if(length(unique(Dat_tmp_grp$group)) ==1 ) {Anova.p = NA; Kruskal.p = NA; next}
  
  model_tmp <- aov(distKb~group,data=Dat_tmp_grp)
  summary(model_tmp)
  Anova.p = summary(model_tmp)[[1]][["Pr(>F)"]][1]
  
  model_tmp <- kruskal.test(distKb~group,data=Dat_tmp_grp)
  #summary(model_tmp)
  model_tmp$p.value
  Kruskal.p = model_tmp$p.value
  
  
  row_tmp <- cbind.data.frame(group2=grp2,
                              anova.p = Anova.p,
                              kruskal.p = Kruskal.p,
                              stringsAsFactors=F)
  
  in.group2_p.df <- rbind.data.frame(row_tmp,
                                     in.group2_p.df,
                                     stringsAsFactors = F)
  
}
in.group2_p.df


ARGtype.grpInfo <- Dat_tmp %>% select(group,group1,group2) %>% unique() %>% arrange(desc(group1))
ARGtype.grpInfo

group_selected = "group2"
grpPair_Df <- data.frame(t(combn(unique(Dat_tmp[,group_selected]), 2)))
grpPair_Df$X1 <- as.character(grpPair_Df$X1)
grpPair_Df$X2 <- as.character(grpPair_Df$X2)
grpPair_Df$Ttest.p <- NA
grpPair_Df$Wilcox.p <- NA
for(i in c(1:nrow(grpPair_Df))){
  #i=1
  X1 = grpPair_Df$X1[i]
  X2 = grpPair_Df$X2[i]
  
  X1_value <- Dat_tmp$distKb[which(Dat_tmp[,group_selected] == X1)]
  X2_value <- Dat_tmp$distKb[which(Dat_tmp[,group_selected] == X2)]
  
  grpPair_Df$Ttest.p[i] <-  t.test(X1_value,X2_value)$p.value
  grpPair_Df$Wilcox.p[i] <-  wilcox.test(X1_value,X2_value)$p.value
}

grpPair_Df$Ttest.sig <- cut(grpPair_Df$Ttest.p, breaks = c(-Inf,0.001,0.05,0.01,Inf), labels = c("***","**","*",NA))
grpPair_Df

assign(paste("p.inGroup1_df_",dc,sep = ""), in.group1_p.df, envir = .GlobalEnv)
assign(paste("p.inGroup2_df_",dc,sep = ""), in.group2_p.df, envir = .GlobalEnv)
assign(paste("p.grpPair_df_",dc,sep = ""), grpPair_Df, envir = .GlobalEnv)

assign( paste("ARGtype_groupInfo_",dc,sep = ""),
        ARGtype.grpInfo,
        envir = .GlobalEnv )


# plot ARG-MRG type pairs in MAGs- -------------------------------
rm(list = ls())
ARG.MRG_pairs_df <- fread("_data/nonRegARG-MRG_pairs_inBins.csv")

load("_data/otherARGtypes.RData")
load("_data/otherMetalRGtypes.RData")

x_rank_df_Tailing <- c("sulfonamide","multidrug","tetracycline","other","phenicol","glycopeptide","bacitracin","MLS",
                       "beta-lactam","aminoglycoside","mupirocin")
x_rank_df_AMDSedi <- c("sulfonamide","multidrug","phenicol","tetracycline","other","aminoglycoside","MLS","bacitracin",
                       "glycopeptide","beta-lactam","mupirocin")
x_rank_df_AMDSedi[!(x_rank_df_AMDSedi %in% ARG.MRG_pairs_df$ARG_type)]
x_rank_df_Tailing[!(x_rank_df_Tailing %in% ARG.MRG_pairs_df$ARG_type)]


library(ggsci)
library(scales)
ARG.MRG_pairs_df <- ARG.MRG_pairs_df %>% filter(NearestDistance < 100000)
ARG.MRG_pairs_df$ARG_type[which(ARG.MRG_pairs_df$ARG_type %in% otherARGtypes)] <- "other"
ARG.MRG_pairs_df$MRG_type[which(ARG.MRG_pairs_df$MRG_type %in% otherMetals)] <- "other_metal_resistance"


dat_pairs <- ARG.MRG_pairs_df %>% 
  group_by(spCollection, ARG_type, MRG_type) %>% 
  summarise(n=n()) %>% 
  as.data.frame() %>% 
  merge((ARG.MRG_pairs_df %>% group_by(spCollection,ARG_type) %>% summarise(sum=n())),
        by=c("spCollection","ARG_type")) %>%
  mutate(freq_inDrug = n/sum)

MRGtype_rank <- c("Arsenic_resistance", "Chromium_resistance", "Cobalt_resistance","Copper_resistance",
                  "Iron_resistance","Lead_resistance","Mercury_resistance","Multi-metal_resistance",
                  "Nickel_resistance","Zinc_resistance","other_metal_resistance") #从_colors.RData拷过来的

dat_pairs$spCollection <- factor(dat_pairs$spCollection,levels = c("Dump","Tailing","AMDSedi"))
dat_pairs$MRG_type <- factor(dat_pairs$MRG_type, levels = MRGtype_rank)


MetalType_colorsDf2 <- cbind.data.frame(
  Colors=c("#EFC000FF","#4DBBD5FF","#00A087FF","steelblue3","#F39B7FFF","#8491B4FF","#91D1C2FF","#E64B35FF" ,
           "#7E6148FF","#B09C85FF","gray"),
  DrugTypes=c("Arsenic_resistance", "Chromium_resistance", "Cobalt_resistance","Copper_resistance",
              "Iron_resistance","Lead_resistance","Mercury_resistance","Multi-metal_resistance",
              "Nickel_resistance","Zinc_resistance","other_metal_resistance"),
  stringsAsFactors=F
)

Colors = sapply(MRGtype_rank, function(x)MetalType_colorsDf2$Colors[which(MetalType_colorsDf2$DrugTypes == x)])


library(ggplot2)

for(dc in levels(dat_pairs$spCollection)[2:3]){
  #  dc = levels(dat_pairs$spCollection)[2]
  x_rank_df <- eval(parse(text = paste("x_rank_df_", dc, sep = "")))
  
  dat_pairs_tmp <- dat_pairs %>% filter(spCollection == dc) %>% filter(ARG_type != "unclassified")
  dat_pairs_tmp$ARG_type <- factor(dat_pairs_tmp$ARG_type,
                                   levels = x_rank_df)
  P <- ggplot(dat_pairs_tmp,aes(x=ARG_type, y=freq_inDrug, fill=MRG_type)) +
    geom_col() + 
    scale_fill_manual(values = Colors) +
    theme_bw() +
    theme (axis.text.x = element_text (angle = 90, hjust = 1 ) ) + xlab("") +ylab("Nearest MRG composition") +
    #guides(fill = guide_legend(nrow = 1)) +   
    #theme(legend.position="none") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  P
  assign(paste("P_",dc,sep = ""),P, envir = .GlobalEnv)
}

P_Tailing

library(ggpubr)
ggsave(ggarrange(P_Tailing+theme(legend.position="none") ,P_Tailing+theme(legend.position="none") ,P_AMDSedi+theme(legend.position="none") ,
                 nrow = 1,ncol = 3),
       filename = "FigS9.c-d.MRGcompP_drugs.pdf", width = 9,height = 4) 



# plot in overall -------


dat_pairs_overall <- ARG.MRG_pairs_df %>% 
  group_by(spCollection, MRG_type) %>% 
  summarise(n=n()) %>% as.data.frame() %>% 
  merge((ARG.MRG_pairs_df %>% group_by(spCollection) %>% summarise(sum=n())),
        by=c("spCollection")) %>%
  mutate(freq=n/sum)

dat_pairs_overall$spCollection <- factor(dat_pairs_overall$spCollection,levels = c("Dump","AMDSedi","Tailing"))

dat_pairs_overall$MRG_type <- factor(dat_pairs_overall$MRG_type,
                                     levels = MRGtype_rank)


MRGcompP_overall <- ggplot(dat_pairs_overall,aes(x=1, y=freq,fill =MRG_type)) +
  geom_col() + 
  scale_fill_manual(values = Colors) +
  #guides(fill = guide_legend(nrow = 1)) +   
  theme(legend.position="none") +
  facet_wrap(vars(spCollection),ncol = 3,nrow = 1) +
  theme ( axis.text.x = element_blank (), axis.text.y = element_blank() ) + xlab("") +ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

MRGcompP_overall

ggsave(MRGcompP_overall+theme(legend.position="none"),
       filename = "FigS9.c-d.Arg-MrgComposition-overall.pdf",
       device = "pdf", width = 3,height = 3.5)
