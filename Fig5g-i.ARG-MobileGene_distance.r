setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")

library(dplyr)
library(ggplot2)
library(data.table)

plotDat.l <- fread("_data/nonRegARG-MobileGene_leastDistance.txt")

# plot ARG-MGEnzyme nearest distances based on mediam value  (in drug types)------------------------------------------
for(dc in unique(plotDat.l$sampleCollection)){
  #dc = levels(plotDat.l$sampleCollection)[3]
  
  plotDat_tmp <- plotDat.l %>% filter(sampleCollection == dc) %>% filter(core_label_1 != "unclassified")
  
  
  x_rank_df <- plotDat_tmp %>%  group_by(core_label_1) %>% summarise(meanDist = mean(distKb)) %>% arrange(meanDist)
  assign(paste("x_rank_df_",dc,sep = ""),x_rank_df,envir = .GlobalEnv)
  plotDat_tmp$core_label_1 <- factor(plotDat_tmp$core_label_1,levels = x_rank_df$core_label_1)
  
  
  
  P <- plotDat_tmp %>%
    ggplot(aes(x = core_label_1, y = distKb)) +
    stat_summary(fun.y = "mean", geom = "bar", position = "identity",fill="dodgerblue3", alpha = 0.5, width = 0.5) + #,color="black"
    #stat_summary(geom = "errorbar", fun.data = mean_sdl, position = position_dodge(0.9), width = 0.2) +
    geom_point(shape=1,color="darkgray",alpha=0.6,size=1) +
    #scale_y_log10() +
    theme_bw() + theme(axis.text.x = element_text (angle = 90, hjust = 1 )) +
    xlab("") + ylab("ARG-MGE enzymes nearest distance (kb)") +
    theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) # +
  #stat_compare_means(method = "anova", label.y = -10) +      # Add global p-value
  #geom_hline(yintercept = mean(plotDat_tmp$distKb), linetype = 2) + # Add horizontal line at base mean
  #stat_compare_means(label = "p.signif", method = "t.test",  ref.group = dt) + # Pairwise comparison against reference
  #annotate("text", x = i, y=-1, label = "ref")
  
  ?assign
  assign(x=paste("P_",dc,sep = ""), value = P, envir = .GlobalEnv)
}





library(ggpubr)
ggsave(ggarrange(P_pubAMD + coord_cartesian(ylim=c(0,150)),
                 P_Tailing + coord_cartesian(ylim=c(0,50)),
                 P_AMDSedi + coord_cartesian(ylim=c(0,80)),
                 nrow = 1,ncol = 3),
       filename = "Fig6g-i.Arg-MGEnzyme_nearestDistanceInDrugs.pdf", width = 9,height = 4)


# plot ARG-MGEnzyme nearest distances based on mediam value  (overall)------------------------------------------
for(dc in unique(plotDat.l$sampleCollection)){
  plotDat_tmp <-  plotDat.l %>% filter(sampleCollection == dc) %>% filter(core_label_1 != "unclassified")
  
  P_overall <- plotDat_tmp %>%
    ggplot(aes(x = 0, y = distKb)) +
    stat_summary(fun.y = "mean", geom = "bar", position = "identity",fill="dodgerblue3", alpha = 0.5, width = 0.5) + 
    #stat_summary(geom = "errorbar", fun.data = mean_sdl, position = position_dodge(0.9), width = 0.2) +
    geom_point(shape=1,color="darkgray",alpha=0.6,size=1) +
    #scale_y_log10() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) +
    # facet_wrap(vars(sampleCollection),nrow = 1,ncol = 3,scales = "free_y") +
    theme(axis.text.x = element_blank(), axis.ticks.x  = element_blank()) +  xlab("") + ylab("") +
    scale_x_discrete(labels= dc)
  
  assign(paste("P.overall_",dc,sep = ""), P_overall, envir = .GlobalEnv)
}
 

ggsave(ggarrange(P.overall_pubAMD + coord_cartesian(ylim=c(0,150)),
                 P.overall_Tailing + coord_cartesian(ylim=c(0,50)),
                 P.overall_AMDSedi + coord_cartesian(ylim=c(0,80)),
                 nrow = 1,ncol = 3),
       filename = "Fig6g-i.Arg-MGEnzyme_nearestDistanceOverall.pdf",  width = 3,height = 3.5)


# 计算不同drug type的再分组， pvalue等 --------------------------------------------------------
library(agricolae)

Dat <- plotDat.l %>% filter(core_label_1 != "unclassified") %>% 
  mutate(group= core_label_1) %>% select(sampleCollection, group, distKb)


dc=unique(Dat$sampleCollection)[2]

Dat_tmp <- Dat %>% filter(sampleCollection == dc)

model_anov <- aov(distKb~group,data=Dat_tmp)
com_anov = duncan.test(model_anov,"group",alpha=0.05,group=TRUE)

group1_df <- com_anov$group #group1是按com_anov给的分组
group1_df

Dat_tmp$group1 <- sapply(Dat_tmp$group, 
                         function(x) group1_df$groups[which(rownames(group1_df)==x)])

# group2是看图自己定义的分组
# group2_list <- list(grp1 = c("a","ab"),  grp2 = c("bc","c")) # dc == "pubAMD"
 group2_list <- list(grp1 = c("a"), grp2 = c("b","bc","bcd","cd"), grp3=c( "cde", "de","e")) #dc == "Tailing"
# group2_list <- list(grp1 = c("a"), grp2 = c("b"), grp3=c("bc"), grp4=c("cd","d")) #dc == "AMDSedi"



Dat_tmp$group2 <- sapply(Dat_tmp$group1,
                         function(x) names(which(sapply(group2_list, function(y) x %in% y)))) #根据group1找出对应的group2

#按group1分组，看组内不同drug type的对比显著性####
p.df_inGroup1 <- NULL
for(grp1 in unique(Dat_tmp$group1)){
  #  grp1=unique(Dat_tmp$group1)[2]
  Dat_tmp_grp <- Dat_tmp %>% filter(group1 == grp1)
  
  if(length(unique(Dat_tmp_grp$group)) ==1 ) {Anova.p = NA; Kruskal.p = NA; next}
  
  model_tmp <- aov(distKb~group,data=Dat_tmp_grp)
  summary(model_tmp)
  Anova.p = summary(model_tmp)[[1]][["Pr(>F)"]][1]
  
  row_tmp <- cbind.data.frame(group1=grp1,
                              anova.p = Anova.p,
                              stringsAsFactors=F)
  
  p.df_inGroup1 <- rbind.data.frame(row_tmp,
                                    p.df_inGroup1,
                                    stringsAsFactors = F)
  
}
p.df_inGroup1


#按group2分组，看组内不同drug type的对比显著性(如果组内显著性都不高，p>0.05，就按照group2来画图) ####
p.df_inGroup2 <- NULL
for(grp2 in unique(Dat_tmp$group2)){
  #  grp2=unique(Dat_tmp$group2)[2]
  Dat_tmp_grp <- Dat_tmp %>% filter(group2 == grp2)
  
  if(length(unique(Dat_tmp_grp$group)) ==1 ) {Anova.p = NA; Kruskal.p = NA; next}
  
  model_tmp <- aov(distKb~group,data=Dat_tmp_grp)
  summary(model_tmp)
  Anova.p = summary(model_tmp)[[1]][["Pr(>F)"]][1]

  
  row_tmp <- cbind.data.frame(group2=grp2,
                              anova.p = Anova.p,
                              stringsAsFactors=F)
  
  p.df_inGroup2 <- rbind.data.frame(row_tmp,
                                    p.df_inGroup2,
                                    stringsAsFactors = F)
  
}
p.df_inGroup2
ARGType_groupDf <- Dat_tmp %>% select(group,group1,group2) %>% unique() %>% arrange(desc(group1))
ARGType_groupDf


# 先看看group3组内p值是否>0.05,如果>0.05,则计算group3的两两组间pvalue
group_selected = "group2"
Dat_tmp <- as.data.frame(Dat_tmp)
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
grpPair_Df$Wilcox.sig <- cut(grpPair_Df$Wilcox.p, breaks = c(-Inf,0.001,0.05,0.01,Inf), labels = c("***","**","*",NA))
grpPair_Df

assign(paste("p.df_inGroup1_",dc,sep = ""), p.df_inGroup1, envir = .GlobalEnv)
assign(paste("p.df_inGroup2_",dc,sep = ""), p.df_inGroup2, envir = .GlobalEnv)
assign(paste("grpPair_Df_",dc,sep = ""), grpPair_Df, envir = .GlobalEnv)

assign( paste("ARGType_groupDf_",dc,sep = ""),
        ARGType_groupDf,
        envir = .GlobalEnv )
