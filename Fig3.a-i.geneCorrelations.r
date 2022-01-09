library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer) 
library(data.table)

# read data -----------------------------------------
allARG_df.l <- fread("_data/nonRegARG.abund.csv")
mobileGene_df <- fread("_data/MobilityGene.abund.csv")
B124_df <- fread("_data/B124.abund.csv")
crAssPhage_df <- fread("_data/crAssPhage.abund.csv")
metalRG_df <- fread("_data/metalRG.abund.csv")

# calculate total ARG  #####
ARG_df <- allARG_df.l %>% 
  dplyr::group_by(sampleID) %>%
  dplyr::summarise(totalARG=sum(DepthPG),sampleTypes=unique(sampleTypes)) %>%
  as.data.frame(stringsAsFactors=F) %>%
  dplyr::select(sampleID,totalARG)


# mobility genes #####
IntTnp_depth_wideDf <- 
  mobileGene_df %>% 
  reshape2::dcast(sampleID ~ target,value.var =  "DepthPG")


# metal resistance genes #####
metals_df <- metalRG_df %>%
  group_by(sampleID) %>% summarise(Metals = sum(DepthPG))

# merge all the data ------------------------------------
dat <- merge(merge(merge(ARG_df, metals_df, by="sampleID"), 
                   merge(B124_df, crAssPhage_df, by="sampleID"),
                   by="sampleID"),
             IntTnp_depth_wideDf, by="sampleID")

dat$sampleTypes <- sapply(dat$sampleID,
                          function(x) allARG_df.l$sampleTypes[which(allARG_df.l$sampleID == x)[1]])
dataCollection <- vector("character", length = nrow(dat))
for(i in c(1:nrow(dat))){
  if(dat$sampleTypes[i] == "Tailing"){
    dataCollection[i] <- "Tailing"
  }else if(grepl("^LD", dat$sampleID[i], perl = T)){
    dataCollection[i] <- "AMDSedi"
  }else{
    dataCollection[i] <- "pubAMD"
  }
  
}
dat$dataCollection <- dataCollection



# plotting -----------------------------------------
lm_rp_plain <- function(df, y_column, x_column){
  y = y_column
  x = x_column
  m <- lm(y ~ x, df);
  corrTest=cor.test(y, x);
  a = round(unname(coef(m)[1]), digits = 2)
  b = round(unname(coef(m)[2]), digits = 2)
  R = round(corrTest$estimate,digits = 3)
  pvalue = ifelse(corrTest$p.value == 0, 
                  "< 2.2e-16", 
                  paste("= ",formatC(corrTest$p.value, format = "e", digits = 3),sep = ""))
  
  return(paste("y = ",a," + ",b,"x, r = ",R,", p ", pvalue, sep = ""))
}

dat$dataCollection <- factor(dat$dataCollection,levels = c("pubAMD","Tailing","AMDSedi"))
dat$sampleTypes <- factor(dat$sampleTypes, levels = c("AMD","AMDSediment","biofilm","Tailing"))



# equations ----------------
plotDat <- dat
equationByDatacoll_log_df<-NULL
for(dc in unique(plotDat$dataCollection)){
  # dc=unique(plotDat$dataCollection)[1]
  
  plotDat_dc <- plotDat %>% filter(dataCollection == dc)
  
  
  
  for(yVar in colnames(plotDat_dc)[3:11]){
    
    #yVar = colnames(plotDat_dc)[3]
    plotDat_tmp <- plotDat_dc
    colnames(plotDat_tmp)[which(colnames(plotDat_tmp) == yVar)] <- "yVar"
    
    plotDat_tmp$yVar <- log1p(plotDat_tmp$yVar) #把y进行log转换
    plotDat_tmp$totalARG <- log1p(plotDat_tmp$totalARG)  #把x进行log转换
    eqt_lm <- lm_rp_plain(df=plotDat_tmp, plotDat_tmp$yVar, plotDat_tmp$totalARG)
    
    equationDf_c <- cbind.data.frame(dataCollection = dc, 
                                     y = yVar,
                                     eqation = eqt_lm,
                                     stringsAsFactors=F)
    
    equationByDatacoll_log_df <- bind_rows(equationByDatacoll_log_df, equationDf_c)
  }
  
}

equationByDatacoll_log_df$r = sub("^.*r\\s\\=\\s(.*)\\,.*$","\\1",equationByDatacoll_log_df$eqation)



# equations------------------------

plotDat <- dat 
equationByDatacoll_df<-NULL
for(dc in unique(plotDat$dataCollection)){
  #dc=unique(plotDat$dataCollection)[1]
  
  plotDat_dc <- plotDat %>% filter(dataCollection == dc)
  
  
  
  for(yVar in colnames(plotDat_dc)[3:11]){
    
    #yVar = colnames(plotDat_dc)[3]
    plotDat_tmp <- plotDat_dc
    colnames(plotDat_tmp)[which(colnames(plotDat_tmp) == yVar)] <- "yVar"
    
    #plotDat_tmp$yVar <- log1p(plotDat_tmp$yVar) #把y进行log转换
    #plotDat_tmp$totalARG <- log1p(plotDat_tmp$totalARG)  #把x进行log转换
    eqt_lm <- lm_rp_plain(df=plotDat_tmp, plotDat_tmp$yVar, plotDat_tmp$totalARG)
    
    equationDf_c <- cbind.data.frame(dataCollection = dc, 
                                     y = yVar,
                                     eqation = eqt_lm,
                                     stringsAsFactors=F)
    
    equationByDatacoll_df <- bind_rows(equationByDatacoll_df, equationDf_c)
  }
  
}

equationByDatacoll_df$r = sub("^.*r\\s\\=\\s(.*)\\,.*$","\\1",equationByDatacoll_df$eqation)

test<- merge(equationByDatacoll_df, equationByDatacoll_log_df, by=c("dataCollection","y"))


# plotting Fig3a-c.---------------

log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}


# plot for metalRG and class 1 integron integrase
library(ggh4x)
datColl_colorDf <- cbind.data.frame(dataColl=c("pubAMD","AMDSedi","Tailing"),
                                    Color=c("#6A6599FF", "#EFC000FF", "#CD534CFF"),
                                    stringsAsFactors=F)


Colors = sapply(levels(dat$dataCollection),
                function(x) datColl_colorDf$Color[which(datColl_colorDf$dataColl == x)])
  
for(yVar in c("Metals","class1Int")){  # colnames(dat)[3:(ncol(dat)-2)]
  plotDat <- dat 
  colnames(plotDat)[which(colnames(plotDat) == yVar)] <- "yVar"
  
  p <- ggplot(plotDat,aes(x=totalARG, y=yVar)) + 
    geom_point(aes(fill = dataCollection), shape=21, color = "black",size=3, show.legend = FALSE)+
    scale_fill_manual(values = Colors) +
    geom_smooth(method = "lm",na.rm = T,aes(color = dataCollection),se = F) +
    scale_color_manual(values = Colors) +
    facet_wrap(vars(dataCollection), 
               scales = "free", 
               nrow = 1,ncol = 3) +
    xlab("Total ARG Abundance (coverage ×/Gb)") + ylab(strsplit(yVar,"_",fixed = T)[[1]][1])  +
    theme_bw() + theme(panel.grid = element_blank())
  p
  
  ggsave(p + 
           scale_x_log10(guide = "axis_minor", minor_breaks=log10_minor_break()) +
           scale_y_log10(guide = "axis_minor", minor_breaks=log10_minor_break()), 
         filename =  paste("Fig3.scatterplots_ncol3_3colors_",yVar,"_scalelog.pdf",sep = ""), 
         device = "pdf", width = 8.5, height = 2.5)
 }


# plot for B124 and crAssPhage separately ####
plotDat <- dat %>% 
  select(sampleID, totalARG, B124_DepthPG, crAssPhage_DepthPG, sampleTypes, dataCollection) %>%
  reshape2::melt(id.vars=c("sampleID","totalARG","sampleTypes","dataCollection"),value.name = "phage")


yVar = "phage"
colnames(plotDat)[which(colnames(plotDat) == yVar)] <- "yVar"


plotDat$variable

p<-ggplot() + 
  geom_point(data= plotDat %>% filter(variable == "crAssPhage_DepthPG"), shape=0, # empty square
             aes(color = dataCollection, x=totalARG, y=yVar),
             size=3) +
  geom_point(data= plotDat %>% filter(variable == "B124_DepthPG"), shape=2, # empty triangle
             aes(color = dataCollection, x=totalARG, y=yVar),
             size=3, show.legend = FALSE) +
  scale_fill_manual(values = Colors) +
  scale_fill_manual(values = Colors) +
  scale_shape_manual(values = c(1,21)) +
  scale_color_manual(values = Colors) +
  facet_wrap(vars(dataCollection), scales = "free", nrow = 1,ncol = 3) +
  xlab("TotalARG abundance (coverage ×/Gb)") + ylab(strsplit(yVar,"_",fixed = T)[[1]][1])  +
  theme_bw() + theme(panel.grid = element_blank()) 


p
ggsave(p + 
         scale_x_log10(guide = "axis_minor", minor_breaks=log10_minor_break()) + 
         scale_y_log10(guide = "axis_minor", minor_breaks=log10_minor_break()), 
       filename =  paste("Fig3d-f.",yVar,"_scalelog.pdf",sep = ""), 
       device = "pdf", width = 8.5, height = 2.5)



# plotting Fig6d-f ---------------

for(yVar in c("transposase")){  # colnames(dat)[3:(ncol(dat)-2)]
  plotDat <- dat 
  colnames(plotDat)[which(colnames(plotDat) == yVar)] <- "yVar"
  
  p <- ggplot(plotDat,aes(x=totalARG, y=yVar)) + 
    geom_point(aes(fill = dataCollection), shape=21, color = "black",size=3, show.legend = FALSE)+
    scale_fill_manual(values = Colors) +
    geom_smooth(method = "lm",na.rm = T,aes(color = dataCollection),se = F) +
    scale_color_manual(values = Colors) +
    facet_wrap(vars(dataCollection), 
               scales = "free", 
               nrow = 1,ncol = 3) +
    xlab("Total ARG Abundance (coverage ×/Gb)") + ylab(strsplit(yVar,"_",fixed = T)[[1]][1])  +
    theme_bw() + theme(panel.grid = element_blank())
  p
  
  ggsave(p + 
           scale_x_log10(guide = "axis_minor", minor_breaks=log10_minor_break()) +
           scale_y_log10(guide = "axis_minor", minor_breaks=log10_minor_break()), 
         filename =  paste("Fig6d-f.scatterplots_ncol3_3colors_",yVar,"_scalelog.pdf",sep = ""), 
         device = "pdf", width = 8.5, height = 2.5)
}

# plotting FigS13 ---------------

for(yVar in c("integrase", "recombinase", "resolvase")){  # colnames(dat)[3:(ncol(dat)-2)]
  plotDat <- dat 
  colnames(plotDat)[which(colnames(plotDat) == yVar)] <- "yVar"
  
  p <- ggplot(plotDat,aes(x=totalARG, y=yVar)) + 
    geom_point(aes(fill = dataCollection), shape=21, color = "black",size=3, show.legend = FALSE)+
    scale_fill_manual(values = Colors) +
    geom_smooth(method = "lm",na.rm = T,aes(color = dataCollection),se = F) +
    scale_color_manual(values = Colors) +
    facet_wrap(vars(dataCollection), 
               scales = "free", 
               nrow = 1,ncol = 3) +
    xlab("Total ARG Abundance (coverage ×/Gb)") + ylab(strsplit(yVar,"_",fixed = T)[[1]][1])  +
    theme_bw() + theme(panel.grid = element_blank())
  p
  
  ggsave(p, # + 
          # scale_x_log10(guide = "axis_minor", minor_breaks=log10_minor_break()) +
          # scale_y_log10(guide = "axis_minor", minor_breaks=log10_minor_break()), 
         filename =  paste("FigS13.plots/scatterplots_ncol3_3colors_",yVar,"_scalelog.pdf",sep = ""), 
         device = "pdf", width = 8.5, height = 2.5)
}
