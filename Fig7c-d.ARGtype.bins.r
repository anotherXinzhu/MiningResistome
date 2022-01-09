rm(list = ls())
library(data.table)
library(ggplot2)
library(ggsci)
library(scales)
library(dplyr)


ARG.bins <- fread("_data/nonRegARG.bins.txt")
load("_data/otherARGtypes.RData")

for(ecosystem in c("AMDSedi","Tailing")){
  
  # ecosystem = "AMDSedi"
  dat <- ARG.bins %>% filter(grepl(ecosystem, binFileName))
  colnames(dat)
  
  
  taxonLevel = "Phylum"
  
  dat_ArgComposition <- dat; colnames(dat_ArgComposition)[which(colnames(dat_ArgComposition) == taxonLevel)] <- "taxon"
  corrected_taxon <- vector("character",nrow(dat_ArgComposition))
  
  for(i in c(1:nrow(dat_ArgComposition))){
    #i=3
    taxon <- dat_ArgComposition$taxon[i]
    
    if(taxon=="") next
    parts <- strsplit(taxon," ",fixed = T)[[1]]
    
    for(pt in parts){
      j=which(parts == pt)
      
      parts[j] <- sub("_[A-Z]+","",pt)
      
    }
    
    corrected_taxon[i] <- paste(parts,collapse = " ")
  }
  dat_ArgComposition$taxon <- corrected_taxon
  
  # calculate top taxons ----------------------------------------------
  Taxon_rankDf <- merge(merge(dat_ArgComposition %>% select(taxon, binFileName) %>% unique() %>% group_by(taxon) %>% summarise(numBins = n()),
                              dat_ArgComposition %>% select(taxon, ARG) %>% unique() %>% group_by(taxon) %>% summarise(numArgSubtype = n()),
                              by="taxon"),
                        dat_ArgComposition %>% select(taxon, drugType) %>% unique() %>% group_by(taxon) %>% summarise(numArgDrugType = n()),
                        by="taxon") %>% filter(taxon != "")  
  
  topTaxonBy <- "numBins" 
  Taxon_rankDf <- Taxon_rankDf %>% arrange(across(all_of(topTaxonBy),desc))
  topTaxons <- Taxon_rankDf$taxon[1:10]
  
  # define other ARG types -------------------------------------------
  MinorArgType_df <- dat_ArgComposition %>% 
    filter(taxon != "") %>% # remove unidentified taxonomic members
    filter(taxon %in% topTaxons) %>% 
    group_by(taxon, drugType) %>% summarise(count=n()) %>% mutate(freq = count/sum(count)) %>% 
    group_by(drugType) %>% summarise(avgRelAbund= mean(freq), maxRelAbund=max(freq)) %>% 
    arrange(desc(maxRelAbund)) #%>% filter(avgRelAbund <= 0.02)
  
  #MinorDrugTypes_list[[taxonLevel]] <- MinorArgType_df$drugType[which(MinorArgType_df$maxRelAbund < 0.02)] #不好对N30和N50进行循环。不要了。
  
  #MinorArgTypes <- MinorArgType_df$drugType[which(MinorArgType_df$maxRelAbund < 0.02)] #按照每个图不同情况来定义other drug types
  #dat_ArgComposition$drugType[which(dat_ArgComposition$drugType %in% MinorArgTypes)]  <- "other"
  
  dat_ArgComposition$drugType[which(dat_ArgComposition$drugType %in% otherARGtypes)] <- "other" 
  
  
  plotDat <- dat_ArgComposition %>% filter(taxon %in% topTaxons) %>% filter(drugType != "unclassified") %>%
    group_by(taxon, drugType) %>% summarise(count=n()) %>% mutate(freq = count/sum(count)) %>% 
    mutate(taxon_fct = factor(taxon, levels = topTaxons)) 
  
  plotDat$drugType <- factor(plotDat$drugType, 
                             levels = c("aminoglycoside", "bacitracin", "beta-lactam", "glycopeptide", "multidrug", "MLS", "mupirocin",
                                        "phenicol","sulfonamide" , "tetracycline", "other" )) #是为了把others放在最后，multidrug排到红色。
  
  library(ggsci)
  drugType_colorsDf <- cbind.data.frame(
    Colors=c(pal_npg("nrc")(10), "gray","gold3","orchid"),
    DrugTypes=c("aminoglycoside", "bacitracin", "beta-lactam","glycopeptide",
                "MLS","tetracycline","phenicol","multidrug",
                "mupirocin","sulfonamide","other","fosmidomycin","rifamycin"), # 后来添加了fosmidomycin 和 rifamycin是SARG
    stringsAsFactors=F
  )
  colors <- sapply(levels(plotDat$drugType), 
                   function(x) drugType_colorsDf$Colors[which(drugType_colorsDf$DrugTypes== x)])
  
 
  
  #横着的bar
  plotDat$drugType <- factor(plotDat$drugType, 
                             levels = c("other", "tetracycline","sulfonamide","phenicol",  
                                        "mupirocin",  "MLS", "multidrug", "glycopeptide", "beta-lactam", "bacitracin","aminoglycoside"
                             )) #是为了把others放在最后，multidrug排到红色。
  plotDat$taxon_fct_rev <- factor(plotDat$taxon,
                                  levels = rev(topTaxons))
  library(ggplot2)
  p<-ggplot(plotDat,aes(fill=drugType, y=taxon_fct_rev , x=count))+ 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values = colors) +
    theme_bw() + theme(panel.grid = element_blank()) + xlab("") + ylab("") +
    # scale_x_discrete(labels = paste(topTaxons,"\n#Bins = ", Taxon_rankDf$numBins[1:N], sep = "")) +
    theme(axis.text.x = element_text(angle = 90)) 
  p
  
  pdfName <- paste("Fig7c.ArgCompositionPlot_",ecosystem,"_top10",taxonLevel,"wMost",topTaxonBy,".pdf",sep = "")
  ggsave(p, filename = pdfName, device = "pdf", width = 8, height = 3)
  
}
