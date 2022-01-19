setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")

library(ggplot2)
library(ggsci)
library(scales)
library(dplyr)
library(data.table)

dat.original <- fread("_data/nonRegARG.bins.txt",data.table=F)
load("_data/otherARGtypes.RData")

# colors -----------------------------
drugType_colorsDf <- cbind.data.frame(
  Colors=c(pal_npg("nrc")(10), "gray"),
  DrugTypes=c("aminoglycoside", "bacitracin", "beta-lactam","glycopeptide",
              "MLS","tetracycline","phenicol","multidrug",
              "mupirocin","sulfonamide","other"), 
  stringsAsFactors=F
)


# for loops and plotting ----------------------------------

for(ecosystem in c("AMDSedi","Tailing")){
  
  
  # ecosystem = "Tailing"
  dat <- dat.original %>% filter(grepl(ecosystem, binFileName))
  colnames(dat)
  
  
  for(taxonLevel in c("Family","Genus")){  #c("Phylum","Family","Genus")
    
    # taxonLevel = "Family"
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
    
    N = unname(sapply(taxonLevel, function(x) if(x == "Family") 20 else if(x == "Genus") 50))
    
    topTaxons <- Taxon_rankDf$taxon[1:N]
    
    # define other ARG types -------------------------------------------
    dat_ArgComposition$drugType[which(dat_ArgComposition$drugType %in% otherARGtypes)] <- "other" 
    
    
    
    plotDat <- dat_ArgComposition %>% filter(taxon %in% topTaxons) %>% filter(drugType != "unclassified") %>%
      group_by(taxon, drugType) %>% summarise(count=n()) %>% mutate(freq = count/sum(count)) %>% 
      mutate(taxon_fct = factor(taxon, levels = topTaxons)) 
    
    plotDat$drugType <- factor(plotDat$drugType, 
                               levels = c("aminoglycoside", "bacitracin", "beta-lactam", "glycopeptide", "multidrug", "MLS", "mupirocin",
                                          "phenicol","sulfonamide" , "tetracycline", "other" )) 
    
    
    colors <- sapply(levels(plotDat$drugType), 
                     function(x) drugType_colorsDf$Colors[which(drugType_colorsDf$DrugTypes== x)])
    
    
    
    plotDat$drugType <- factor(plotDat$drugType, 
                               levels = c("other", "tetracycline","sulfonamide","phenicol",  
                                          "mupirocin",  "MLS", "multidrug", "glycopeptide", "beta-lactam", "bacitracin","aminoglycoside"
                               ))
    plotDat$taxon_fct_rev <- factor(plotDat$taxon,
                                    levels = rev(topTaxons))
    p<-ggplot(plotDat,aes(fill=drugType, y=taxon_fct_rev , x=count))+ 
      geom_bar(position="fill", stat="identity") +
      scale_fill_manual(values = colors) +
      theme_bw() + theme(panel.grid = element_blank()) + xlab("") + ylab("") +
      theme(axis.text.x = element_text(angle = 90)) 
    
    
    
    p
    pdfName <- paste("FigS17.ArgComposition_",ecosystem,"_top",N,taxonLevel,".pdf",sep = "")
    pdfName
    pdfWidth <- unname(sapply(N,function(x) if(x==10) 4 else if(x==20) 7 else if(x==30) 10 else if(x==50) 13))
    ggsave(p, filename = pdfName, device = "pdf", width = 7, height = pdfWidth)
    
    
    
  }# loop through taxonlevel
  
}# loop through ecosystems

