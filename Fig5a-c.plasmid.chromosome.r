setwd("E:/yxz_nuts/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")


library(dplyr)
library(ggplot2)
load("_data/plasmid.sourceData.nonRegARGs.RData")


# pubAMD 画图------------------------------------------------------------------------------
pubAMD_ARGQuant_df1 <- pubAMD_df %>% 
  filter(plasmid != "unclassified") %>% 
  filter(drug_type != "unclassified") %>%
  dplyr::group_by( drug_type, plasmid) %>% 
  dplyr::summarise(n=n()) %>% 
  mutate(Freq=n/sum(n))


plotDat_pubAMD <- pubAMD_ARGQuant_df1 


# arrange x by desc plasmid proportion
plotDat_pubAMD$drug_type <- factor(plotDat_pubAMD$drug_type, 
                                   as.character((plotDat_pubAMD %>% filter(plasmid == "plasmid") %>% arrange(desc(Freq)))$drug_type))


P_pubAMD  <- ggplot(plotDat_pubAMD, aes(x = drug_type, y = Freq, fill = plasmid)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "Freq_byNum", x = "", title = "") +
  scale_fill_manual(values = c("#0073C299", "#EFC00099")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) +
  theme(axis.text.x = element_text(angle = 90))


P_pubAMD


# Tailing 画图------------------------------------------------------------------------------
Tailing_ARGQuant_df1 <- Tailing_df %>%
  filter(plasmid != "unclassified") %>% 
  filter(drug_type != "unclassified") %>%
  dplyr::group_by(drug_type, plasmid) %>%
  dplyr::summarise(n=n()) %>%
  mutate(Freq=n/sum(n))


plotDat_Tailing <- Tailing_ARGQuant_df1  # freq by number
#plotDat_Tailing <- Tailing_ARGQuant_df2  # freq by depth

# arrange x by desc plasmid proportion 
plotDat_Tailing$drug_type <- factor(plotDat_Tailing$drug_type, 
                                    as.character((plotDat_Tailing %>% filter(plasmid == "plasmid") %>% arrange(desc(Freq)))$drug_type))


P_Tailing  <- ggplot(plotDat_Tailing, aes(x = drug_type, y = Freq, fill = plasmid)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "Freq_byNum", x = "", title = "") +
  scale_fill_manual(values = c("#0073C299", "#EFC00099")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) +
  theme(axis.text.x = element_text(angle = 90))

P_Tailing



# AMDSedi 画图------------------------------------------------------------------------------
AMDSedi_ARGQuant_df1 <- AMDSedi_df %>% 
  filter(plasmid != "unclassified") %>% 
  filter(drug_type != "unclassified") %>%
  dplyr::group_by(drug_type, plasmid) %>% 
  dplyr::summarise(n=n()) %>% 
  mutate(Freq=n/sum(n))


plotDat_AMDSedi <- AMDSedi_ARGQuant_df1  # freq by number


# arrange x by desc plasmid proportion 
plotDat_AMDSedi$drug_type <- factor(plotDat_AMDSedi$drug_type, 
                                    as.character((plotDat_AMDSedi %>% filter(plasmid == "plasmid") %>% arrange(desc(Freq)))$drug_type))

P_AMDSedi  <- ggplot(plotDat_AMDSedi, aes(x = drug_type, y = Freq, fill = plasmid)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "Freq_byNum", x = "", title = "") +
  scale_fill_manual(values = c("#0073C299", "#EFC00099")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) +
  theme(axis.text.x = element_text(angle = 90))


P_AMDSedi


library(ggpubr)
ggsave(ggarrange(P_pubAMD + theme(legend.position = "none"),
                 P_Tailing + theme(legend.position = "none"),
                 P_AMDSedi + theme(legend.position = "none"),
                 nrow = 1,ncol = 3),
       filename = "Fig6a-c.plasmChrom_inDrugType.pdf", width = 9,height = 4)



# plasmid vs chromosome in overall -----------------------
plotDatoverall_pubAMD <- pubAMD_df %>%
  filter(plasmid != "unclassified") %>%
  filter(drug_type != "unclassified") %>%
  dplyr::group_by( plasmid) %>% 
  dplyr::summarise(n=n()) %>% 
  mutate(Freq=n/sum(n))

plotDatoverall_Tailing <- Tailing_df %>% 
  filter(plasmid != "unclassified") %>%
  filter(drug_type != "unclassified") %>%
  dplyr::group_by( plasmid) %>%
  dplyr::summarise(n=n()) %>% 
  mutate(Freq=n/sum(n))


plotDatoverall_AMDSedi <- AMDSedi_df %>% 
  filter(plasmid != "unclassified") %>%
  filter(drug_type != "unclassified") %>%
  dplyr::group_by( plasmid) %>% 
  dplyr::summarise(n=n()) %>% 
  mutate(Freq=n/sum(n))

plotDatoverall <- rbind.data.frame(plotDatoverall_pubAMD %>% mutate(dataCollection = "pubAMD"),
                                   plotDatoverall_Tailing %>% mutate(dataCollection = "Tailing"),
                                   plotDatoverall_AMDSedi %>% mutate(dataCollection="AMDSedi"),
                                   stringsAsFactors = F)
plotDatoverall$dataCollection <- factor(plotDatoverall$dataCollection, 
                                        levels = c("pubAMD","Tailing","AMDSedi"))

P_overall <- ggplot(plotDatoverall, aes(x = 0, y = Freq, fill = plasmid)) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_bw() +
  facet_wrap(vars(dataCollection)) +
  labs(y = "Freq_byNum", x = "", title = "") +
  scale_fill_manual(values = c("#0073C299", "#EFC00099")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) +
  theme(axis.text.x = element_text(angle = 90))


ggsave(P_overall,
       filename = "Fig6a-c.plasmChrom_overall.pdf", width = 3.5,height = 3.2)

