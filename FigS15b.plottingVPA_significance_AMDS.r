setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")
library(dplyr)
library(vegan)
library(xlsx)
library(tibble)

# read physico-chemical properties --------------
grp1.metals_df <- read.xlsx("_data/physico_chemical.amdsediments.xlsx", sheetName = "metals.log")
grp2.phychem_df <- read.xlsx("_data/physico_chemical.amdsediments.xlsx", sheetName = "other.physico.chem.log")
grp3.geo_df <-  read.xlsx("_data/physico_chemical.amdsediments.xlsx", sheetName = "geo")
grp4.climate_df <-   read.xlsx("_data/physico_chemical.amdsediments.xlsx", sheetName = "climate")

load("_data/1.GeoInformaton.RData")

# load ARG data --------------
ARG_df.l <- fread("_data/nonRegARG.abund.csv", data.table = F) %>% 
  filter(grepl("^LD", sampleID)) %>%
  mutate(sampleName = sapply(sampleID,function(x) AMDsedi_geo_df$sample[AMDsedi_geo_df$sampleID == x]))

all(ARG_df.l$sampleName %in% grp1.metals_df$sample)

arg <- ARG_df.l %>% 
  reshape2::dcast(sampleName~ARG, value.var = "DepthPG") %>%
  column_to_rownames("sampleName") 


rownames(arg) == grp1.metals_df$sample
rownames(arg) == grp2.phychem_df$sample
rownames(arg) == grp3.geo_df$sample
rownames(arg) == grp4.climate_df$sample


# env vars kept after optimization are shown below ---------------
optim.rdaVars_env_4grps<-c("THg","MeHg","MeHg_HgT","Fe2_", "Fe3_","Fe","Cu","Pb","avail_Cu",
                           "avail_Cd","avail_Mn","avail_Pb","Nemerow","PLI", "TImetal_avail","Nemerow_cat",
                           "PLI_cat", "EcoRI_cat", "TImetal_cat","TImetal_avail_cat1","pH","TN",
                           "C_N","lng","lat","position","MAT")


# VPA plot -------------- -------------- --------------
metal.inVpa <- colnames(grp1.metals_df)[colnames(grp1.metals_df) %in% optim.rdaVars_env_4grps]
phychem.inVpa <- colnames(grp2.phychem_df)[colnames(grp2.phychem_df) %in% optim.rdaVars_env_4grps]
geo.inVpa <- colnames(grp3.geo_df)[colnames(grp3.geo_df) %in% optim.rdaVars_env_4grps]
climate.inVpa <- colnames(grp4.climate_df)[colnames(grp4.climate_df) %in% optim.rdaVars_env_4grps]


env_4grps <- merge(merge(grp1.metals_df, grp2.phychem_df,by="sample"),
                   merge(grp3.geo_df, grp4.climate_df, by="sample"),
                   by="sample")
colnames(env_4grps)[which(sapply(env_4grps,class) == "character")]
#convert character to factor 
env_4grps$Nemerow_cat <- factor(env_4grps$Nemerow_cat)
env_4grps$PLI_cat <- factor(env_4grps$PLI_cat)
env_4grps$EcoRI_cat <- factor(env_4grps$EcoRI_cat)
env_4grps$TImetal_avail_cat1 <- factor(env_4grps$TImetal_avail_cat1,
                                       levels = c("l1", "l2", "l3", "l4", "l5", "l6", "l7", "l9", "l10"))
env_4grps$TImetal_cat  <- factor(env_4grps$TImetal_cat)
env_4grps$position <- factor(env_4grps$position)


env_4grps <- env_4grps %>% tibble::column_to_rownames("sample")

# match env and arg matrix
rownames(env_4grps) == rownames(arg)
#env_4grps <- env_4grps[match(rownames(arg), rownames(env_4grps)),]


# vpa and plot
vpt <- varpart(arg, 
               env_4grps[,metal.inVpa], 
               env_4grps[,phychem.inVpa], 
               env_4grps[,c(geo.inVpa,climate.inVpa)])

pdf("FigS14b.vpa.nonRegARGcomposition.amdsediments.pdf")
plot(
  vpt,
  bg = 2:5,
  id.size = 1.1,
  cex = 1.2,
  Xnames = c('metals', "phychem","geo")
)
dev.off()



# significance -----------------------
# calculate significance of each individual env type, including the co-linear part ####
env <- env_4grps 
fml <- paste("arg ~", paste(metal.inVpa, collapse = " + "), sep = "");fml
rda.metal <- rda(as.formula(fml), env)
anova(rda.metal)$`Pr(>F)`  # 0.001

fml <- paste("arg ~", paste(phychem.inVpa, collapse = " + "), sep = "")
rda.phychem <- rda(as.formula(fml), env)
set.seed(10); anova(rda.phychem)$`Pr(>F)`  #  0.443

fml <- paste("arg ~", paste(geo.inVpa, collapse = " + "), sep = "")
rda.geo <- rda(as.formula(fml), env)
set.seed(10); anova(rda.geo)$`Pr(>F)`  # 0.001


# calculate significance of each individual env type, excluding the co-linear part ####
fml <-  paste("arg ~",
              paste(metal.inVpa, collapse = " + "), 
              " + Condition(",
              paste(phychem.inVpa, collapse = ") + Condition("), 
              ") + Condition(",
              paste(geo.inVpa, collapse = ") + Condition("),
              ")",
              sep = ""); fml
rda.metal.conditioned <- rda(as.formula(fml), data=env)
set.seed(10); anova(rda.metal.conditioned)$`Pr(>F)`  # 0.001


fml <-  paste("arg ~",
              paste(phychem.inVpa, collapse = " + "), 
              " + Condition(",
              paste(metal.inVpa, collapse = ") + Condition("), 
              ") + Condition(",
              paste(geo.inVpa, collapse = ") + Condition("),
              ")",
              sep = ""); fml
rda.phychem.conditioned <- rda(as.formula(fml), data=env)
set.seed(10); anova(rda.phychem.conditioned)$`Pr(>F)`  # 0.189


fml <-  paste("arg ~",
              paste(geo.inVpa, collapse = " + "), 
              " + Condition(",
              paste(metal.inVpa, collapse = ") + Condition("), 
              ") + Condition(",
              paste(phychem.inVpa, collapse = ") + Condition("),
              ")",
              sep = ""); fml
rda.geo.conditioned <- rda(as.formula(fml), data=env)
set.seed(10); anova(rda.geo.conditioned)$`Pr(>F)`  # 0.001


