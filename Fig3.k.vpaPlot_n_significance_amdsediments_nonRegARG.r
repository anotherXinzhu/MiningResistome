setwd( "E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes/")

# read physico-chemical properties --------------
library(xlsx)
grp1.metals_df <- read.xlsx("_data/physico_chemical.amdsediments.xlsx", sheetName = "metals.log")
grp2.phychem_df <- read.xlsx("_data/physico_chemical.amdsediments.xlsx", sheetName = "other.physico.chem.log")
grp3.geo_df <-  read.xlsx("_data/physico_chemical.amdsediments.xlsx", sheetName = "geo")
grp4.climate_df <-   read.xlsx("_data/physico_chemical.amdsediments.xlsx", sheetName = "climate")

load("_data/1.GeoInformaton.RData")

# load ARG data --------------
library(data.table)
ARG_df.l <- fread("_data/nonRegARG.abund.csv") %>% 
  group_by(sampleID) %>% summarise(DepthPG= sum(DepthPG)) %>%
  filter(grepl("^LD", sampleID)) %>%
  mutate(sampleName = sapply(sampleID,function(x) AMDsedi_geo_df$sample[AMDsedi_geo_df$sampleID == x]))

all(ARG_df.l$sampleName %in% grp1.metals_df$sample)

arg <- ARG_df.l %>% column_to_rownames("sampleName") %>% select(-sampleID)


# env vars kept after optimization are shown below ---------------
optim.rdaVars_env_4grps<-c("THg","MeHg","MeHg_HgT","Fe3_","Fe","Pb","avail_Cu","avail_Mn","TImetal_avail","Nemerow_cat",
                           "TImetal_cat","TImetal_avail_cat1","C_N","lng","position","MAP")



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
env_4grps$position <- factor(env_4grps$position)
env_4grps$TImetal_avail_cat1 <- factor(env_4grps$TImetal_avail_cat1,
                                       levels = c("l1", "l2", "l3", "l4", "l5", "l6", "l7", "l9", "l10"))
env_4grps$TImetal_cat  <- factor(env_4grps$TImetal_cat)


env_4grps <- env_4grps %>% tibble::column_to_rownames("sample")

# match env and arg matrix
rownames(env_4grps) == rownames(arg)
env_4grps <- env_4grps[match(rownames(arg), rownames(env_4grps)),]


# vpa and plot ####
vpt <- varpart(arg, 
               env_4grps[,metal.inVpa], 
               env_4grps[,phychem.inVpa], 
               env_4grps[,c(geo.inVpa,climate.inVpa)])


pdf("Fig3j.vpa.amdsediments.pdf")
plot(
  vpt,
  bg = 2:5,
  id.size = 1.1,
  cex = 1.2,
  Xnames = c('metals', "phychem","geo")
)
title('Env vars explain ARG composition \nin AMD sediments')
dev.off()



# significance -----------------------
env <- env_4grps
geo.inVpa <- c(geo.inVpa, climate.inVpa)

# calculate significance of each individual env type, including the co-linear part ####

fml <- paste("arg ~", paste(metal.inVpa, collapse = " + "), sep = "");fml
rda.metal <- rda(as.formula(fml), env)
anova(rda.metal)$`Pr(>F)`  # 0.001

fml <- paste("arg ~", paste(phychem.inVpa, collapse = " + "), sep = "")
rda.phychem <- rda(as.formula(fml), env)
set.seed(10); anova(rda.phychem)$`Pr(>F)`  # 0.036

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
set.seed(10); anova(rda.phychem.conditioned)$`Pr(>F)`  # 0.019


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

# calculate significance of the colineared portion (not sure whether correct or not) #####
fml <-  paste("arg ~",
              paste(geo.inVpa, collapse = ":"), 
              ":",
              paste(metal.inVpa, collapse = ":"), 
              sep = "");fml
rda.metal.geo.colinear <- rda(as.formula(fml), data=env)
set.seed(10); anova(rda.metal.geo.colinear)$`Pr(>F)`  # 0.13 



fml <-  paste("arg ~",
              paste(geo.inVpa, collapse = ":"), 
              ":",
              paste(phychem.inVpa, collapse = ":"), 
              sep = "");fml
rda.geo.phychem.colinear <- rda(as.formula(fml), data=env)
set.seed(10); anova(rda.geo.phychem.colinear)$`Pr(>F)`  # 0.017

fml <-  paste("arg ~",
              paste(metal.inVpa, collapse = ":"), 
              ":",
              paste(phychem.inVpa, collapse = ":"), 
              sep = "");fml
rda.metal.phychem.colinear <- rda(as.formula(fml), data=env)
set.seed(10); anova(rda.geo.phychem.colinear)$`Pr(>F)`  # 0.017



# partial mantel test -----------------------------------------------------
# calculate distance 
rownames(arg) == rownames(env_4grps)
arg.dist <- vegdist(arg, method = 'bray', upper = T ,diag = T)
metal.dist <- vegdist(env_4grps[,metal.inVpa] %>% select_if(is.numeric),
                      method = 'euclidean',upper = T ,diag = T)
library(geosphere)
geo.dist <-  distm(env_4grps[,c("lng","lat")])

# effects of metal parameters on ARG composition, controlling geographic 
mantel.partial(arg.dist, metal.dist, as.dist(geo.dist), method = "spearman", permutations = 9999) 
'Mantel statistic r: 0.2139  ; Significance: 1e-04 ' 

# effects of geographic parameters on ARG composition, controlling metals 
mantel.partial(arg.dist, as.dist(geo.dist), metal.dist, method = "spearman", permutations = 999) 
'Mantel statistic r:   ; Significance:  ' 
