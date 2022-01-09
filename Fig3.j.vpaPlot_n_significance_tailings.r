setwd( "E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes/")

library(data.table)
env.complete<-fread("_data/physico-chemical.tailings.csv")


# correct env --------------------------
library(dplyr)
library(tibble)
env.corrected <- env.complete %>% column_to_rownames("V1")

# Total N has negative values, convert to min/2
env.corrected$TN[which(env.corrected$TN <= 0)] <- 
  min(env.corrected$TN[which(env.corrected$TN > 0)])/2 # change 0 or minus value to min/2

# log transfer the variables with max/min > 100 
sapply(env.corrected, class)
env.conti <- env.corrected %>% select(where(is.numeric))

bigRangeEnvVars <- colnames(env.conti)[which(sapply(env.conti, function(x) max(x)/min(x)) > 100)]
sapply(env.conti[,bigRangeEnvVars],min)
for(brVars in bigRangeEnvVars){
  #brVars = bigRangeEnvVars[1]
  env.corrected[,brVars]  <- log1p(env.corrected[,brVars])
}

env.corrected <- env.corrected[complete.cases(env.corrected),] # remove NA values

# PCA first  ---------------------------
group.metals <- c("avail_K","K", "DTPA_Cd","Cd", "DTPA_Cu", "Cu", "DTPA_Pb","Pb","DTPA_Zn","Zn",
                  "Nemerow","Nemerow_cat","PLI","PLI_cat", "EcoRI", "EcoRI_cat",
                  "TImetal", "TImetal_cat", "TImetal_avail","TImetal_avail_cat1",  "multiMetalIndex_total", "multiMetalIndex_avail")
group.phychems <- c("TP","pH", "EC", "TN", "TC")
group.climate <- c("MAT","MAP")
group.geo <- c("longitude", "latitude","lon_cat", "lat_cat","position5","EMW","NS")

groups_tPCA <- c("group.metals","group.phychems")

for(g in groups_tPCA){
  # g=groups_tPCA[1]
  envVars <- eval(parse(text = g))
  
  datPCA <- env.corrected %>% select(all_of(envVars)) %>% select(where(is.numeric))
  
  varType = sub("group\\.(\\w+)","\\1",g)
  PCs <- prcomp(scale(datPCA))$x
  colnames(PCs) <- paste(varType,colnames(PCs),sep = "_")
  assign(paste( "PCs_",varType,sep = ""), PCs, envir = .GlobalEnv)
}


rownames(PCs_metals) == rownames(env.corrected)
rownames(PCs_metals) == rownames(PCs_phychems)
env_PCs <- cbind.data.frame(PCs_metals,PCs_phychems, 
                            env.corrected %>% select(where(is.character)),
                            env.corrected %>% select(longitude, latitude, MAP, MAT))


# read ARG  ---------------------------
arg <- fread("_data/nonRegARG.abund.csv") %>%
  filter(sampleID %in% rownames(env.corrected)) %>%
  group_by(sampleID) %>%
  summarise(DepthPG=sum(DepthPG))

arg <- arg[match(rownames(env.corrected), arg$sampleID),]



# env vars kept after optimization are shown below ---------------
metal.inVpa <- c( "metals_PC1",  "metals_PC2" , "metals_PC3" , "metals_PC5" , "metals_PC7" , "metals_PC12", "PLI_cat")
geo.inVpa <- c( "longitude", "lon_cat",   "lat_cat",   "EMW" )
phychem.inVpa <- "phychems_PC4"

# VPA plot -------------- -------------- --------------
library(vegan)
rownames(env_PCs) == arg$sampleID
vpt <- varpart(arg$DepthPG, 
               env_PCs[,metal.inVpa], 
               env_PCs[,phychem.inVpa], 
               env.corrected[,geo.inVpa])

pdf("Fig3j.vpa.tailings.pdf")
plot(
  vpt,
  bg = 2:5,
  id.size = 1.1,
  cex = 1.2,
  Xnames = c('metals', "phychem","geo")
)
title('Env vars explain ARG composition \nin tailings')
dev.off()


# significance -----------------------
env <- env_PCs
arg <- arg %>% column_to_rownames("sampleID")

# calculate significance of each individual env type, including the co-linear part ####
fml <- paste("arg ~", paste(metal.inVpa, collapse = " + "), sep = "");fml
rda.metal <- rda(as.formula(fml), env)
anova(rda.metal)$`Pr(>F)`  # 0.001

fml <- paste("arg ~", paste(phychem.inVpa, collapse = " + "), sep = "")
rda.phychem <- rda(as.formula(fml), env)
set.seed(10); anova(rda.phychem)$`Pr(>F)`  # 0.046

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
set.seed(10); anova(rda.phychem.conditioned)$`Pr(>F)`  # 0.126


fml <-  paste("arg ~",
              paste(geo.inVpa, collapse = " + "), 
              " + Condition(",
              paste(metal.inVpa, collapse = ") + Condition("), 
              ") + Condition(",
              paste(phychem.inVpa, collapse = ") + Condition("),
              ")",
              sep = ""); fml
rda.geo.conditioned <- rda(as.formula(fml), data=env)
set.seed(10); anova(rda.geo.conditioned)$`Pr(>F)`  # 0.031

# calculate significance of the colineared portion (not sure whether correct or not) #####
fml <-  paste("arg ~",
              paste(geo.inVpa, collapse = ":"), 
              ":",
              paste(metal.inVpa, collapse = ":"), 
              sep = "");fml
rda.metal.geo.colinear <- rda(as.formula(fml), data=env)
set.seed(10); anova(rda.metal.geo.colinear)$`Pr(>F)`  # 0.004



fml <-  paste("arg ~",
              paste(geo.inVpa, collapse = ":"), 
              ":",
              paste(metal.inVpa, collapse = ":"),
              ":",
              paste(phychem.inVpa, collapse = ":"),
              sep = "");fml
rda.metal.geo.phychem.colinear <- rda(as.formula(fml), data=env)
set.seed(10); anova(rda.metal.geo.phychem.colinear)$`Pr(>F)`  # 0.002


# partial mantel test -----------------------------------------------------
# calculate distance 
arg <- arg %>% column_to_rownames("sampleID")
rownames(arg) == rownames(env_PCs)
arg.dist <- vegdist(arg, method = 'euclidean',upper = T ,diag = T)
metal.dist <- vegdist(env_PCs[,metal.inVpa] %>% select_if(is.numeric),
                      method = 'euclidean',upper = T ,diag = T)
library(geosphere)
geo.dist <-  distm(env_PCs[,c("longitude","latitude")] ,fun = distHaversine)

# effects of metal parameters on ARG composition, controlling geographic 
mantel.partial(arg.dist, metal.dist, as.dist(geo.dist), method = "spearman", permutations = 9999) 
'Mantel statistic r: 0.1848  ; Significance: 1e-04 ' 

# effects of geographic parameters on ARG composition, controlling metals 
mantel.partial(arg.dist, as.dist(geo.dist), metal.dist, method = "spearman", permutations = 9999) 
'Mantel statistic r:  0.07761 ; Significance: 0.0097  ' 
