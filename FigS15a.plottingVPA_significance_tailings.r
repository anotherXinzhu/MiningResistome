setwd("E:/NutSync/papers/AMD sedi/_EDA_and_Figures/_Figures_withNonRegGenes")
library(dplyr)
library(vegan)
library(xlsx)
library(tibble)

# read physico-chemical properties --------------
env.complete<-fread("_data/physico-chemical.tailings.csv")

# correct env --------------------------

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


# read ARG  ---------------------------
arg <- fread("_data/nonRegARG.abund.csv",data.table = F) %>%
  filter(sampleID %in% rownames(env.corrected)) %>%
  reshape2::dcast(sampleID~ARG, value.var = "DepthPG") %>%
  column_to_rownames("sampleID")

rownames(arg) == rownames(env.corrected)


# env vars kept after optimization are shown below ---------------
optim.rdaVars_env_4grps <- c("avail_K","K","DTPA_Cu","Cu","DTPA_Pb","Pb","PLI_cat","EcoRI_cat", 
                             "pH", "EC","longitude","latitude","position5",
                             "NS","lon_cat","lat_cat","EWM","MAP")

# VPA plot -------------- -------------- --------------
metal.inVpa <- c("avail_K","K","DTPA_Cu","Cu","DTPA_Pb","Pb","PLI_cat","EcoRI_cat")
phychem.inVpa <- c( "pH", "EC")
geo.inVpa <- c("longitude","latitude","position5", "NS","lon_cat","lat_cat")
climate.inVpa <- c("MAP")

env_4grps <- env.corrected
colnames(env_4grps)[which(sapply(env_4grps,class) == "character")]
#convert character to factor 
env_4grps$position5 <- factor(env_4grps$position5)
env_4grps$NS <- factor(env_4grps$NS)
env_4grps$EMW <- factor(env_4grps$EMW)
env_4grps$lon_cat <- factor(env_4grps$lon_cat)
env_4grps$lat_cat  <- factor(env_4grps$lat_cat)
env_4grps$Nemerow_cat <- factor(env_4grps$Nemerow_cat)
env_4grps$PLI_cat <- factor(env_4grps$PLI_cat)
env_4grps$EcoRI_cat <- factor(env_4grps$EcoRI_cat)
env_4grps$TImetal_cat <- factor(env_4grps$TImetal_cat)
env_4grps$TImetal_avail_cat1 <- factor(env_4grps$TImetal_avail_cat1)


# vpa and plot #########
vpt <- varpart(arg, 
               env_4grps[,metal.inVpa], 
               env_4grps[,phychem.inVpa], 
               env_4grps[,c(geo.inVpa,climate.inVpa)])

pdf("FigS14b.vpa.nonRegARGcomposition.tailings.pdf")
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
env <- env.corrected 
fml <- paste("arg ~", paste(metal.inVpa, collapse = " + "), sep = "");fml
rda.metal <- rda(as.formula(fml), env)
anova(rda.metal)$`Pr(>F)`  # 0.001

fml <- paste("arg ~", paste(phychem.inVpa, collapse = " + "), sep = "")
rda.phychem <- rda(as.formula(fml), env)
set.seed(10); anova(rda.phychem)$`Pr(>F)`  #  0.001

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
set.seed(10); anova(rda.phychem.conditioned)$`Pr(>F)`  # 0.014


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


