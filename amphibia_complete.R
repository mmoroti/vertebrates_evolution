# Wed Nov 23 16:53:35 2022 ------------------------------
# Amphibia

# I modified the function fitGeospiza so that it receives the two arguments (phy, trait) that will go in the model. The function's new name is fit_modified2. In fit_modified2, I included the aicw function from the geiger package to also return us the AICw in the model comparison.
#---- CREATE a FUNCTION to COMPARE EVOLUTION MODELS
fit_modified2=function(phy, trait){
  #trait=match.arg(trait, c("body_size","litter_size"))
  
  # define set of models to compare
  models=c("BM", "OU", "EB", "white")
  summaries=c("diffusion", "Ornstein-Uhlenbeck", "early burst", "white noise")
  
  ## ESTIMATING measurement error ##
  aic.se=numeric(length(models))
  lnl.se=numeric(length(models))
  w.se=numeric(length(models))
  for(m in 1:length(models)){
    cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep=""), models[m],
        " with SE *** \n", sep="")
    tmp=fitContinuous(phy, trait,SE=NA, model=models[m],
                      bounds=list(SE=c(0,0.5)), ncores=2)
    print(tmp)
    aic.se[m]=tmp$opt$aicc
    lnl.se[m]=tmp$opt$lnL
    w.se[m]=tmp$opt$aic
  }
  
  
  ## ASSUMING no measurement error ##
  aic=numeric(length(models))
  lnl=numeric(length(models))
  w.se=numeric(length(models))
  
  for(m in 1:length(models)){
    cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep=""), models[m],
        " *** \n", sep="")
    tmp=fitContinuous(phy, trait,SE=0,model=models[m], ncores=2)
    print(tmp)
    aic[m]=tmp$opt$aicc
    lnl[m]=tmp$opt$lnL
    w.se[m]=tmp$opt$aic
  }
  
  ## COMPARE AIC ##
  names(aic.se)<-names(lnl.se)<-names(aic)<-names(lnl)<-models
  delta_aic<-function(x) x-x[which(x==min(x))]
  
  # no measurement error
  daic=delta_aic(aic)
  cat("\n\n\n\t\t\t\t*** MODEL COMPARISON: "," *** \n",sep="")
  cat("\tdelta-AIC values for models assuming no measurement error
    \t\t\t\t zero indicates the best model\n\n")
  print(daic, digits=2)
  
  # measurement error
  daic.se=delta_aic(aic.se)
  cat("\n\n\n\n\t\t\t\t*** MODEL COMPARISON: "," ***\n",sep="")
  cat("\t\t   delta-AIC values for models estimating SE
    \t\t\t\t zero indicates the best model\n\n")
  print(daic.se, digits=2)
  cat("\n\n\n")
  
  # wight akaike
  w=aicw(w.se)
  cat("\n\n\n\n\t\t\t\t*** MODEL COMPARISON: "," ***\n",sep="")
  cat("\t\t   w-AIC values for models estimating SE
    \t\t\t\t zero indicates the best model\n\n")
  print(w, digits=2)
  cat("\n\n\n")
  
  res_aicc=rbind(aic, aic.se, daic, daic.se, w$w)
  rownames(res_aicc)=c("AICc","AICc_SE","dAICc", "dAICc_SE", "AICw")
  #print(w)
  return(res_aicc)
}

# Packages
library(picante)
library(RRphylo)
library(Rphylopars)
library(phytools)
library(ape)
library(letsR)
library(raster)
library(visdat)
library(tidyverse)
library(ggExtra)
library(cowplot)
library(geiger)
library(rgdal)

# phylogeny 
setwd("C:/Users/Cliente/Desktop/vertebrates_evolution/phylogeny")
anura.phy <- read.tree("amph_shl_new_Consensus_7238.tre")
anura.phy

# traits amphibio
setwd("C:/Users/Cliente/Desktop/vertebrates_evolution/traits/amphibia")
anura.traits <- read.csv2("AmphiBIO_v1.csv", sep=',')
anura.traits <- anura.traits %>% select("Species", "Body_size_mm")
vis_miss(anura.traits) # 22.8% missing data in body mass
# +75% missing data in "Litter_size_min_n","Litter_size_max_n". We decided to remove this trait

# remove space with underscore in species name
anura.traits$Species <- sub(" ", "_", anura.traits$Species)

# Transposition
amphibia_trans <- t(anura.traits)
colnames(amphibia_trans) <- amphibia_trans[1,]
amphibia_trans_2 <- t(as.data.frame(amphibia_trans[-1, ])) #remove sp line 

# Datasets
ncol(amphibia_trans) # 6776
anura.phy # 7239

# Match traits with phylogeny
match.phylo.comm(anura.phy, amphibia_trans)
amphibia_phy <- prune.sample(amphibia_trans, anura.phy)
amphibia_phy # we lost 829 spp

# Traits with position in phylogeny
names_phy <- setNames(data.frame(amphibia_phy$tip.label), "Species") 
nrow(names_phy) # 6410

traits_anura_phy <- left_join(names_phy,anura.traits,
                        by="Species") %>% 
  distinct(.keep_all = TRUE)

# Confering dataset
nrow(traits_anura_phy) #6410
amphibia_phy # 6410
vis_miss(traits_anura_phy) # ~21% missing data

# Thu Nov 24 14:37:12 2022 ------------------------------
# Herein, I started the imputation process
# First, we need to remove 
body_amph <- remove_missing(traits_anura_phy)
nrow(body_amph) # 5030 spp

# Transposition
body_amph_trans <- t(body_amph)
View(body_amph_trans)
colnames(body_amph_trans) <- body_amph_trans[1,]
#amphibia_trans_2 <- t(as.data.frame(amphibia_trans[-1, ])) #remove sp line 

# Datasets
ncol(body_amph_trans) # 5030spp
anura.phy # 7239
# Match traits with phylogeny without missing data
amphibia_phy_body <- prune.sample(body_amph_trans, 
                                  anura.phy)#5030spp

# Thu Nov 24 14:50:35 2022 ------------------------------
# Calculing the best fit evolution model
# Function fit_modified2 is required
###--- Anura
# The Phylogeny needs to rooted and ultrametric
amphibia_phy_body <- force.ultrametric(amphibia_phy_body)
amphibia_phy_body <- ape::multi2di(amphibia_phy_body)

# Preparing the data
# Anura
svl_amph_sem <- as.numeric(body_amph$Body_size_mm)
names(svl_amph_sem) <- body_amph$Species

# Body models
body_anura_model <- fit_modified2(amphibia_phy_body, svl_amph_sem)
body_anura_model # dAICc OU model - AICw 1.0 OU
?fit_modified

###
# Input data
traits_anura_phy$Body_size_mm <- log10(as.numeric(traits_anura_phy$Body_size_mm))

input_amphibia <- phylopars(trait_data = traits_anura_phy, tree = teste,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE, model="mvOU")

###
# Calculing without NA's
svl_amph_sem # without NA's
amphibia_phy_body # without NA's

rates.amphibia <- RRphylo(tree= amphibia_phy_body, y= svl_amph_sem)
rates_amphibia_body <- rates.amphibia$rates[5030:10059,]

# Salvando dataset
#salvando
getwd()
write.csv2(rates_amphibia_body,'rates_amphibia_body.csv')
hist(log10(rates_amphibia_body**2))     
