# Sat Jan 07 09:39:00 2023 ------------------------------
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

# Squamata dataset
dir()

# Load traits 
# We used the data of Guedes et al., 2022 Ecography
# Guedes used the dataset of Feldman et al. (2015)

# We used mass as our measure of body size instead of other measures, such as snout–vent length or total length, because these cannot be compared easily between clades that differ greatly in their body plan (see, e.g., figure S2c of Feldman et al., 2016, where squamates of similar length differ by two orders of magnitude in mass) (Slavenko et al., 2019)
setwd("C:/Users/Cliente/Desktop/vertebrates_evolution/traits/squamata")
reptile_traits <- read.csv2("polyBasedData.csv", sep=',', na.strings=c("","NA")) # we need to replace blank spaces in the name of species for '_'
names(reptile_traits)

# Selecting traits
reptile_traits_short <- reptile_traits %>% select("Comb_sciname", "BodyMass")
vis_miss(reptile_traits_short) # 6.53% missing data
nrow(reptile_traits_short) # 10531 lines

# Now, we need to replace the blank spaces of body mass for NA's
reptile_traits_short$Comb_sciname <- sub(" ", "_", reptile_traits_short$Comb_sciname)

# Load phylogeny
squamata_tree <- read.tree("squam_shl_new_Consensus_9755.tre")

# remover NA's
<<<<<<< HEAD
body_reptile <- remove_missing(reptile_traits_short)
# 688 spp. removed

# transpor tabela 
body_rep_trans <- t(body_reptile)
View(body_rep_trans)
colnames(body_rep_trans) <- body_rep_trans[1,]

# check datasets
ncol(body_rep_trans) # 9843spp

# dar match com filogenia
squamata_tree # 9755

# Match traits with phylogeny without missing data
squamata_phy_body <- prune.sample(body_rep_trans, 
                                  squamata_tree) # 9487 spp.

# Dataset needs adjustments 
squamata_phy_body # 9487 spp - Perdemos 356 spp na árvore
nrow(body_reptile) # 9843spp - dataset

# remove species in body_reptile 
names_phy <- as.data.frame(squamata_phy_body$tip.label)
names(names_phy) = "Comb_sciname"

# join with names present in phylogeny
body_squa_tree <- left_join(names_phy, body_reptile, by = "Comb_sciname")

nrow(body_squa_tree) # 9487 spp.
squamata_phy_body # 9487 spp.

# Preparing object to calculing trait evolution rate
# Mammals
mass_squamata <- as.numeric(body_squa_tree$BodyMass)
names(mass_squamata) <- body_squa_tree$Comb_sciname
mass_squamata

# calcular taxa de evolução (dataset without NA's)
rates.squamata <- RRphylo(tree= squamata_phy_body, y= mass_squamata)
###
rates_squamata_body <- rates.squamata$rates[5030:10059,]

# Salvando dataset
getwd()
write.csv2(rates_amphibia_body,'rates_amphibia_body.csv')
hist(log10(rates_amphibia_body**2))     

###---
# Imputation script
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


