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
# transpor tabela
# dar match com filogenia
# calcular taxa de evolução 