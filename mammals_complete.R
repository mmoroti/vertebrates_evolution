# Tue Oct 25 12:46:26 2022 ------------------------------
# package
library(tidyverse)
library(phytools)
library(picante)
library(RRphylo)

# trait data
dir()
setwd("C:/Users/Cliente/Desktop/vertebrates_evolution/traits/mammals")
trait_mammals <- as_tibble(read_delim("trait_data_imputed.csv"))
View(trait_mammals)
trait_mammals <- trait_mammals %>% select("iucn2020_binomial","adult_mass_g", "litter_size_n")
visdat::vis_miss(trait_mammals)
trait_mammals

# phylogeny data
mammals.phy <- read.nexus("https://raw.githubusercontent.com/n8upham/MamPhy_v1/master/_DATA/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_MCC_v2_target.tre")

#--- Checking phylogeny
# In mammals, the tree proposed by Upham needs some nomenclature adjustments to proceed with the analysis.
mammals.upham.drop <- drop.tip(mammals.phy, "_Anolis_carolinensis") #remove Anolis
sp_names_full <- mammals.upham.drop$tip.label
split_names <- sapply(strsplit(sp_names_full, split='_'), function(x) (c(x[1], x[2])))
split_names_t <- data.frame(t(split_names))
new_names <- paste(split_names_t$X1, split_names_t$X2, sep="_")

# check if it is correct
data.frame(original=mammals.upham.drop$tip.label[-1][1:10], fixed=new_names[1:10])
mammals.upham.drop$tip.label <- new_names
mammals.upham.drop #4175 tips

#Trocando . por _ 
trait_mammals$iucn2020_binomial <- gsub(" ", "_",trait_mammals$iucn2020_binomial)

#Removing missing data. This dataset is imputed.
trait_mammals_all <- remove_missing(trait_mammals)

# Transposition
mammals_trans <- t(trait_mammals_all)
colnames(mammals_trans) <- mammals_trans[1,]
mammals_trans <- mammals_trans[-1, ] #remove sp line 
ncol(mammals_trans) # 6028
mammals.upham.drop # 4175
View(mammals_trans)
#--- Match with phylogeny

match.phylo.comm(mammals.upham.drop, mammals_trans)
mammals_phy <- prune.sample(mammals_trans, mammals.upham.drop)
# 4175 - 3769 = 406 spp perdidas

# Traits with position in phylogeny
names_phy <- setNames(data.frame(mammals_phy$tip.label), "iucn2020_binomial") 
nrow(names_phy) # 3769

traits_phy <- left_join(names_phy,trait_mammals,
                        by="iucn2020_binomial") %>% 
                        distinct(.keep_all = TRUE)

teste <- concat(names_phy,trait_mammals,
      by.x="iucn2020_binomial",all.x = TRUE, sort = FALSE)

nrow(traits_phy) # 3772
mammals_phy # precisa ter 3769

###------------
# Tue Nov 15 16:45:18 2022 ------------------------------

# Mammals
mass_input_mammals <- as.numeric(traits_phy$adult_mass_g)
litter_size_input_mammals <- as.numeric(traits_phy$litter_size)

names(mass_input_mammals) <- traits_phy$iucn2020_binomial
names(litter_size_input_mammals) <- traits_phy$iucn2020_binomial
View(litter_size_input_mammals)

# taxas
# body size
rates.mammals.mass <- RRphylo(tree= mammals_phy, y= mass_input_mammals)
rates_mammals_body <- rates.mammals.mass$rates[3769:7537,]
head(rates_mammals_body)
# litter size
rates.mammals.litter <- RRphylo(tree= mammals_phy, y= litter_size_input_mammals)
rates_mammals_litter <- rates.mammals.litter$rates[3769:7537,]
head(rates_mammals_litter)
#salvando
getwd()
write.csv2(rates_mammals_body,'rates_mammals_body.csv')
write.csv2(rates_mammals_litter,'rates_mammals_litter.csv')