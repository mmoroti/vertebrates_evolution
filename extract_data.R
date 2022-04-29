# Mon Mar 07 14:32:23 2022 ------------------------------
# Loading packages
library(picante)
library(RRphylo)
library(Rphylopars)
library(ape)
library(letsR)
library(visdat)
library(tidyverse)
library(gapminder)

#--- Extract data ---#
# phylogeny 
setwd("E:/OneDrive/2019 - Moroti/Tese/Capitulo 1/RMarkdown/phylogeny")
anura.phy <- read.tree("phy_amphibia.tre")
reptile.phy <- read.tree("phy_reptile.tre")
birds.phy <- read.tree("phy_birds.tre")
mammals.phy <- 

# functional traits
setwd("E:/OneDrive/2019 - Moroti/Tese/Capitulo 1/RMarkdown/trait_data")
amphibia_full_traits <- as.tibble(read.csv2("traits_amphibia.csv"))
reptile_full_traits <- as.tibble(read.csv2("traits_squamatas.csv"))
birds_full_traits <- as.tibble(read.csv2("traits_birds.csv"))
mammals_full_traits <- as.tibble(read.csv2("traits_mammals.csv"))

# Fri Apr 29 12:56:25 2022 ------------------------------
#--- Transform ---#

#--- Adjusting in nomenclature
amphibia_full_traits$sp <- amphibia_full_traits$sp %>% str_replace_all( "Boana_claresignata", "Bokermannohyla_claresignata") %>% str_replace_all("Boana_clepsydra", "Bokermannohyla_clepsydra") %>% str_replace_all("Proceratophrys_salvatori", "Odontophrynus_salvatori") %>%  str_replace_all("Scinax_v_signatus", "Scinax_v-signatus") %>% str_replace_all("Scinax_x_signatus", "Scinax_x-signatus") %>% str_replace_all("Boana_bandeirantes", "Hypsiboas_bandeirantes") %>% str_replace_all("Boana_poaju", "Hypsiboas_poaju") %>% str_replace_all("Boana_paranaiba", "Hypsiboas_paranaiba") %>% str_replace_all("Physalaemus_nattereri", "Eupemphix_nattereri") %>% str_replace_all("Agalychnis_aspera", "Hylomantis_aspera") %>% str_replace_all("Agalychnis_granulosa", "Hylomantis_granulosa") %>% str_replace_all("Lithobates_palmipes", "Rana_palmipes") %>% str_replace_all("Pithecopus_azureus", "Phyllomedusa_azurea") %>% str_replace_all("Pithecopus_ayeaye", "Phyllomedusa_ayeaye") %>% str_replace_all("Pithecopus_nordestinus", "Phyllomedusa_nordestina") %>% str_replace_all("Pithecopus_rohdei", "Phyllomedusa_rohdei") %>% str_replace_all("Pithecopus_rohdei", "Phyllomedusa_rohdei") %>% str_replace_all("Aparasphenodon_ararapa", "Aparasphenodon_arapapa") %>% str_replace_all("Ololygon_tripui", "Scinax_tripui") %>% str_replace_all("Ololygon_cosenzai", "Scinax_cosenzai") %>% str_replace_all("Ololygon_strigilata", "Scinax_strigilatus") %>% str_replace_all("Pithecopus_hypochondrialis", "Phyllomedusa_hypochondrialis") %>% str_replace_all("Pithecopus_megacephalus","Phyllomedusa_megacephala") %>% str_replace_all("Pithecopus_hypochondrialis", "Phyllomedusa_hypochondrialis") 

#--- Selecting continuous variables
amphibia_traits <- amphibia_full_traits %>% dplyr::select(c(sp, Body_size_mm,Brood.Size)) %>% rename(body_size = Body_size_mm, litter_size = Brood.Size)
vis_miss(amphibia_traits) # 64.64% NA in brood size ~ 3.72% body size

reptile_traits <- reptile_full_traits %>% dplyr::select(c(sp, SVL..mm.,litter_or_clutch_size_n)) %>% rename(body_size = SVL..mm., litter_size = litter_or_clutch_size_n)
vis_miss(reptile_traits) # 57.91% NA in brood size ~ 14.84% body size

birds_traits <- birds_full_traits %>% dplyr::select(c(sp, BodyMass.Value,litter_or_clutch_size_n)) %>% rename(body_mass = BodyMass.Value, litter_size = litter_or_clutch_size_n)
vis_miss(birds_traits) # 29.51% NA in brood size ~ 9.62% body size

mammals_traits <- mammals_full_traits %>% dplyr::select(c(sp, body_mass_g.avg.,Littersize)) %>% rename(body_mass = body_mass_g.avg., litter_size = Littersize)
vis_miss(mammals_traits) # 0.92% NA in brood size ~ 22.46% body size

#--- Checking phylogenis
# Transposition
amphibia_trans <- t(amphibia_traits_short)
colnames(amphibia_trans) <- amphibia_trans[1,]
amphibia_trans <- amphibia_trans[-1, ] #remove sp line 
amphibia_trans <- amphibia_trans[,-34 ]#Aplastodiscus_weygoldti two times
#View(amphibia_trans)

#
match.phylo.comm(anura.phy, amphibia_trans)
anura_phy <- prune.sample(amphibia_trans, anura.phy)
ncol(amphibia_trans) #562 cols or spp
anura_phy #562 tips

#
rem.col.phy <- c("Brachycephalus_sulfuratus","Chiasmocleis_alagoana","Crossodactylus_fransciscanus","Crossodactylus_timbuhy","Crossodactylus_werneri","Dendropsophus_bromeliaceus","Eleutherodactylus_bilineatus","Melanophryniscus_milanoi","Melanophryniscus_xanthostomus","Scinax_strigilatus","Ololygon_tupinamba","Phyllodytes_megatympanum","Physalaemus_atim","Proceratophrys_mantiqueira","Proceratophrys_phyllostoma","Proceratophrys_pombali","Proceratophrys_tupinamba","Pseudopaludicola_atragula","Pseudopaludicola_jaredi","Pseudopaludicola_pocoto","Rhinella_sebbeni","Scinax_caissara","Scinax_melanodactylus","Scinax_rossaferesae","Scinax_skuki","Sphaenorhynchus_canga","Trachycephalus_typhonius","Vitreorana_baliomma")

amphibia_traits_short <- amphibia_traits[!amphibia_traits$sp %in% rem.col.phy,]
amphibia_traits_short <- amphibia_traits_short[-34,]
#nrow(amphibia_traits_short)

# Fri Apr 29 17:06:10 2022 ------------------------------
### Testando a imputação 
test <- as.data.frame(amphibia_traits_short)

test <- test %>% rename(species = sp) 
test$body_size <- as.numeric(test$body_size)
test$litter_size <- as.numeric(test$litter_size)
 View(test)

#list1 <- list(anura_phy)
#list2 <- list(test)
#appended_list <- append(list1, list2)
#View(appended_list)

# Imputação dados faltantes - SVL/Body mass - Número de prole
#Brownian motion
p_BM <- phylopars(trait_data = test, tree = anura_phy,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE)

p_BM$anc_recon[1:20,] # Data with imputed species means
p_BM$anc_var[1:20,] # Variances for each estimate
p_BM$anc_recon[1:20,] - sqrt(p_BM$anc_var[1:20,])*1.96 # Lower 95% CI
p_BM$anc_recon[1:20,] + sqrt(p_BM$anc_var[1:20,])*1.96 # Upper 95% CI

p_OU <- phylopars(trait_data = test, tree = anura_phy,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE, model="mvOU")



###
#--- Taxa de evolução dos atributos
svl <- log10(as.numeric(amphibia_traits_short$body_size))
names(svl) <- amphibia_traits_short$sp
svl
# Extrair váriaveis climáticas - LetsR
# modelos PGLS