# Mon Mar 07 14:32:23 2022 ------------------------------
# Loading packages
library(picante)
library(RRphylo)
library(Rphylopars)
library(phytools)
library(ape)
library(letsR)
library(visdat)
library(tidyverse)
library(ggExtra)

#--- Extract data ---#
# phylogeny 
setwd("E:/OneDrive/2019 - Moroti/Tese/Capitulo 1/RMarkdown/phylogeny")
anura.phy <- read.tree("phy_amphibia.tre")
reptile.phy <- read.tree("phy_reptile.tre")
birds.phy <- read.tree("phy_birds.tre")
mammals.phy <- read.nexus("https://raw.githubusercontent.com/n8upham/MamPhy_v1/master/_DATA/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_MCC_v2_target.tre") 

# functional traits
setwd("E:/OneDrive/2019 - Moroti/Tese/Capitulo 1/RMarkdown/trait_data")
amphibia_full_traits <- as_tibble(read.csv2("traits_amphibia.csv"))
reptile_full_traits <- as_tibble(read.csv2("traits_squamatas.csv"))
birds_full_traits <- as_tibble(read.csv2("traits_birds.csv"))
mammals_full_traits <- as_tibble(read.csv2("traits_mammals.csv"))

# Fri Apr 29 12:56:25 2022 ------------------------------
#--- Transform ---#
#--- Adjusting in nomenclature
amphibia_full_traits$sp <- amphibia_full_traits$sp %>% str_replace_all( "Boana_claresignata", "Bokermannohyla_claresignata") %>% str_replace_all("Boana_clepsydra", "Bokermannohyla_clepsydra") %>% str_replace_all("Proceratophrys_salvatori", "Odontophrynus_salvatori") %>%  str_replace_all("Scinax_v_signatus", "Scinax_v-signatus") %>% str_replace_all("Scinax_x_signatus", "Scinax_x-signatus") %>% str_replace_all("Boana_bandeirantes", "Hypsiboas_bandeirantes") %>% str_replace_all("Boana_poaju", "Hypsiboas_poaju") %>% str_replace_all("Boana_paranaiba", "Hypsiboas_paranaiba") %>% str_replace_all("Physalaemus_nattereri", "Eupemphix_nattereri") %>% str_replace_all("Agalychnis_aspera", "Hylomantis_aspera") %>% str_replace_all("Agalychnis_granulosa", "Hylomantis_granulosa") %>% str_replace_all("Lithobates_palmipes", "Rana_palmipes") %>% str_replace_all("Pithecopus_azureus", "Phyllomedusa_azurea") %>% str_replace_all("Pithecopus_ayeaye", "Phyllomedusa_ayeaye") %>% str_replace_all("Pithecopus_nordestinus", "Phyllomedusa_nordestina") %>% str_replace_all("Pithecopus_rohdei", "Phyllomedusa_rohdei") %>% str_replace_all("Pithecopus_rohdei", "Phyllomedusa_rohdei") %>% str_replace_all("Aparasphenodon_ararapa", "Aparasphenodon_arapapa") %>% str_replace_all("Ololygon_tripui", "Scinax_tripui") %>% str_replace_all("Ololygon_cosenzai", "Scinax_cosenzai") %>% str_replace_all("Ololygon_strigilata", "Scinax_strigilatus") %>% str_replace_all("Pithecopus_hypochondrialis", "Phyllomedusa_hypochondrialis") %>% str_replace_all("Pithecopus_megacephalus","Phyllomedusa_megacephala") %>% str_replace_all("Pithecopus_hypochondrialis", "Phyllomedusa_hypochondrialis") 

# Squamata
# Birds
# Mammals
mammals_full_traits$sp <- str_replace_all(mammals_full_traits$sp,"Lycalopex", "Pseudalopex") %>% str_replace_all("Puma_yagouaroundi","Herpailurus_yagouaroundi")

#--- Selecting continuous variables
amphibia_traits <- amphibia_full_traits %>% dplyr::select(c(sp, Body_size_mm,Brood.Size)) %>% rename(body_size = Body_size_mm, litter_size = Brood.Size)
vis_miss(amphibia_traits) # 64.64% NA in brood size ~ 3.72% body size

reptile_traits <- reptile_full_traits %>% dplyr::select(c(sp, SVL..mm.,litter_or_clutch_size_n)) %>% rename(body_size = SVL..mm., litter_size = litter_or_clutch_size_n)
vis_miss(reptile_traits) # 57.91% NA in brood size ~ 14.84% body size

birds_traits <- birds_full_traits %>% dplyr::select(c(sp, BodyMass.Value,litter_or_clutch_size_n)) %>% rename(body_mass = BodyMass.Value, litter_size = litter_or_clutch_size_n)
vis_miss(birds_traits) # 29.51% NA in litter_size ~ 9.62% body_mass

mammals_traits <- mammals_full_traits %>% dplyr::select(c(sp, body_mass_g.avg.,Littersize)) %>% rename(body_mass = body_mass_g.avg., litter_size = Littersize) %>% distinct(sp, .keep_all = TRUE)
vis_miss(mammals_traits) # 1.2% NA in body_mass ~ 26.91% litter_size

#--- Checking phylogenis
# Transposition
amphibia_trans <- t(amphibia_traits_short)
colnames(amphibia_trans) <- amphibia_trans[1,]
amphibia_trans <- amphibia_trans[-1, ] #remove sp line 
amphibia_trans <- amphibia_trans[,-34 ]#Aplastodiscus_weygoldti two times
#View(amphibia_trans)
View(amphibia_trans)
#--- Match with phylogeny
match.phylo.comm(anura.phy, amphibia_trans)
anura_phy <- prune.sample(amphibia_trans, anura.phy)

# Check data
ncol(amphibia_trans) #562 cols or spp
anura_phy #562 tips

#--- Remove species in anura traits list
rem.col.phy <- c("Brachycephalus_sulfuratus","Chiasmocleis_alagoana","Crossodactylus_fransciscanus","Crossodactylus_timbuhy","Crossodactylus_werneri","Dendropsophus_bromeliaceus","Eleutherodactylus_bilineatus","Melanophryniscus_milanoi","Melanophryniscus_xanthostomus","Scinax_strigilatus","Ololygon_tupinamba","Phyllodytes_megatympanum","Physalaemus_atim","Proceratophrys_mantiqueira","Proceratophrys_phyllostoma","Proceratophrys_pombali","Proceratophrys_tupinamba","Pseudopaludicola_atragula","Pseudopaludicola_jaredi","Pseudopaludicola_pocoto","Rhinella_sebbeni","Scinax_caissara","Scinax_melanodactylus","Scinax_rossaferesae","Scinax_skuki","Sphaenorhynchus_canga","Trachycephalus_typhonius","Vitreorana_baliomma")

# Remove
amphibia_traits_short <- amphibia_traits[!amphibia_traits$sp %in% rem.col.phy,]
# Remove Aplastodiscus_weygoldti (appearing twice)
amphibia_traits_short <- amphibia_traits_short[-34,]
class(amphibia_traits_short)
#nrow(amphibia_traits_short)

# Squamatas
# Tue May 17 16:38:00 2022 ------------------------------
reptile_trans <- t(reptile_traits_short)
colnames(reptile_trans) <- reptile_trans[1,]
reptile_trans <- reptile_trans[-1, ] #remove sp line 

#--- Match with phylogeny
match.phylo.comm(reptile.phy, reptile_trans)
squamata_phy <- prune.sample(reptile_trans, reptile.phy)

# Check data
ncol(reptile_trans) #405
squamata_phy #405 tips

#--- Remove species in squamata traits list
rem.col.phy.rep <- c("Leposternon_polystegum ","Amerotyphlops_arenensis","Amphisbaena_metallurga", "Kinosternon_scorpioides","Paleosuchus_palpebrosus","Tomodon_dorsatus")
nrow(reptile_traits_short)

# Remove
reptile_traits_short <- reptile_traits[!reptile_traits$sp %in% rem.col.phy.rep,]
reptile_traits_short <- reptile_traits_short[-36,]# Remove Leposternon_polystegum

# Tue May 17 18:02:55 2022 ------------------------------
# Birds
birds_traits
birds.phy

#--- Checking phylogenis
# Transposition
birds_trans <- t(short_birds_traits)
colnames(birds_trans) <- birds_trans[1,]
birds_trans <- birds_trans[-1, ] #remove sp line 
#birds_trans <- birds_trans[,-34 ]#Aplastodiscus_weygoldti two times
View(birds_trans)

#--- Match with phylogeny
match.phylo.comm(birds.phy, birds_trans)
birds_phy <- prune.sample(birds_trans, birds.phy)

View(birds_trans)

# remove names in traits
rem.birds.phy <- c("Anodorhynchus_glaucus","Aramides_cajaneus","Campylopterus_calcirupicola","Cichlocolaptes_mazarbarnetti","Formicivora_acutirostris","Ortalis_squamata","Sporophila_maximiliani","Synallaxis_cinerea","Vireo_chivi","Porphyriops_melanops","Celeus_galeatus","Suiriri_affinis","Paraclaravis_geoffroyi","Chlorestes_notata","Chordeiles_nacunda","Strix_huhula","Strix_virgata","Heliodoxa_rubricauda","Diopsittaca_nobilis","Hydropsalis_anomala","Nycticryphes_semicollaris","Eupsittula_aurea","Eupsittula_cactorum","Pyrrhocoma_ruficeps","Geranoaetus_albicaudatus","Spinus_magellanicus","Pionus_reichenowi","Vanellus_cayanus","Hydropsalis_maculicaudus","Myiothlypis_leucophrys","Mareca_sibilatrix","Procacicus_solitarius","Cercomacroides_laeta","Phylloscartes_eximius","Anumara_forbesi","Cyanoloxia_brissonii","Griseotyrannus_aurantioatrocristatus","Poospiza_thoracica","Pachysylvia_muscicapina","Myrmoderus_loricatus","Rhopias_gularis","Stigmatura_bahiae","Pseudopipra_pipra","Ceratopipra_rubrocapilla","Philohydor_lictor","Microspingus_lateralis","Microspingus_melanoleucus","Tangara_ornata","Tangara_sayaca","Xenops_rutilans","Hydropsalis_forcipata","Systellura_longirostris","Nyctipolus_hirundinaceus","Podicephorus_major","Mustelirallus_albicollis","Psittacara_leucophthalmus","Ramphastos_vitellinus","Laterallus_viridis","Rupornis_magnirostris","Sarkidiornis_sylvicola","Hydropsalis_parvula","Anas_platalea","Anas_versicolor","Stephanoxis_loddigesii","Sternula_superciliaris","Formicivora_paludicola","Elaenia_sordida","Myiodynastes_maculatus","Myrmotherula_axillaris","Tyranniscus_burmeisteri","Sclerurus_cearensis","Scytalopus_gonzagai","Scytalopus_petrophilus","Sporophila_pileata","Lipaugus_conditus","Tityra_cayana","Xiphorhynchus_atlanticus","Microspingus_cinereus","Tangara_brasiliensis","Pipraeidea_bonariensis","Myrmoderus_squamosus","Anabacerthia_lichtensteini","Ortalis_araucuan","Orthopsittaca_manilatus","Parabuteo_leucorrhous","Pseudastur_polionotus","Xiphorhynchus_guttatoides","Agelaioides_fringillarius","Basileuterus_auricapilla","Myiothlypis_rivularis","Tangara_cyanomelas","Tangara_flava","Asthenes_moreirae","Automolus_lammi","Icterus_pyrrhopterus","Calidris_subruficollis","Amazilia_sapphirina","Antrostomus_sericocaudatus","Buteogallus_coronatus","Amadonastur_lacernulatus","Heterospizias_meridionalis","Celeus_tinnunculus","Synallaxis_hellmayri","Herpsilochmus_scapularis","Clibanornis_rectirostris","Sporophila_beltoni","Celeus_ochraceus","Pygochelidon_melanoleuca","Cyanocorax_caeruleus","Atticora_tibialis","Hirundinea_bellicosa","Serpophaga_griseicapilla","Nannochordeiles_pusillus","Antrostomus_rufus","Myiothlypis_flaveola","Myiothlypis_leucoblephara","Spinus_yarrellii","Sturnella_superciliaris","Islerothraupis_cristata","Tangara_palmarum","Cantorchilus_leucotis","Cantorchilus_longirostris","Asemospiza_fuliginosa","Microspingus_cabanisi","Pheugopedius_genibarbis","Saltatricula_atricollis","Lipaugus_ater","Myrmoderus_ruficauda","Phalcoboenus_chimango","Psittacara_acuticaudatus","Sporophila_angolensis")

# Remove
short_birds_traits <- birds_traits[!birds_traits$sp %in% rem.birds.phy,]
nrow(short.birds.traits)
nrow(birds_traits)
View(short_birds_traits)

#Removendo nomes em duplicada
short.birds.traits <- short_birds_traits %>% distinct(sp, .keep_all = TRUE)
View(short.birds.traits)

# Mammals
#--- Checking phylogenis
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

# Transposition
mammals_trans <- t(mammals_traits)
colnames(mammals_trans) <- mammals_trans[1,]
mammals_trans <- mammals_trans[-1, ] #remove sp line 
#birds_trans <- birds_trans[,-34 ]#Aplastodiscus_weygoldti two times
#View(mammals_trans)

#--- Match with phylogeny
match.phylo.comm(mammals.upham.drop, mammals_trans)
mammals_phy <- prune.sample(mammals_trans, mammals.upham.drop)

#check names
#name <- as.data.frame(mammals.upham.drop$tip.label)
#View(name)

# remove names in traits
rem.mammals.phy <- c("Abrawayaomys_chebezi","Abrawayaomys_ruschii","Akodon_sanctipaulensis","Brucepattersonius_griserufescens","Brucepattersonius_guarani","Brucepattersonius_misionensis","Brucepattersonius_paradisus","Cabassous_tatouay","Callicebus_barbarabrownae","Callicebus_melanochir","Callithrix_flaviceps","Chacodelphys_formosa","Dasyprocta_iacki","Dasyprocta_prymnolopha","Dasypus_hybridus","Dasypus_septemcinctus","Graomys_chacoensis","Hylaeamys_oniscus","Leontopithecus_caissara","Leontopithecus_chrysopygus","Leopardus_geoffroyi","Leopardus_guttulus","Marmosa_paraguayana","Marmosops_bishopi","Monodelphis_unistriata","Nectomys_rattus","Oxymycterus_angularis","Oxymycterus_caparaoe","Oxymycterus_hispidus","Oxymycterus_roberti","Phaenomys_ferrugineus","Phyllomys_kerri","Phyllomys_medius","Phyllomys_thomasi","Phyllomys_unicolor","Reithrodon_typicus","Sapajus_cay","Sylvilagus_tapetillus","Tolypeutes_tricinctus","Trinomys_mirapitanga","Wilfredomys_oenax")

# Remove
short_mammals_traits <- mammals_traits[!mammals_traits$sp %in% rem.mammals.phy,]
nrow(short_mammals_traits)

# Fri Apr 29 17:06:10 2022 ------------------------------
### Preparing the data
# Anura
short_amphibia_traits <- as.data.frame(amphibia_traits_short)
short_amphibia_traits <- short_amphibia_traits %>% rename(species = sp) 
short_amphibia_traits$body_size <- log10(as.numeric(short_amphibia_traits$body_size))
short_amphibia_traits$litter_size <- log10(as.numeric(short_amphibia_traits$litter_size)) 

# Plot data before imputation 
phylo.heatmap(amphibia_phy_rooted,as.matrix(short_amphibia_traits),standardize=TRUE,lwd=3,pts=FALSE)
View(as.matrix(test))

#--- Imputing missing data - SVL/Body mass - litter size
# Preparing the phylogenetic data
anura_phy_ultra <- force.ultrametric(anura_phy)
amphibia_phy_rooted <- ape::multi2di(anura_phy_ultra)
is.ultrametric(amphibia_phy_rooted)

#--- Brownian motion
p_BM <- phylopars(trait_data = short_amphibia_traits, tree = amphibia_phy_rooted,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE)

# Variation checking
#p_BM$anc_recon[1:20,] # Data with imputed species means
#p_BM$anc_var[1:20,] # Variances for each estimate
#p_BM$anc_recon[1:20,] - sqrt(p_BM$anc_var[1:20,])*1.96 # Lower 95% CI
#$anc_recon[1:20,] + sqrt(p_BM$anc_var[1:20,])*1.96 # Upper 95% CI

# Checking with imputed data
# Traits amphibia
traits_amphy_input <- p_BM$anc_recon[1:562,]
traits_amphy_input <- as.data.frame(traits_amphy_input)

phylo.heatmap(amphibia_phy_rooted,as.matrix(traits_amphy_input),standardize=TRUE,lwd=3,pts=FALSE)

#--- OU motion
anura_phy_ultra <- force.ultrametric(anura_phy)
amphibia_phy_rooted <- ape::multi2di(anura_phy_ultra)
is.ultrametric(amphibia_phy_rooted)

p_OU <- phylopars(trait_data = test, tree = amphibia_phy_rooted,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE, model="mvOU")

amphibia_phy_rooted

p_OU$anc_recon[1:20,] # Data with imputed species means
p_OU$anc_var[1:20,] # Variances for each estimate

p_OU$anc_recon[1:20,] - sqrt(p_BM$anc_var[1:20,])*1.96 # Lower 95% CI
p_OU$anc_recon[1:20,] + sqrt(p_BM$anc_var[1:20,])*1.96 # Upper 95% CI
View(p_OU$anc_recon)

###
# Put species name in rows
#rownames(short_amphibia_traits) <- short_amphibia_traits$species
#short_amphibia_traits <- short_amphibia_traits[,-1]

# imputed traits amphibia
traits_amphy_input <- p_BM$anc_recon[1:562,]
traits_amphy_input <- as.data.frame(traits_amphy_input)

# Exploring relationship in body size and litter size
anura_traits_fig <- ggplot(traits_amphy_input, aes(x=svl, y=litter_size)) +
  geom_point() +
  labs(x = "SVL", y= "Litter size", title = "Relationship size and litter size") +
  stat_smooth(method = "loess") +
  theme_gray(base_size = 10)

p1 <- ggMarginal(anura_traits_fig, type="histogram")
p2 <- ggMarginal(anura_traits_fig, type="density")

###---------------------------------------
# Tue May 17 17:13:16 2022 ------------------------------
# Preparing the data of Squamata imputation
short_reptile_traits <- as.data.frame(reptile_traits_short)
short_reptile_traits <- short_reptile_traits %>% rename(species = sp) 

short_reptile_traits$body_size <- log10(as.numeric(short_reptile_traits$body_size))
short_reptile_traits$litter_size <- log10(as.numeric(short_reptile_traits$litter_size)) 

# Preparing the phylogenetic data
squamata_phy_ultra <- force.ultrametric(squamata_phy)
squamata_phy_rooted <- ape::multi2di(squamata_phy_ultra)
is.ultrametric(squamata_phy_rooted)

# Plot data before imputation 
phylo.heatmap(squamata_phy_rooted,as.matrix(short_reptile_traits),standardize=TRUE,lwd=3,pts=FALSE)

View(as.matrix(test))

#--- Imputing missing data - SVL/Body mass - litter size
#--- Brownian motion
p_BM_squamata <- phylopars(trait_data = short_reptile_traits, tree = squamata_phy_rooted,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE)
View(p_BM_squamata$anc_recon)

# Traits amphibia
traits_squamata_input <- as.data.frame(p_BM_squamata$anc_recon[1:405,])

# Checking with imputed data
phylo.heatmap(amphibia_phy_rooted,as.matrix(traits_amphy_input),standardize=TRUE,lwd=3,pts=FALSE)

# OU motion
p_OU_squamata <- phylopars(trait_data =short_reptile_traits, tree = squamata_phy_rooted,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE, model="mvOU")

traits_squamata_input_ou <- as.data.frame(p_OU_squamata$anc_recon[1:405,]) # Data with imputed species means

# Exploring relationship in body size and litter size
# brownian motion
rep_traits_fig <- ggplot(traits_squamata_input, aes(x=body_size, y=litter_size)) +
  geom_point() +
  labs(x = "SVL", y= "Litter size", title = "Relationship size and litter size") +
  stat_smooth(method = "loess") +
  theme_gray(base_size = 10)

p1 <- ggMarginal(rep_traits_fig, type="histogram")
p2 <- ggMarginal(rep_traits_fig, type="density")

# ou motion
rep_traits_fig_ou <- ggplot(traits_squamata_input_ou, aes(x=body_size, y=litter_size)) +
  geom_point() +
  labs(x = "SVL", y= "Litter size", title = "Relationship size and litter size") +
  stat_smooth(method = "loess") +
  theme_gray(base_size = 10)

p3 <- ggMarginal(rep_traits_fig_ou, type="histogram")
p4 <- ggMarginal(rep_traits_fig_ou, type="density")

# Birds
short_birds_traits <- as.data.frame(short.birds.traits)
short_birds_traits <- short_birds_traits %>% rename(species = sp) 

short_birds_traits$body_mass <- log10(as.numeric(short_birds_traits$body_mass))
short_birds_traits$litter_size <- log10(as.numeric(short_birds_traits$litter_size)) 

# Preparing the phylogenetic data
birds_phy_ultra <- force.ultrametric(birds_phy)
birds_phy_rooted <- ape::multi2di(birds_phy_ultra)
is.ultrametric(birds_phy_rooted)

# Plot data before imputation 
phylo.heatmap(birds_phy_rooted,as.matrix(short_birds_traits),standardize=TRUE,lwd=3,pts=FALSE)

View(short_birds_traits)

#--- Imputing missing data - SVL/Body mass - litter size
#--- Brownian motion
p_BM_birds <- phylopars(trait_data = short_birds_traits, tree = birds_phy_rooted,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE)

View(p_BM_birds$anc_recon)

# Traits amphibia
traits_birds_input <- as.data.frame(p_BM_birds$anc_recon[1:702,])

View(traits_birds_input)

# Checking with imputed data
phylo.heatmap(birds_phy_rooted,as.matrix(traits_birds_input),standardize=TRUE,lwd=3,pts=FALSE)

# OU motion
p_OU_birds <- phylopars(trait_data =short_birds_traits, tree = birds_phy_rooted,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE, model="mvOU")

traits_birds_input_ou <- as.data.frame(p_OU_birds$anc_recon[1:702,]) # Data with imputed species means

# Exploring relationship in body size and litter size
# brownian motion
birds_traits_fig <- ggplot(traits_birds_input, aes(x=body_mass, y=litter_size)) +
  geom_point() +
  labs(x = "Body mass", y= "Litter size", title = "Relationship mass and litter size") +
  stat_smooth(method = "loess") +
  theme_gray(base_size = 10)
?stat_smooth
p5 <- ggMarginal(birds_traits_fig, type="histogram")
p6 <- ggMarginal(birds_traits_fig, type="density")
p6

# ou motion
birds_traits_fig_ou <- ggplot(traits_birds_input_ou, aes(x=body_mass, y=litter_size)) +
  geom_point() +
  labs(x = "Body mass", y= "Litter size", title = "Relationship mass and litter size") +
  stat_smooth(method = "loess") +
  theme_gray(base_size = 10)

p7 <- ggMarginal(birds_traits_fig_ou, type="histogram")
p8 <- ggMarginal(birds_traits_fig_ou, type="density")
p8

View(traits_birds_input)

# Mammals
short.mammals.traits <- as.data.frame(short_mammals_traits)
short.mammals.traits <- short.mammals.traits %>% rename(species = sp) 

short.mammals.traits$body_mass <- log10(as.numeric(short.mammals.traits$body_mass))
short.mammals.traits$litter_size <- log10(as.numeric(short.mammals.traits$litter_size)) 

# Preparing the phylogenetic data
mammals_phy_ultra <- force.ultrametric(mammals_phy)
mammals_phy_rooted <- ape::multi2di(mammals_phy_ultra)
is.ultrametric(mammals_phy_rooted)

# Plot data before imputation 
phylo.heatmap(mammals_phy_rooted,as.matrix(short.mammals.traits),standardize=TRUE,lwd=3,pts=FALSE)

#--- Imputing missing data - SVL/Body mass - litter size
#--- Brownian motion
p_BM_mammals <- phylopars(trait_data = short.mammals.traits, tree = mammals_phy_rooted,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE)

View(p_BM_mammals$anc_recon)

# Traits amphibia
traits_mammals_input <- as.data.frame(p_BM_mammals$anc_recon[1:208,])

# Checking with imputed data
phylo.heatmap(mammals_phy_rooted,as.matrix(traits_mammals_input),standardize=TRUE,lwd=3,pts=FALSE)

# OU motion
p_OU_mammals <- phylopars(trait_data =short.mammals.traits, tree = mammals_phy_rooted,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE, model="mvOU")

traits_mammals_input_ou <- as.data.frame(p_OU_mammals$anc_recon[1:208,]) # Data with imputed species means

# Exploring relationship in body size and litter size
# brownian motion
mammals_traits_fig <- ggplot(traits_mammals_input, aes(x=body_mass, y=litter_size)) +
  geom_point() +
  labs(x = "Body mass", y= "Litter size", title = "Relationship mass and litter size") +
  stat_smooth(method = "loess") +
  theme_gray(base_size = 10)

p9 <- ggMarginal(mammals_traits_fig, type="histogram")
p10 <- ggMarginal(mammals_traits_fig, type="density")
p10

# ou motion
mammals_traits_fig_ou <- ggplot(traits_mammals_input_ou, aes(x=body_mass, y=litter_size)) +
  geom_point() +
  labs(x = "Body mass", y= "Litter size", title = "Relationship mass and litter size") +
  stat_smooth(method = "loess") +
  theme_gray(base_size = 10)

p11 <- ggMarginal(mammals_traits_fig_ou, type="histogram")
p12 <- ggMarginal(mammals_traits_fig_ou, type="density")
p11
p12

View(traits_birds_input)

#------------------------------------
#--- Calculing trait evolution rate
# Preparing the data
# Anura
svl_input <- as.numeric(traits_amphy_input$body_size)
litter_size_input <- as.numeric(traits_amphy_input$litter_size)

names(svl_input) <- rownames(traits_amphy_input)
names(litter_size_input) <- rownames(traits_amphy_input)

# Squamata
svl_input_squamata <- as.numeric(traits_squamata_input$body_size)
litter_size_input_squamata <- as.numeric(traits_squamata_input$litter_size)

names(svl_input_squamata) <- rownames(traits_squamata_input)
names(litter_size_input_squamata) <- rownames(traits_squamata_input)

# Birds
mass_input_birds <- as.numeric(traits_birds_input$body_mass)
litter_size_input_birds <- as.numeric(traits_birds_input$litter_size)

names(mass_input_birds) <- rownames(traits_birds_input)
names(litter_size_input_birds) <- rownames(traits_birds_input)

# Mammals
mass_input_mammals <- as.numeric(traits_mammals_input$body_mass)
litter_size_input_mammals <- as.numeric(traits_mammals_input$litter_size)

names(mass_input_mammals) <- rownames(traits_mammals_input)
names(litter_size_input_mammals) <- rownames(traits_mammals_input)

#--- Evolution rate with imputed data
# With data imputation
# Anura
rates.anura.svl <- RRphylo(tree= amphibia_phy_rooted, y=svl_input)
rates.anura.litter <- RRphylo(tree= amphibia_phy_rooted, y=litter_size_input)
# Visualizing Anura rates
rates_anura_svl <- rates.anura.svl$rates[546:1091,]
rates_anura_litter <- rates.anura.litter$rates[546:1091,]

# Squamata
rates.squa.svl <- RRphylo(tree= squamata_phy_rooted, y=svl_input_squamata)
rates.squa.litter <- RRphylo(tree= squamata_phy_rooted, y=litter_size_input_squamata)
# Visualizing squamata rates
rates_squa_svl <- rates.squa.svl$rates[405:809,]
rates_squa_litter <- rates.squa.litter$rates[405:809,]

# Birds
rates.bird.mass <- RRphylo(tree= birds_phy_rooted, y=mass_input_birds)
rates.bird.litter <- RRphylo(tree= birds_phy_rooted, y=litter_size_input_birds)
# Visualizing birds rates
rates_bird_svl <- rates.bird.mass$rates[405:809,]
rates_bird_litter <- rates.bird.litter$rates[405:809,]

# Mammals
rates.mammals.mass <- RRphylo(tree= mammals_phy, y= mass_input_mammals)
rates.mammals.litter <- RRphylo(tree=  mammals_phy, y=litter_size_input_mammals)

View(rates.mammals.mass$rates)

# Visualizing mammals rates
rates_mammals_svl <- rates.mammals.mass$rates[208:415,]
rates_mammals_litter <- rates.mammals.litter$rates[208:415,]

###---------------------------------------------------###
# Fri May 20 10:45:41 2022 ------------------------------
### Calculing trait evolution rate without imputation

#--- traits
short_amphibia_traits
vis_miss(short_amphibia_traits) # 64.64% NA in brood size ~ 3.72% body size
short_reptile_traits
vis_miss(short_reptile_traits) # 58.02% NA in brood size ~ 16.54% body size
short_birds_traits
vis_miss(short_birds_traits) # 26.5% NA in litter_size ~ 6.13% body_mass
short.mammals.traits
vis_miss(short.mammals.traits) # 0.96% NA in body_mass ~ 23.08%

#---Phylogenies
anura.phy 
reptile.phy 
birds.phy
mammals.phy

#--- Cut traits
litter_amph <- short_amphibia_traits[,-2]
litter_squa <- short_reptile_traits[,-2]
litter_bird <- short_birds_traits[,-2]
litter_mamm <- short.mammals.traits[,-2]

body_amph <- short_amphibia_traits[,-3]
body_squa <- short_reptile_traits[,-3]
body_bird <- short_birds_traits[,-3]
body_mamm <- short.mammals.traits[,-3]

# Removing missing data
litter_amph <- remove_missing(litter_amph, vars=names(litter_amph))
nrow(litter_amph) # 207 spp
litter_squa <- remove_missing(litter_squa, vars=names(litter_squa))
nrow(litter_squa) #170 spp
litter_bird <- remove_missing(litter_bird, vars=names(litter_bird))
nrow(litter_bird) # 516 spp
litter_mamm <- remove_missing(litter_mamm, vars=names(litter_mamm))
nrow(litter_mamm) # 160 spp

body_amph <- remove_missing(body_amph, vars=names(body_amph))
nrow(body_amph) # 546 spp
body_squa <- remove_missing(body_squa, vars=names(body_squa))
nrow(body_squa) # 338 spp
body_bird <- remove_missing(body_bird, vars=names(body_bird))
nrow(body_bird) # 659 spp
body_mamm <- remove_missing(body_mamm, vars=names(body_mamm))
nrow(body_mamm) # 206 spp

#--- Match with phylogeny
# Transposition
# Anura
body_amphi_trans <- t(body_amph)
colnames(body_amphi_trans) <- body_amph[,1]
#body_amphi_trans <- as.matrix(body_amphi_trans[-1, ]) #remove sp line 

litter_amph_trans <- as.matrix(t(litter_amph))
colnames(litter_amph_trans) <- litter_amph[,1]
#litter_mamm_trans_t <- litter_mamm_trans[-1, ]

# Squamata
# Birds

# Mammals
body_mamm_trans <- as.matrix(t(body_mamm))
colnames(body_mamm_trans) <- body_mamm[,1]
#body_mamm_trans <- body_mamm_trans[-1, ] #remove sp line 

litter_mamm_trans <- as.matrix(t(litter_mamm))
colnames(litter_mamm_trans) <- litter_mamm[,1]
#litter_mamm_trans_t <- litter_mamm_trans[-1, ]

#--- Match with phylogeny
#match.phylo.comm(mammals.upham.drop, body_mamm_trans)
# Anura
match.phylo.comm(anura_phy,litter_amph_trans)
body_amph_ <- body_amph[-34,] 

anura_phy_body <- prune.sample(body_amphi_trans, 
                                 anura_phy)
anura_phy_litter <- prune.sample(litter_amph_trans, 
                                 anura_phy)

#Squamata
#Birds

# Mammals
mammals_phy_body <- prune.sample(body_mamm_trans, 
                                 mammals.upham.drop)
mammals_phy_litter <- prune.sample(litter_mamm_trans, 
                                 mammals.upham.drop)

# Calculando taxa de evolução do atributo
#--- Calculing trait evolution rate
# Preparing the data
# Anura
svl_amph_sem <- as.numeric(body_amph_$body_size)
litter_size_amph_sem <- as.numeric(litter_amph$litter_size)

names(svl_amph_sem) <- body_amph_$species
names(litter_size_amph_sem) <- litter_amph$species

# Squamata
#svl_input_squamata <- as.numeric(traits_squamata_input$body_size)
#litter_size_input_squamata <- as.numeric(traits_squamata_input$litter_size)

#names(svl_input_squamata) <- rownames(traits_squamata_input)
#names(litter_size_input_squamata) <- rownames(traits_squamata_input)

# Birds
#mass_input_birds <- as.numeric(traits_birds_input$body_mass)
#litter_size_input_birds <- as.numeric(traits_birds_input$litter_size)

#names(mass_input_birds) <- rownames(traits_birds_input)
#names(litter_size_input_birds) <- rownames(traits_birds_input)

# Mammals
mass_mammals_sem <- as.numeric(body_mamm$body_mass)
litter_mammals_sem <- as.numeric(litter_mamm$litter_size)

names(mass_mammals_sem) <- body_mamm$species
names(litter_mammals_sem) <- litter_mamm$species

#--- Evolution rate with imputed data
# With data imputation
# Anura
svl_amph_sem 
litter_size_amph_sem

rates.anura.svl.sem <- RRphylo(tree= anura_phy_body, y=svl_amph_sem)
rates.anura.litter.sem <- RRphylo(tree= anura_phy_litter, y=litter_size_amph_sem)

View(rates.anura.svl.sem$rates)
# Creating objetc Anura rates
rates_anura_svl_sem <- rates.anura.svl.sem$rates[545:1089,]
rates_anura_litter <- rates.anura.litter.sem$rates[207:413,]

# Squamata
#rates.squa.svl <- RRphylo(tree= squamata_phy_rooted, y=svl_input_squamata)
#rates.squa.litter <- RRphylo(tree= squamata_phy_rooted, y=litter_size_input_squamata)
# Visualizing squamata rates
#rates_squa_svl <- rates.squa.svl$rates[405:809,]
#rates_squa_litter <- rates.squa.litter$rates[405:809,]

# Birds
#rates.bird.mass <- RRphylo(tree= birds_phy_rooted, y=mass_input_birds)
#rates.bird.litter <- RRphylo(tree= birds_phy_rooted, y=litter_size_input_birds)
# Visualizing birds rates
#rates_bird_svl <- rates.bird.mass$rates[405:809,]
#rates_bird_litter <- rates.bird.litter$rates[405:809,]

# Mammals
rates.mammals.mass.sem <- RRphylo(tree= mammals_phy_body, y= mass_mammals_sem)

rates.mammals.litter.sem <- RRphylo(tree=  mammals_phy_litter, y=litter_mammals_sem)
View(rates.mammals.mass.sem$rates)
rates.mammals.litter.sem$rates

# Visualizing mammals rates
rates_mammals_svl <- rates.mammals.mass.sem$rates[206:411,]
rates_mammals_litter <- rates.mammals.litter.sem$rates[160:319,]

###---------------------------------------------------###
# Fri May 20 10:47:12 2022 ------------------------------
# Extrair váriaveis climáticas - LetsR

###---------------------------------------------------###
# Fri May 20 10:47:12 2022 ------------------------------
# Modelos PGLS
