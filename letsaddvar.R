# Wed Jun 01 15:35:55 2022 ------------------------------
# Herein, we extracted the median of the climate niche of species range.
# Packages
library(letsR)
library(raster)
library(rgdal)
library(tidyverse)

# Function
# The function lets.summarise returns NA's in climate values. We remove modified the funcion. and remove the last 'for'
lets.summarizer2 <- function (x, pos, xy = TRUE, fun = mean) 
{
  var <- x[, pos, drop = FALSE]
  sp <- x[, -pos, drop = FALSE]
  if (xy) {
    sp <- sp[, -(1:2), drop = FALSE]
  }
  Species <- colnames(sp)
  n <- length(Species)
  lpos <- length(pos)
  resum <- matrix(NA, nrow = n, ncol = lpos)
  colnames(resum) <- colnames(var)
  for (i in 1:n) {
    vari <- var[(sp[, i] == 1), , drop = FALSE]
    is_all_na <- apply(vari, 2, function(x) {
      all(is.na(x))
    })
    if (nrow(vari) == 0 | all(is_all_na)) {
      resum[i, ] <- rep(NA, lpos)
    }
    else {
      if (any(is_all_na)) {
        resum[i, is_all_na] <- NA
      }
      resum[i, !is_all_na] <- apply(vari[, !is_all_na, 
                                         drop = FALSE], 2, fun, na.rm = TRUE)
    }
  }
  resul <- as.data.frame(cbind(Species, resum))
  return(resul)
}

# Climate data from CHELSA
Bio02 <- raster("C:/Users/Cliente/Desktop/vertebrates_evolution/CHELSA_bio2_1981-2010_V.2.1.tif") #Mean Diurnal Temperature Range
Bio03 <- raster("C:/Users/Cliente/Desktop/vertebrates_evolution/CHELSA_bio3_1981-2010_V.2.1.tif") #Isothermality
Bio04 <- raster("C:/Users/Cliente/Desktop/vertebrates_evolution/CHELSA_bio4_1981-2010_V.2.1.tif") #Temperature Seasonality
Bio07 <- raster("C:/Users/Cliente/Desktop/vertebrates_evolution/CHELSA_bio7_1981-2010_V.2.1.tif") #Temperature Annual Range
Bio11 <- raster("C:/Users/Cliente/Desktop/vertebrates_evolution/CHELSA_bio11_1981-2010_V.2.1.tif") #Mean Temperature of Coldest Quarter
Bio15 <- raster("C:/Users/Cliente/Desktop/vertebrates_evolution/CHELSA_bio15_1981-2010_V.2.1.tif") #Precipitation Seasonality
Bio17 <- raster("C:/Users/Cliente/Desktop/vertebrates_evolution/CHELSA_bio17_1981-2010_V.2.1.tif") #Precipitation of Driest Quarter

#--- Anura
# Directory
setwd("E:/OneDrive/2019 - Moroti/Dados/Shape/Amphibia")
dir()
anura_ma <- readOGR("ANURA_SA_clip.shp")
#anura_ma <- sf::st_transform(teste, "+proj=longlat +datum=WGS84")

# Projection
#projection(Bio02) <- projection(anura_ma)

# Creating PresenceAusence object for lets.addvar
anura_presaus <- lets.presab(anura_ma , xmn = -180, xmx = 180, ymn = -60, ymx = 90)

#projection(anura_presaus$Richness_Raster) <- projection(Bio02)

# lets.addvar(PresenceAbsence object, raster variables, onlyvar = If TRUE only the matrix object will be returned., fun = mean)
anura_temp_median_02 <- lets.addvar(anura_presaus, Bio02, fun = median, onlyvar = FALSE)
anura_temp_median_03 <- lets.addvar(anura_presaus, Bio03, fun = median, onlyvar = TRUE)
anura_temp_median_04 <- lets.addvar(anura_presaus, Bio04, fun = median, onlyvar = TRUE) 
anura_temp_median_07 <- lets.addvar(anura_presaus, Bio07, fun = median, onlyvar = TRUE)
anura_temp_median_11 <- lets.addvar(anura_presaus, Bio11, fun = median, onlyvar = TRUE)
anura_temp_median_15 <- lets.addvar(anura_presaus, Bio15, fun = median, onlyvar = TRUE)
anura_temp_median_17 <- lets.addvar(anura_presaus, Bio17, fun = median, onlyvar = TRUE)

# Grouping data set
anura_temp_median <- cbind(anura_temp_median_02, anura_temp_median_03, anura_temp_median_04, anura_temp_median_07, anura_temp_median_11, anura_temp_median_15, anura_temp_median_17)

# Now, we need to summarise the climate variables per specie
temp_median_anura <- lets.summarizer2(anura_temp_median,  pos = 2625:2631, fun = median) 
glimpse(temp_median_anura)

# We need to adjust some names in this dataset, and change climate variables to 'numeric' 
# Preparing the data for Phylogenetic Generalized Last Squared
temp_median_anura$Species <- temp_median_anura$Species %>% 
  str_replace_all( "Brachycephalus_margariatus", "Brachycephalus_margaritatus") %>% 
  str_replace_all("Agalychnis_granulosa", "Hylomantis_granulosa") %>%
  str_replace_all("Agalychnis_aspera", "Hylomantis_aspera") %>% 
  str_replace_all("Pithecopus_ayeaye", "Phyllomedusa_ayeaye")%>%
  str_replace_all("Pithecopus_megacephalus", "Phyllomedusa_megacephala")%>%
  str_replace_all("Pithecopus_rohdei", "Phyllomedusa_rohdei")%>%
  str_replace_all("Pithecopus_azureus", "Phyllomedusa_azurea")%>%
  str_replace_all("Pithecopus_hypochondrialis", "Phyllomedusa_hypochondrialis") %>%
  str_replace_all("Pithecopus_nordestinus", "Phyllomedusa_nordestina")%>%
  str_replace_all("Boana_paranaiba", "Hypsiboas_paranaiba")%>%
  str_replace_all("Boana_poaju", "Hypsiboas_poaju") %>%
  str_replace_all("Boana_bandeirantes", "Hypsiboas_bandeirantes")%>%
  str_replace_all("Ololygon_tripui", "Scinax_tripui")%>%
  str_replace_all("Ololygon_cosenzai", "Scinax_cosenzai")%>%
  str_replace_all("Aparasphenodon_ararapa", "Aparasphenodon_arapapa")%>%
  str_replace_all("Physalaemus_nattereri", "Eupemphix_nattereri") %>%
  str_replace_all("Lithobates_palmipes", "Rana_palmipes")

# Grouping data
litter_anura <- as.data.frame(rates_anura_litter)
body_anura <- as.data.frame(rates_anura_svl)

sp_anura_rates <- cbind(Species = rownames(body_anura), body_anura, litter_anura)

anura_dataset <- left_join(sp_anura_rates,temp_median_anura,  by="Species")
glimpse(anura_dataset)

#--- Squamata
# Directory 
setwd("C:/Users/Cliente/Desktop/vertebrates_evolution/shapefiles_ranges/repteis")

#Shapefiles
reptiles_ma <- readOGR("modeled_reptiles.shp")
#reptiles_ma <- sf::st_transform(teste, "+proj=longlat +datum=WGS84")

# Projection
projection(reptiles_ma)

# Creating PresenceAusence object for lets.addvar
repteis_presaus <- lets.presab(reptiles_ma , xmn = -180, xmx = 180, ymn = -60, ymx = 90)

# lets.addvar(PresenceAbsence object, raster variables, onlyvar = If TRUE only the matrix object will be returned., fun = mean)
squa_temp_median_02 <- lets.addvar(repteis_presaus, Bio02, fun = median, onlyvar = FALSE)
squa_temp_median_03 <- lets.addvar(repteis_presaus, Bio03, fun = median, onlyvar = TRUE)
squa_temp_median_04 <- lets.addvar(repteis_presaus, Bio04, fun = median, onlyvar = TRUE) 
squa_temp_median_07 <- lets.addvar(repteis_presaus, Bio07, fun = median, onlyvar = TRUE)
squa_temp_median_11 <- lets.addvar(repteis_presaus, Bio11, fun = median, onlyvar = TRUE)
squa_temp_median_15 <- lets.addvar(repteis_presaus, Bio15, fun = median, onlyvar = TRUE)
squa_temp_median_17 <- lets.addvar(repteis_presaus, Bio17, fun = median, onlyvar = TRUE)

# Grouping data set
squa_temp_median <- cbind(squa_temp_median_02, squa_temp_median_03, squa_temp_median_04, squa_temp_median_07, squa_temp_median_11, squa_temp_median_15,squa_temp_median_17)
glimpse(squa_temp_median)

# Now, we need to summaries the climate variables per specie
ncol(squa_temp_median) 
temp_median_squa <- lets.summarizer2(squa_temp_median,  pos = 1049:1055, fun = median) 
glimpse(temp_median_squa)

# Change '//s' for '_'
temp_median_squa$Species <- str_replace_all(temp_median_squa$Species, '\\s', '_')

# Grouping data
litter_squa <- as.data.frame(rates_squa_litter)
body_squa <- as.data.frame(rates_squa_svl)

sp_squa_rates <- cbind(Species = rownames(body_squa), 
                       body_squa, litter_squa)

squa_dataset <- left_join(sp_squa_rates,temp_median_squa,  by="Species")

nrow(sp_squa_rates) #405 species
nrow(squa_dataset)  #405 species
glimpse(squa_dataset)

#--- Birds
# Directory
setwd("C:/Users/Cliente/Desktop/vertebrates_evolution/shapefiles_ranges/aves")
birds_ma <- readOGR("birds.shp")

#anura_ma <- sf::st_transform(teste, "+proj=longlat +datum=WGS84")
#names(birds_ma)
# Projection
#projection(temp) <- projection(birds_ma)
projection(birds_ma)
#projection(temp)

# Creating PresenceAusence object for lets.addvar
birds_presaus <- lets.presab(birds_ma, xmn = -180, xmx = 180, ymn = -60, ymx = 90)

# lets.addvar(PresenceAbsence object, raster variables, onlyvar = If TRUE only the matrix object will be returned., fun = mean)
#birds_temp_median <- lets.addvar(birds_presaus, temp, fun = median, onlyvar = FALSE)
birds_temp_median_02 <- lets.addvar(birds_presaus, Bio02, fun = median, onlyvar = FALSE)
birds_temp_median_03 <- lets.addvar(birds_presaus, Bio03, fun = median, onlyvar = TRUE)
birds_temp_median_04 <- lets.addvar(birds_presaus, Bio04, fun = median, onlyvar = TRUE) 
birds_temp_median_07 <- lets.addvar(birds_presaus, Bio07, fun = median, onlyvar = TRUE)
birds_temp_median_11 <- lets.addvar(birds_presaus, Bio11, fun = median, onlyvar = TRUE)
birds_temp_median_15 <- lets.addvar(birds_presaus, Bio15, fun = median, onlyvar = TRUE)
birds_temp_median_17 <- lets.addvar(birds_presaus, Bio17, fun = median, onlyvar = TRUE)

# Grouping data set
birds_temp_median <- cbind(birds_temp_median_02, birds_temp_median_03, birds_temp_median_04, birds_temp_median_07, birds_temp_median_11, birds_temp_median_15, birds_temp_median_17)
glimpse(birds_temp_median)

# Now, we need to summaries the climate variables per specie
ncol(birds_temp_median) 
temp_median_birds <- lets.summarizer2(birds_temp_median,  pos = 1526:1532, fun = median) 
View(temp_median_birds)

# Change '//s' for '_'
temp_median_birds$Species <- str_replace_all(temp_median_birds$Species, '\\s', '_')

# Grouping data
litter_bird <- as.data.frame(rates_bird_litter)
body_bird <- as.data.frame(rates_bird_svl)

sp_bird_rates <- cbind(Species = rownames(body_bird), 
                       body_bird, litter_bird)

bird_dataset <- left_join(sp_bird_rates,temp_median_birds,  by="Species")

nrow(sp_bird_rates) #701 species
nrow(bird_dataset)  #701 species
glimpse(squa_dataset)

#--- Mammals
# Directory
setwd("C:/Users/Cliente/Desktop/vertebrates_evolution/shapefiles_ranges/mamiferos")
mammals_ma <- readOGR("data_0.shp")
#anura_ma <- sf::st_transform(teste, "+proj=longlat +datum=WGS84")

# Projection
#projection(temp) <- projection(mammals_ma)
#projection(mammals_ma)
#projection(temp)

# Creating PresenceAusence object for lets.addvar
mammals_presaus <- lets.presab(mammals_ma , xmn = -180, xmx = 180, ymn = -60, ymx = 90)

# lets.addvar(PresenceAbsence object, raster variables, onlyvar = If TRUE only the matrix object will be returned., fun = mean)
#birds_temp_median <- lets.addvar(birds_presaus, temp, fun = median, onlyvar = FALSE)
mammals_temp_median_02 <- lets.addvar(mammals_presaus, Bio02, fun = median, onlyvar = FALSE)
mammals_temp_median_03 <- lets.addvar(mammals_presaus, Bio03, fun = median, onlyvar = TRUE)
mammals_temp_median_04 <- lets.addvar(mammals_presaus, Bio04, fun = median, onlyvar = TRUE) 
mammals_temp_median_07 <- lets.addvar(mammals_presaus, Bio07, fun = median, onlyvar = TRUE)
mammals_temp_median_11 <- lets.addvar(mammals_presaus, Bio11, fun = median, onlyvar = TRUE)
mammals_temp_median_15 <- lets.addvar(mammals_presaus, Bio15, fun = median, onlyvar = TRUE)
mammals_temp_median_17 <- lets.addvar(mammals_presaus, Bio17, fun = median, onlyvar = TRUE)

# Grouping data set
mammals_temp_median <- cbind(mammals_temp_median_02, mammals_temp_median_03, mammals_temp_median_04, mammals_temp_median_07, mammals_temp_median_11, mammals_temp_median_15, mammals_temp_median_17)
glimpse(mammals_temp_median)

# Now, we need to summaries the climate variables per specie
ncol(mammals_temp_median) 
temp_median_mammals <- lets.summarizer2(mammals_temp_median,  pos = 1255:1261, fun = median) 
glimpse(temp_median_mammals)

# Change '//s' for '_'
temp_median_mammals$Species <- str_replace_all(temp_median_mammals$Species, '\\s', '_')

# Grouping data
body_mammal <- as.data.frame(rates_mammals_svl)
litter_mammal <- as.data.frame(rates_mammals_litter)

sp_mammal_rates <- cbind(Species = rownames(body_mammal), 
                         body_mammal, litter_mammal)

mammals_dataset <- left_join(sp_mammal_rates,temp_median_mammals,  by="Species")

nrow(sp_mammal_rates)  #208 species
nrow(mammals_dataset)  #208 species
glimpse(mammals_dataset)

###
# Salvando os dataset com nicho climático e taxa de evolução das espécies
getwd()
write.csv(mammals_dataset, "mammals_dataset.csv", row.names = FALSE)
write.csv(bird_dataset, "birds_dataset.csv", row.names = FALSE)
write.csv(squa_dataset, "squamata_dataset.csv", row.names = FALSE)
write.csv(anura_dataset, "anura_dataset.csv", row.names = FALSE)

read.csv2("anura_dataset.csv", sep=",")
read.csv2("squamata_dataset.csv", sep=",")
read.csv2("birds_dataset.csv", sep=",")
read.csv2("mammals_dataset.csv", sep=",")

###
# Testando GLS
names(anura_dataset)
library(nlme)
library(visdat)
vis_miss(mammals_dataset)

pgls<- gls(rates_anura_svl ~ CHELSA_bio2_1981.2010_V.2.1_median, data=anura_dataset)

summary(pgls)
