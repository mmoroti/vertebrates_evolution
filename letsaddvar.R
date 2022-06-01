# Wed Jun 01 15:35:55 2022 ------------------------------
# Objeto de presença ausência
#sf::st_crs(reptiles_ma)
#sf_read("modeled_reptiles.shp")
###
# Packages
library(letsR)

# Directory 
setwd("E:/OneDrive/2019 - Moroti/Dados/Shape/Répteis")

#Shapefiles
reptiles_ma <- readOGR("modeled_reptiles.shp")
reptiles_ma <- sf::st_transform(teste, "+proj=longlat +datum=WGS84")

# Projection
projection(temp) <- projection(reptiles_ma)
projection(reptiles_ma)
projection(temp)

# Creating PresenceAusence object for lets.addvar
teste <- lets.presab(reptiles_ma , xmn = -180, xmx = 180, ymn = -60, ymx = 90)
# lets.addvar(PresenceAbsence object, raster variables, onlyvar = If TRUE only the matrix object will be returned., fun = mean)
teste_temp_median <- lets.addvar(teste, temp, fun = median, onlyvar = FALSE)

#--- Anura
# Directory
setwd("E:/OneDrive/2019 - Moroti/Dados/Shape/Amphibia")
dir()
anura_ma <- readOGR("ANURA_SA_clip.shp")
#anura_ma <- sf::st_transform(teste, "+proj=longlat +datum=WGS84")

# Projection
projection(temp) <- projection(anura_ma)
projection(anura_ma)
projection(temp)

# Creating PresenceAusence object for lets.addvar
anura_presaus <- lets.presab(anura_ma , xmn = -180, xmx = 180, ymn = -60, ymx = 90)
# lets.addvar(PresenceAbsence object, raster variables, onlyvar = If TRUE only the matrix object will be returned., fun = mean)
anura_temp_median <- lets.addvar(anura_presaus, temp, fun = median, onlyvar = FALSE)
View(anura_temp_median)

