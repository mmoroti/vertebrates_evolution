# Wed Jun 01 15:35:55 2022 ------------------------------
# Objeto de presença ausência
#sf::st_crs(reptiles_ma)
#sf_read("modeled_reptiles.shp")
###
# Packages
library(letsR)
library(raster)

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

# Climate data
r <- raster(getData("worldclim", var = "bio", 
                    res = 10), 1) / 10

#--- Anura
# Directory
setwd("E:/OneDrive/2019 - Moroti/Dados/Shape/Amphibia")
dir()
anura_ma <- readOGR("ANURA_SA_clip.shp")
#anura_ma <- sf::st_transform(teste, "+proj=longlat +datum=WGS84")

# Projection
projection(r) <- projection(anura_ma)
projection(anura_ma)

# Creating PresenceAusence object for lets.addvar
anura_presaus <- lets.presab(anura_ma , xmn = -180, xmx = 180, ymn = -60, ymx = 90)

projection(anura_presaus$Richness_Raster) <- projection(r)

# lets.addvar(PresenceAbsence object, raster variables, onlyvar = If TRUE only the matrix object will be returned., fun = mean)
anura_temp_median <- lets.addvar(anura_presaus, r, fun = median, onlyvar = FALSE)
View(anura_temp_median)
# Now, we need to summarise the climate variables per specie
temp_median_anura <- lets.summarizer2(anura_temp_median,  pos = ncol(anura_temp_median), fun = median) 

temp_median_anura
hist(as.numeric(temp_median_anura$bio1_median))

###
#--- Squamata
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

#--- Birds
# Directory
setwd("E:/OneDrive/2019 - Moroti/Dados/Shape/Aves")
dir()
birds_ma <- readOGR("birds.shp")
#anura_ma <- sf::st_transform(teste, "+proj=longlat +datum=WGS84")
names(birds_ma)
# Projection
projection(temp) <- projection(birds_ma)
projection(birds_ma)
projection(temp)

# Creating PresenceAusence object for lets.addvar
birds_presaus <- lets.presab(birds_ma , xmn = -180, xmx = 180, ymn = -60, ymx = 90)

# lets.addvar(PresenceAbsence object, raster variables, onlyvar = If TRUE only the matrix object will be returned., fun = mean)
birds_temp_median <- lets.addvar(birds_presaus, temp, fun = median, onlyvar = FALSE)
View(birds_temp_median)

#--- Mammals
# Directory
setwd("E:/OneDrive/2019 - Moroti/Dados/Shape/Mamíferos")
dir()
mammals_ma <- readOGR("data_0.shp")
#anura_ma <- sf::st_transform(teste, "+proj=longlat +datum=WGS84")

# Projection
projection(temp) <- projection(mammals_ma)
projection(mammals_ma)
projection(temp)

# Creating PresenceAusence object for lets.addvar
mammals_presaus <- lets.presab(mammals_ma , xmn = -180, xmx = 180, ymn = -60, ymx = 90)
# lets.addvar(PresenceAbsence object, raster variables, onlyvar = If TRUE only the matrix object will be returned., fun = mean)
mammals_temp_median <- lets.addvar(mammals_presaus, temp, fun = median, onlyvar = FALSE)
View(mammals_temp_median)