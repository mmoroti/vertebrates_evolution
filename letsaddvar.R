library(geiger)

?fitContinuous()

## Not run:
data(temp) # Global mean temperature
data(PAM) # Phyllomedusa presence-absence matrix
# Mean temperature
PAM_temp_mean <- lets.addvar(PAM, temp)
# Standard deviation of temperature
PAM_temp_sd <- lets.addvar(PAM, temp, fun = sd, onlyvar = TRUE)
# Mean and SD in the PAM
PAM_temp_mean_sd <- cbind(PAM_temp_mean, PAM_temp_sd)
## End(Not run)


class(PAM_temp_mean_sd) # matrix array
class(temp) # rasterLayer
View(PAM$Presence_and_Absence_Matrix) # matrix array

