# Mon May 30 11:31:00 2022 ------------------------------

# I modified the function fitGeospiza so that it receives the two arguments (phy, trait) that will go in the model. The function's new name is fit_modified

#---- CREATE a FUNCTION to COMPARE EVOLUTION MODELS
fit_modified=function(phy, trait){
  #trait=match.arg(trait, c("body_size","litter_size"))
  
  # define set of models to compare
  models=c("BM", "OU", "EB", "white")
  summaries=c("diffusion", "Ornstein-Uhlenbeck", "early burst", "white noise")
  
  ## ESTIMATING measurement error ##
  aic.se=numeric(length(models))
  lnl.se=numeric(length(models))
  
  for(m in 1:length(models)){
    cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep=""), models[m],
        " with SE *** \n", sep="")
    tmp=fitContinuous(phy, trait,SE=NA, model=models[m],
                      bounds=list(SE=c(0,0.5)), ncores=2)
    print(tmp)
    aic.se[m]=tmp$opt$aicc
    lnl.se[m]=tmp$opt$lnL
  }
  
  
  ## ASSUMING no measurement error ##
  aic=numeric(length(models))
  lnl=numeric(length(models))
  
  for(m in 1:length(models)){
    cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep=""), models[m],
        " *** \n", sep="")
    tmp=fitContinuous(phy, trait,SE=0,model=models[m], ncores=2)
    print(tmp)
    aic[m]=tmp$opt$aicc
    lnl[m]=tmp$opt$lnL
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
  
  res_aicc=rbind(aic, aic.se, daic, daic.se)
  rownames(res_aicc)=c("AICc","AICc_SE","dAICc", "dAICc_SE")
  
  return(res_aicc)
}

###--- Anura
# The Phylogeny needs to rooted and ultrametric
anura_phy_ultra_2 <- force.ultrametric(anura_phy_body)
amphibia_phy_rooted_2 <- ape::multi2di(anura_phy_ultra_2)

# Body models
body_anura_model <- fit_modified(amphibia_phy_rooted_2, svl_amph_sem)
body_anura_model # OU model

# Ultrametric and rooted
anura_phy_ultra_3 <- force.ultrametric(anura_phy_litter)
amphibia_phy_rooted_3 <- ape::multi2di(anura_phy_ultra_3)

# Litter size models
litter_anura_model <- fit_modified(amphibia_phy_rooted_3, litter_size_amph_sem)
litter_anura_model # OU model

###--- Squamata
# The Phylogeny needs to rooted and ultrametric
squa_phy_ultra_2 <- force.ultrametric(squa_phy_body)
squa_phy_rooted_2 <- ape::multi2di(squa_phy_ultra_2)

# Body size Model 
body_squa_model <- fit_modified(squa_phy_rooted_2, svl_squamata)
body_squa_model # OU model

# Ultrametric and rooted
squa_phy_ultra_3 <- force.ultrametric(squa_phy_litter)
squa_phy_rooted_3 <- ape::multi2di(squa_phy_ultra_3)

# Litter size model
litter_squa_model <- fit_modified(squa_phy_rooted_3, litter_size_squamata)
litter_squa_model # OU model

###--- Birds
# The Phylogeny needs to rooted and ultrametric
bird_phy_ultra_2 <- force.ultrametric(bird_phy_body)
bird_phy_rooted_2 <- ape::multi2di(bird_phy_ultra_2)

# Body size Model 
body_bird_model <- fit_modified(bird_phy_rooted_2, mass_birds)
body_bird_model # EB model

# Ultrametric and rooted
bird_phy_ultra_3 <- force.ultrametric(bird_phy_litter)
bird_phy_rooted_3 <- ape::multi2di(bird_phy_ultra_3)

# Litter size model
litter_bird_model <- fit_modified(bird_phy_rooted_3, litter_size_birds)
litter_bird_model # OU model

###--- Mammals
# The Phylogeny needs to rooted and ultrametric
mam_phy_ultra_2 <- force.ultrametric(mammals_phy_body)
mam_phy_rooted_2 <- ape::multi2di(mam_phy_ultra_2)

# Body size Model 
body_mam_model <- fit_modified(mam_phy_rooted_2, mass_mammals_sem)
body_mam_model # OU model

# Ultrametric and rooted
mam_phy_ultra_3 <- force.ultrametric(mammals_phy_litter)
mam_phy_rooted_3 <- ape::multi2di(mam_phy_ultra_3)

# Litter size model
litter_mam_model <- fit_modified(mam_phy_rooted_3, litter_mammals_sem)
litter_mam_model # OU model
