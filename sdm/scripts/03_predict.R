# sdm - 03 predict

## setwd: sdm directory
#setwd("C:/Users/a1235304/Dropbox/Short-nosed and Leaf-scaled sea snake TSSC/atm-analysis/sdm")
#
## install/load packages
#packages <- c("sf","leaflet","readr","janitor","dplyr",
#              "mapview","spatstat","tidyverse","raster",
#              "dismo","lubridate","SDMtune","readxl",
#              "terra","stars","lwgeom","maptools","ggspatial",
#              "prettymapr","tidyterra","caret","corrplot",
#              "plotROC","ggpubr")
#
## load package list all at once
#invisible(lapply(packages, library, character.only = TRUE))

##############################################################################################################

### predict ###

# data for prediction

### predict coastal
env_predict_coastal <- env_in
# to isolate coastal, make ancestry coefficient values of less than 0.0 (=offshore) a value of 0
env_predict_coastal$Ancestry.coefficient[values(env_predict_coastal$Ancestry.coefficient) < 0.0] <- 0
plot(env_predict_coastal)

### predict offshore
env_predict_offshore <- env_in
# to isolate offshore, make ancestry coefficient values of greater than 0.0 (=coastal) a value of 0
env_predict_offshore$Ancestry.coefficient[values(env_predict_offshore$Ancestry.coefficient) > 0.0] <- 0
plot(env_predict_offshore)

# ancestry layer for prediction
env_predict_coastal
env_predict_offshore

# optimised model for prediction
results_list$short_nosed_opt_mod.rds$best_model
results_list$leaf_scaled_opt_mod.rds$best_model

# put in a species-specific loop

# create a list of data for prediction
env_predict_list <- list(
  coastal = env_predict_coastal,
  offshore = env_predict_offshore
)

# predict: Leaf-scaled sea snake - generates two prediction maps
leaf_scaled_predict_list <- list()

for (location in names(env_predict_list)) {
  prediction <- 
    SDMtune::predict(
      results_list$leaf_scaled_opt_mod.rds$best_model, 
      data = env_predict_list[[location]],
      fun = c("mean", "sd"), 
      type = "cloglog", parallel = T)
  
  string <- location
  last_part <- sub(".*_", "", string)
  
  predict_plot <- ggplot() + 
    geom_spatraster(data = prediction$mean) +
    scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
    labs(title = paste0("Mean spatial prediction for Leaf-scaled sea snake (", last_part, ")")) +
    theme_bw()

   leaf_scaled_predict_list[[last_part]] <- predict_plot
}

# predict: Short-nosed sea snake - generates two prediction maps
short_nosed_predict_list <- list()

for (location in names(env_predict_list)) {
  prediction <- 
    SDMtune::predict(
      results_list$short_nosed_opt_mod.rds$best_model, 
      data = env_predict_list[[location]],
      fun = c("mean", "sd"), 
      type = "cloglog", parallel = T)
  
  string <- location
  last_part <- sub(".*_", "", string)
  
  predict_plot <- ggplot() + 
    geom_spatraster(data = prediction$mean) +
    scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
    labs(title = paste0("Mean spatial prediction for Short-nosed sea snake (", last_part, ")")) +
    theme_bw()

   short_nosed_predict_list[[last_part]] <- predict_plot
}






## predict step
#predict <- 
#  SDMtune::predict(
#    results_list$leaf_scaled_opt_mod.rds$best_model, 
#    data = env_predict_coastal,
#    fun = c("mean", "sd"), 
#    type = "cloglog", parallel = T)
#
#ggplot() + 
#  geom_spatraster(data = predict$mean) +
#  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
#  theme_bw()