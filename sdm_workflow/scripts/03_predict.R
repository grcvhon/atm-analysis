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
  
  # Mean and variance in spatial prediction
  prediction_mean <- prediction$mean
  prediction_sd <- prediction$sd

  # thresholded models
  TSS_me <- prediction_mean > results_list$leaf_scaled_opt_mod.rds$model_evaluation_metrics$TSS
  LPT_me <- prediction_mean > results_list$leaf_scaled_opt_mod.rds$model_evaluation_metrics$LPT
  thresh_me <- raster::cut(raster(prediction_mean), breaks = c(-Inf, 0.25, 0.5, 0.76, Inf))
  
  # Your binary raster
  r <- TSS_me  # or thresh if that’s the name
  
  # Count how many TRUE (or 1) cells
  n_cells <- global(r, "sum", na.rm = TRUE)
  n_cells
  
  # Get cell area in km²
  cell_area <- cellSize(r, unit = "km")
  
  # Multiply cell area by suitability mask
  suitable_area <- mask(cell_area, r, maskvalues = FALSE)
  
  # Sum total area
  total_area <- global(suitable_area, "sum", na.rm = TRUE)
  total_area

  predict_plot <- ggplot() + 
    geom_spatraster(data = prediction$mean) +
    scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
    labs(title = paste0("Mean spatial prediction for Leaf-scaled sea snake (", last_part, ")")) +
    annotation_scale(mapping = aes(location = "br")) +
    theme_bw()

  thresh_me_df <- as.data.frame(thresh_me, xy = TRUE, na.rm = TRUE)

  thresh_me_plot <- ggplot(thresh_me_df, aes(x = x, y = y, fill = layer)) +
    geom_tile() +
    scale_fill_distiller(palette = "Spectral", na.value = "transparent") + 
    labs(title = paste0("Threshold prediction for Leaf-scaled sea snake (", last_part, ")")) +
    annotation_scale(mapping = aes(location = "br")) +
    coord_fixed() +
    theme_bw()

  leaf_scaled_predict_list[[last_part]] <- 
    list(
      prediction_mean = prediction_mean,
      prediction_sd = prediction_sd, 
      TSS_me = TSS_me,
      LPT_me = LPT_me,
      threshold_me = thresh_me,
      binary_raster = r,
      number_TRUE_cells = n_cells,
      cell_area_km2 = cell_area,
      suitable_area = suitable_area,
      total_area = total_area,
      mean_prediction_plot = predict_plot,
      threshold_plot = thresh_me_plot
    )
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
  
  # Mean and variance in spatial prediction
  prediction_mean <- prediction$mean
  prediction_sd <- prediction$sd

  # thresholded models
  TSS_me <- prediction_mean > results_list$short_nosed_opt_mod.rds$model_evaluation_metrics$TSS
  LPT_me <- prediction_mean > results_list$short_nosed_opt_mod.rds$model_evaluation_metrics$LPT
  thresh_me <- raster::cut(raster(prediction_mean), breaks = c(-Inf, 0.25, 0.5, 0.75, Inf))
  
  # Your binary raster
  r <- TSS_me  # or thresh if that’s the name
  
  # Count how many TRUE (or 1) cells
  n_cells <- global(r, "sum", na.rm = TRUE)
  n_cells
  
  # Get cell area in km²
  cell_area <- cellSize(r, unit = "km")
  
  # Multiply cell area by suitability mask
  suitable_area <- mask(cell_area, r, maskvalues = FALSE)
  
  # Sum total area
  total_area <- global(suitable_area, "sum", na.rm = TRUE)
  total_area

  predict_plot <- ggplot() + 
    geom_spatraster(data = prediction$mean) +
    scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
    labs(title = paste0("Mean spatial prediction for Short-nosed sea snake (", last_part, ")")) +
    annotation_scale(mapping = aes(location = "br")) +
    theme_bw()

  thresh_me_df <- as.data.frame(thresh_me, xy = TRUE, na.rm = TRUE)

  thresh_me_plot <- ggplot(thresh_me_df, aes(x = x, y = y, fill = layer)) +
    geom_tile() +
    scale_fill_distiller(palette = "Spectral", na.value = "transparent") + 
    labs(title = paste0("Threshold prediction for Short-nosed sea snake (", last_part, ")")) +
    annotation_scale(mapping = aes(location = "br")) +
    coord_fixed() +
    theme_bw()

  short_nosed_predict_list[[last_part]] <- 
    list(
      prediction_mean = prediction_mean,
      prediction_sd = prediction_sd, 
      TSS_me = TSS_me,
      LPT_me = LPT_me,
      threshold_me = thresh_me,
      binary_raster = r,
      number_TRUE_cells = n_cells,
      cell_area_km2 = cell_area,
      suitable_area = suitable_area,
      total_area = total_area,
      mean_prediction_plot = predict_plot,
      threshold_plot = thresh_me_plot
    )
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