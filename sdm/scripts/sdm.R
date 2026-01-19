# Species distribution modelling

# setwd: sdm directory
setwd("C:/Users/a1235304/Dropbox/Short-nosed and Leaf-scaled sea snake TSSC/atm-analysis/sdm")

# install/load packages
packages <- c("sf","leaflet","readr","janitor","dplyr",
              "mapview","spatstat","tidyverse","raster",
              "dismo","lubridate","SDMtune","readxl",
              "terra","stars","lwgeom","maptools","ggspatial",
              "prettymapr","tidyterra","caret","corrplot",
              "plotROC","ggpubr")

# load package list all at once
invisible(lapply(packages, library, character.only = TRUE))

#### 1. Input model extent ####

# Load nwshelf shapefile
nw_shelf <- st_read("./nw-shelf/NWShelf.shp", quiet = TRUE) %>% 
  st_transform(4326)

# preview nwshelf
plot(nw_shelf)



#### 2. Input occurrence data ####

# Occurrence data sources:
# - `trawled_seasnakes.xlsx` (2009, 2014 - 2021)
# - Atlas of Living Australia records (downloaded 21 January 2025)
# - ATM_DPIRD Fisheries surveys (June 2024 - December 2025)
# - One spreadsheet (KLS lab catalogue records)
#
# Coordinate values: 
# - point coordinates 
# - mid latitude and longitude values (average of start/end trawl coordinates)

occurrence <- read_csv("./occurrence-data/ATM_master-occurrence-dataset.csv")

# [function] filter by species
filter_by_species <- function(species) {
  occurrence %>%
    mutate(
      lat = as.numeric(MidLat),
      lat = if_else(is.na(lat), as.numeric(Latitude), lat),
      long = as.numeric(MidLong),
      long = if_else(is.na(long), as.numeric(Longitude), long)
      ) %>%
    filter(Species == species,
           !is.na(lat),
           !is.na(long)) %>%
    clean_names()
}

# filtered df by species
short_occ <- filter_by_species("apraefrontalis")
leaf_occ <- filter_by_species("foliosquama")

# generate sf from df
short_sf <- 
  short_occ %>% 
  st_as_sf(coords = c("long", "lat"), 
           crs = 4326) %>% 
  distinct(.keep_all = T)

leaf_sf <- 
  leaf_occ %>% 
  st_as_sf(coords = c("long", "lat"), 
           crs = 4326) %>% 
  distinct(.keep_all = T)

# map distribution and colour by species
ggplot() +
  geom_sf(data = nw_shelf, fill = NA) +
  geom_sf(data = short_sf, col = "orange", cex = 3, pch = 16) +
  geom_sf(data = leaf_sf, col = "maroon", cex = 3, pch = 16) + 
  annotation_scale(mapping = aes(location = "br")) +
  theme_bw()



#### 3. Bias layer ####

# Generate bias layer from point occurrences 
# using `2020-05-20_SnakeOcc.csv`

# read bias point occurrences
bias_pts <- read_csv("./bias-layer/2020-05-20_SnakeOcc.csv")

# generate sf object from bias_pts
bias_sf <- bias_pts %>% 
  st_as_sf(coords = c("longitude", "latitude"),
           crs = 4326) %>% 
  distinct(.keep_all = T)

# visualise
ggplot() +
  geom_sf(data = bias_sf)

# vectorise bias_sf
bias_vect <- vect(bias_sf)

# mask bias points with nwshelf boundary
nw_bias_sf <- mask(bias_vect, mask = vect(nw_shelf)) %>% 
  st_as_sf()

# visualise
ggplot() +
  geom_sf(data = nw_shelf, fill = NA) +
  geom_sf(data = nw_bias_sf, aes(col = species))

# generate point pattern object
bias_ppp <- nw_bias_sf %>% 
  st_transform(crs = 3577) %>% 
  st_as_sf() %>% 
  as.ppp()

# calculate probability density of points
bias_prob <- bias_ppp %>% 
  density(., sigma = 0.05) %>% 
  rast()

plot(bias_prob)
crs(bias_prob) <- "EPSG:3577"

bias_prob <- bias_prob %>% 
  project(., y= "EPSG:4326") %>% 
  mask(mask = vect(nw_shelf)) %>% 
  resample(x = ., y = rast(bathymetry)) %>% 
  terra::scale()
  
bias_prob[values(bias_prob) < 0] <- NA

nx_bias <- minmax(bias_prob)
bias_prob <- (bias_prob - nx_bias[1,]) / (nx_bias[2,] - nx_bias[1,])

# visualise bias layer
ggplot() +
  geom_spatraster(data = bias_prob) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = nw_shelf, fill = NA) + 
  annotation_scale(mapping = aes(location = "br")) +
  theme_bw()



#### 4. Background/Pseudoabsence layer ####

# Leaf-scaled sea snake ####

# Generate 1,000 random background points within the bias layer raster 
# and excludes areas where the leaf-scaled sea snake is known to be present
# note: different output every time as seed is not set

leaf_bgpts <- 
  dismo::randomPoints(mask = raster(bias_prob), 
                      n = 1000, 
                      p = as_Spatial(leaf_sf), 
                      prob = TRUE) %>% 
  as_tibble() %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326)

# visualise
ggplot() +
  geom_spatraster(data = bias_prob) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = nw_shelf, fill = NA) + 
  geom_sf(data = leaf_bgpts, col = "aquamarine3", cex = 0.01) +
  annotation_scale(mapping = aes(location = "br")) +
  theme_bw() +
  labs(title = "Background points for Leaf-scaled sea snake within bias layer (excl. presence points)") +
  theme(plot.title = element_text(size = 11))


# Generate 1,000 random points within the nw_shelf
# note: different output every time as seed is not set
leaf_bgpts_ext <-
  spsample(x = as_Spatial(nw_shelf), n = 1000, type = "random") %>%
  st_as_sf()

# visualise
ggplot() +
  geom_spatraster(data = bias_prob) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = nw_shelf, fill = NA) + 
  geom_sf(data = leaf_bgpts_ext, col = "aquamarine3", cex = 0.8) +
  annotation_scale(mapping = aes(location = "br")) +
  theme_bw() +
  labs(title = "Background points for Leaf-scaled sea snake across model extent") +
  theme(plot.title = element_text(size = 11))

# combine two background points layer
leaf_bgpts_comb <- rbind(leaf_bgpts, leaf_bgpts_ext)

# visualise
ggplot() +
  geom_spatraster(data = bias_prob) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = nw_shelf, fill = NA) + 
  geom_sf(data = leaf_bgpts_comb, col = "aquamarine3", cex = 0.8) +
  annotation_scale(mapping = aes(location = "br")) +
  theme_bw() +
  labs(title = "Combined background points for Leaf-scaled sea snake") +
  theme(plot.title = element_text(size = 11))

# Short-nosed sea snake ####

# Generate 1,000 random background points within the bias layer raster 
# and excludes areas where the leaf-scaled sea snake is known to be present
# note: different output every time as seed is not set

short_bgpts <- 
  dismo::randomPoints(mask = raster(bias_prob), 
                      n = 1000, 
                      p = as_Spatial(short_sf), 
                      prob = TRUE) %>% 
  as_tibble() %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326)

# visualise
ggplot() +
  geom_spatraster(data = bias_prob) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = nw_shelf, fill = NA) + 
  geom_sf(data = short_bgpts, col = "aquamarine3", cex = 0.01) +
  annotation_scale(mapping = aes(location = "br")) +
  theme_bw() +
  labs(title = "Background points for Short-nosed sea snake within bias layer (excl. presence points)") +
  theme(plot.title = element_text(size = 11))


# Generate 1,000 random points within the nw_shelf
# note: different output every time as seed is not set
short_bgpts_ext <-
  spsample(x = as_Spatial(nw_shelf), n = 1000, type = "random") %>%
  st_as_sf()

# visualise
ggplot() +
  geom_spatraster(data = bias_prob) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = nw_shelf, fill = NA) + 
  geom_sf(data = short_bgpts_ext, col = "aquamarine3", cex = 0.8) +
  annotation_scale(mapping = aes(location = "br")) +
  theme_bw() +
  labs(title = "Background points for Short-nosed sea snake across model extent") +
  theme(plot.title = element_text(size = 11))

# combine two background points layer
short_bgpts_comb <- rbind(short_bgpts, short_bgpts_ext)

# visualise
ggplot() +
  geom_spatraster(data = bias_prob) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = nw_shelf, fill = NA) + 
  geom_sf(data = short_bgpts_comb, col = "aquamarine3", cex = 0.8) +
  annotation_scale(mapping = aes(location = "br")) +
  theme_bw() +
  labs(title = "Combined background points for Short-nosed sea snake") +
  theme(plot.title = element_text(size = 11))


#### 5. Input predictor rasters ####

# Load initial set of predictor rasters
env_init <- stack(
  
  # environmental variables
  "./predictor-variables/sal_mean.asc",
  "./predictor-variables/sal_amp.asc",
  "./predictor-variables/bathymetry.asc",
  "./predictor-variables/sst_mean.asc",
  "./predictor-variables/sst_amp.asc",
  "./predictor-variables/chlor_mean.asc",
  "./predictor-variables/DistToLand.asc",
  "./predictor-variables/DistToReef.asc",
  "./predictor-variables/DistToFW.asc",
  
  # genetic layer
  "../genetic_layer/laevis/output/genetic_layer.asc",
  
  # passage layer
  "../passage_layer/output/sbs_bearing_seed100_100pts_03h49m31s/passage_layer.asc"
  )



# Find correlated variables

# Test for multicollinearity
env_values <- values(env_init)
env_corr <- cor(env_values, method = "pearson", use = "complete.obs")
env_corr

# List variable names to remove
rm_vars <- findCorrelation(env_corr, cutoff = 0.7, names = T)
rm_vars # "sal_amp"    "sal_mean"   "DistToFW"   "DistToLand"

# Run again to get column number
rm_vars <- findCorrelation(env_corr, cutoff = 0.7)
rm_vars # 2 1 9 7

# List environmental variables
env_pass <- colnames(env_corr[, -rm_vars])
env_pass # "bathymetry" "sst_mean"   "sst_amp"    "chlor_mean" "DistToReef" "K2" "layer"

# Load passed environmental rasters to be used
env_use <- stack("./predictor-variables/bathymetry.asc",
                  "./predictor-variables/sst_mean.asc",
                  "./predictor-variables/sst_amp.asc",
                  "./predictor-variables/chlor_mean.asc",
                  "./predictor-variables/DistToReef.asc",
                  "../genetic_layer/laevis/output/genetic_layer.asc",
                  "../passage_layer/output/sbs_bearing_seed100_100pts_03h49m31s/passage_layer.asc")


#### 6. MaxEnt modelling ####

# prepare variables for MaxEnt modelling

# environmental predictors
env_in <- terra::rast(env_use)

# occurrence data
leaf_in <- leaf_sf %>% as_Spatial()
short_in <- short_sf %>% as_Spatial()

# background layer
leaf_bgpts_in<- leaf_bgpts_comb %>% as_Spatial()
short_bgpts_in<- short_bgpts_comb %>% as_Spatial()

# hyperparameter tuning
model_leaf <- prepareSWD(species = "species",
                         p = coordinates(leaf_in),
                         a = coordinates(bgpts_in),
                         # ensure data generated via `terra`
                         env = env_in)

# create k-folds cross validation
set_k <- 5
folds <- randomFolds(data = model_leaf,
                     k = set_k,
                     only_presence = T)

# train initial cv model
trained_model <- SDMtune::train(method = "Maxent",
                                data = model_leaf,
                                folds = folds)
# *** Warning: File absence has value -9999, treating as no-data value ***

# Assign hyperparameters to tune (replicate ENMevaluate() function parameters)
hyper <- list(reg = seq(0.5, 5, 0.5), 
              fc = c("lqh","lph","lqph"))

# Optimise model using a genetic algorithm (quicker than ENMevaluate())
opt_mod <- optimizeModel(model = trained_model,
                         hypers = hyper,
                         metric = "auc")
# *** Warning: File absence has value -9999, treating as no-data value ***

# Table of tuning results
opt_mod@results

# Select best optimised model
best_mod <- opt_mod@models[[which.max(opt_mod@results$test_AUC)]]



# 7. Model evaluation ####

# ROC curve for best model
plotROC(best_mod@models[[1]])

# Model accuracy metrics
AUC <- auc(best_mod) 
TSS <- tss(best_mod)

# Estimate thresholds
thresh <- thresholds(best_mod@models[[1]], type = "cloglog")

# Model evaulation metrics
mod <-
  opt_mod@results %>%
  slice(which.max(opt_mod@results$test_AUC)) %>% 
  as_tibble() %>%
  mutate(TSS = thresh[3,2],
         LPT = thresh[1,2])
mod

# Variable importance
vi <- maxentVarImp(best_mod)

vi %>%
  ggplot(aes(x = reorder(Variable, Percent_contribution), y = Percent_contribution)) +
  geom_bar(stat = "identity") +
  labs(y = "Variable contribution (%)", x = "", title = expression(italic("Variable importance for Leaf-scaled sea snake "))) +
       #subtitle = "Spatial and temporal scale") +
  coord_flip() +
  theme_bw()

# Response curves
env_vars <- names(env_in)

for(p in 1:length(env_vars)){
  if(p %in% 1){
    plotlist <- list()
    pb <- txtProgressBar(min = 0, max = length(env_vars), style = 3)}
  plotlist[[p]] <- plotResponse(best_mod, var = env_vars[p], rug = T)
  setTxtProgressBar(pb, p)
}

resp_curves <- ggarrange(plotlist = plotlist)
annotate_figure(resp_curves, 
                top = text_grob(paste0("Response curves for Leaf-scaled sea snake"), 
                                face = "bold", size = 14))
resp_curves


# 8. Spatial predictions ####

leaf_predict <- predict(best_mod, data = env_in, 
                        fun = c("mean", "sd"), 
                        type = "cloglog", parallel = T)

# Mean and variance in spatial prediction
leaf_predict_mean <- leaf_predict$mean
leaf_predict_sd <- leaf_predict$sd

# thresholded models
leaf_TSS_me <- leaf_predict_mean > mod$TSS
leaf_LPT_me <- leaf_predict_mean > mod$LPT
leaf_thresh_me <- raster::cut(raster(leaf_predict_mean), breaks = c(-Inf, 0.25, 0.5, 0.76, Inf))

# Your binary raster
leaf_r <- leaf_TSS_me  # or thresh if that’s the name

# Count how many TRUE (or 1) cells
leaf_n_cells <- global(leaf_r, "sum", na.rm = TRUE)
leaf_n_cells

# Get cell area in km²
leaf_cell_area <- cellSize(leaf_r, unit = "km")

# Multiply cell area by suitability mask
leaf_suitable_area <- mask(leaf_cell_area, leaf_r, maskvalues = FALSE)

# Sum total area
leaf_total_area <- global(leaf_suitable_area, "sum", na.rm = TRUE)
leaf_total_area

ggplot() +
  geom_spatraster(data = leaf_predict_mean) + 
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  annotation_scale(mapping = aes(location = "br")) +
  theme_bw() +
  labs(title = "Mean spatial prediction for Leaf-scaled sea snake")



# * MaxEnt workflow function ####

runMaxEnt <- function(name, envlyr, occdat, bgpts, k_folds){
  
  message("Running `runMaxent` function...")
  message("Hyperparameter tuning...")
  
  #hyperparameter tuning
  model <- SDMtune::prepareSWD(species = "species",
                               p = coordinates(occdat),
                               a = coordinates(bgpts),
                               env = envlyr)
  
  message("Creating k-folds cross validation...")
  # create k-folds cross validation
  folds <- randomFolds(data = model,
                       k = k_folds,
                       only_presence = T)
  
  message("Training initial model...")
  # train initial cross-validation model
  # train initial cv model
  trained_model <- SDMtune::train(method = "Maxent",
                                  data = model,
                                  folds = folds)

  # *** Warning: File absence has value -9999, treating as no-data value ***
  
  message("Assigning hyperparameters to tune...")
  # Assign hyperparameters to tune (replicate ENMevaluate() function parameters)
  hyper <- list(reg = seq(0.5, 5, 0.5),
                fc = c("lqh","lph","lqph"))
  
  message("Optimising model using a genetic algorithm...")
  # Optimise model using a genetic algorithm (quicker than ENMevaluate())
  opt_mod <- optimizeModel(model = trained_model,
                           hypers = hyper,
                           metric = "auc")
  # *** Warning: File absence has value -9999, treating as no-data value ***
  
  message("Generating table of tuning results...")
  # Table of tuning results
  tuning_res <- opt_mod@results

  message("Selecting best optimised model...")  
  # Select best optimised model
  best_mod <- opt_mod@models[[which.max(opt_mod@results$test_AUC)]]
  
  message("Generating ROC plot...")
  # model evaluation
  # ROC curve for best model
  roc_plot <- plotROC(best_mod@models[[1]])
  
  message("Calculating AUC and TSS...")
  # Model accuracy metrics
  AUC <- auc(best_mod)
  TSS <- tss(best_mod)
  
  message("Estimating thresholds...")
  # Estimate thresholds
  thresh <- thresholds(best_mod@models[[1]], type = "cloglog")
  
  message("Calculating model evaluation metrics...")
  # Model evaulation metrics
  mod <-
    opt_mod@results %>%
    slice(which.max(opt_mod@results$test_AUC)) %>% 
    as_tibble() %>%
    mutate(TSS = thresh[3,2],
           LPT = thresh[1,2])
  
  message("Determining variable importance...")
  # Variable importance
  vi <- maxentVarImp(best_mod)
  
  vi_plot <- vi %>%
    ggplot(aes(x = reorder(Variable, Percent_contribution), y = Percent_contribution)) +
    geom_bar(stat = "identity") +
    labs(y = "Variable contribution (%)", x = "", 
         title = paste0("Variable importance for ", name)) +
    #subtitle = "Spatial and temporal scale") +
    coord_flip() +
    theme_bw()

  message("Generating response curves...")
  # Response curves
  env_vars <- names(envlyr)
  
  for(p in 1:length(env_vars)){
    if(p %in% 1){
      plotlist <- list()
      pb <- txtProgressBar(min = 0, max = length(env_vars), style = 3)}
    plotlist[[p]] <- plotResponse(best_mod, var = env_vars[p], rug = T)
    setTxtProgressBar(pb, p)
  }
  
  resp_curves <- ggarrange(plotlist = plotlist)
  resp_curv_plot <- 
    annotate_figure(resp_curves, 
                    top = text_grob(paste0("Response curves for ", name), 
                                    face = "bold", size = 14))
  
  message("Performing spatial prediction...")
  predict <- predict(best_mod, data = env_in, 
                     fun = c("mean", "sd"), 
                     type = "cloglog", parallel = T)
  
  # Mean and variance in spatial prediction
  predict_mean <- predict$mean
  predict_sd <- predict$sd
  
  # thresholded models
  TSS_me <- predict_mean > mod$TSS
  LPT_me <- predict_mean > mod$LPT
  thresh_me <- raster::cut(raster(predict_mean), breaks = c(-Inf, 0.25, 0.5, 0.76, Inf))
  
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
  
  mean_pred_plot <- 
    ggplot() +
    geom_spatraster(data = predict_mean) + 
    scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
    annotation_scale(mapping = aes(location = "br")) +
    theme_bw() +
    labs(title = paste0("Mean spatial prediction for ", name))
  
  # output to print
  
  message(paste0("Generating output objects for ", name,"."))
  
  name_friendly <- gsub(" ", "_", name)
  name_friendly <- gsub("-", "_", name_friendly)
  
  assign(x = name_friendly, 
         value = list(best_mod = best_mod,
                      tuning_res = tuning_res,
                      AUC = AUC,
                      TSS = TSS,
                      thresh = thresh,
                      mod_eval_metrics = mod,
                      var_importance = vi,
                      roc_plot = roc_plot,
                      vi_plot = vi_plot,
                      response_curv = resp_curv_plot,
                      mean_predict_plot = mean_pred_plot),
         envir = .GlobalEnv)
  
    message(paste0("Output stored in `", name_friendly, "`. Access output with ", name_friendly, "$"))
    message("Complete.")
}

runMaxEnt(name = "Short-nosed sea snake",
          envlyr = env_in,
          occdat = short_in,
          bgpts = short_bgpts_in,
          k_folds = 5)

runMaxEnt(name = "Leaf-scaled sea snake",
          envlyr = env_in,
          occdat = leaf_in,
          bgpts = leaf_bgpts_in,
          k_folds = 5)



