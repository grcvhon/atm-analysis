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
# using `2020-05-20_SnakeOcc.csv` -- * updated to 2024-12-01_SnakeOcc(02).csv

# read bias point occurrences
bias_pts <- read_csv("./bias-layer/2024-12-01_SnakeOcc(02).csv")
bias_pts$longitude <- as.numeric(bias_pts$longitude)
bias_pts <- na.omit(bias_pts)

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



#### 5. Input predictor rasters ####

# Load initial set of predictor rasters
env_init <- stack(
  
  # environmental variables
  "./predictor-variables/sal_mean.asc",
  "./predictor-variables/sal_amp.asc",
  "./predictor-variables/bathymetry.asc",
  "./predictor-variables/sst_mean.asc",
  #"./predictor-variables/sst_amp.asc",
  "./predictor-variables/chlor_mean.asc",
  "./predictor-variables/DistToLand.asc",
  "./predictor-variables/DistToReef.asc",
  "./predictor-variables/DistToFW.asc",
  
  # genetic layer
  "../genetic_layer/laevis/output/genetic_layer.asc",
  
  # passage layer
  "../passage_layer/output/sbs_bearing_seed100_100pts_03h49m31s/passage_layer.asc",
  
  # conductance layer
  "../passage_layer/circuitscape/cs_output/leaf_cs_main/leaf_cs_main_cum_curmap.asc"
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
rm_vars

# List environmental variables
env_pass <- colnames(env_corr[, -rm_vars])
env_pass # "bathymetry" "sst_mean"   "sst_amp"    "chlor_mean" "DistToReef" "K2" "layer"

# Load passed environmental rasters to be used
env_use <- stack("./predictor-variables/bathymetry.asc",
                 "./predictor-variables/sst_mean.asc",
                 #"./predictor-variables/sst_amp.asc",
                 "./predictor-variables/chlor_mean.asc",
                 "./predictor-variables/DistToReef.asc",
                 "../genetic_layer/laevis/output/genetic_layer.asc",
                 "../passage_layer/output/sbs_bearing_seed100_100pts_03h49m31s/passage_layer.asc",
                 "../passage_layer/circuitscape/cs_output/leaf_cs_main/leaf_cs_main_cum_curmap.asc")



#### 6. MaxEnt modelling ####

# prepare variables for MaxEnt modelling

# environmental predictors
env_in <- terra::rast(env_use)

# occurrence data
leaf_in <- leaf_sf %>% as_Spatial()
#short_in <- short_sf %>% as_Spatial()

# background layer
leaf_bgpts_in <- leaf_bgpts_comb %>% as_Spatial()
#short_bgpts_in<- short_bgpts_comb %>% as_Spatial()




### test function ####

test_runMaxEnt(name = "Leaf-scaled sea snake",
               envlyr = env_in,
               occdat = leaf_in,
               bgpts = leaf_bgpts_in,
               k_folds = 5)





### *** runMaxEnt function *** ####

test_runMaxEnt <- function(name, envlyr, occdat, bgpts, k_folds){
  
  # set workdir
  setwd("C:/Users/a1235304/Dropbox/Short-nosed and Leaf-scaled sea snake TSSC/atm-analysis/sdm/output/")
  
  name_frnd <- gsub(" ", "_", name)
  name_frnd <- gsub("-", "_", name_frnd)
  
  today <- today()
  
  # generate output directory
  today_dirs <- dir(pattern = paste0("^", name_frnd, "_", today, "_"))
  
  n <- length(today_dirs) + 1
  final_dir <- paste0(name_frnd, "_", today, "_", n)
  
  dir.create(final_dir)
  
  sink(file = file.path(final_dir, paste0(final_dir, "_(99)_run_log.txt")),
       split = TRUE)
  
  # start
  print.noquote("Running `runMaxent` function...")
  print.noquote("Hyperparameter tuning...")
  
  #hyperparameter tuning
  model <- SDMtune::prepareSWD(species = "species",
                               p = coordinates(occdat),
                               a = coordinates(bgpts),
                               env = envlyr)
  
  print.noquote("Creating k-folds cross validation...")
  # create k-folds cross validation
  folds <- randomFolds(data = model,
                       k = k_folds,
                       only_presence = T)
  
  print.noquote("Training initial model...")
  # train initial cross-validation model
  # train initial cv model
  trained_model <- SDMtune::train(method = "Maxent",
                                  data = model,
                                  folds = folds)
  
  # *** Warning: File absence has value -9999, treating as no-data value ***
  
  print.noquote("Assigning hyperparameters to tune...")
  # Assign hyperparameters to tune (replicate ENMevaluate() function parameters)
  hyper <- list(reg = seq(0.5, 5, 0.5),
                fc = c("lqh","lph","lqph"))
  
  print.noquote("Optimising model using a genetic algorithm...")
  # Optimise model using a genetic algorithm (quicker than ENMevaluate())
  opt_mod <- optimizeModel(model = trained_model,
                           hypers = hyper,
                           metric = "auc")
  # *** Warning: File absence has value -9999, treating as no-data value ***
  
  print.noquote("Generating table of tuning results...")
  # Table of tuning results
  tuning_res <- opt_mod@results
  
  print.noquote("Selecting best optimised model...")  
  # Select best optimised model
  best_mod <- opt_mod@models[[which.max(opt_mod@results$test_AUC)]]
  
  print.noquote("Generating ROC plot...")
  # model evaluation
  # ROC curve for best model
  roc_plot <- plotROC(best_mod@models[[1]])
  
  print.noquote("Calculating AUC and TSS...")
  # Model accuracy metrics
  AUC <- auc(best_mod)
  TSS <- tss(best_mod)
  
  print.noquote("Estimating thresholds...")
  # Estimate thresholds
  thresh <- thresholds(best_mod@models[[1]], type = "cloglog")
  
  print.noquote("Calculating model evaluation metrics...")
  # Model evaulation metrics
  mod <-
    opt_mod@results %>%
    slice(which.max(opt_mod@results$test_AUC)) %>% 
    as_tibble() %>%
    mutate(TSS = thresh[3,2],
           LPT = thresh[1,2])
  
  print.noquote("Determining variable importance...")
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
  
  print.noquote("Generating response curves...")
  # Response curves
  env_vars <- names(envlyr)
  
  for(p in 1:length(env_vars)){
    if(p %in% 1){
      plotlist <- list()
      pb <- txtProgressBar(min = 0, max = length(env_vars), style = 3)}
    plotlist[[p]] <- plotResponse(best_mod, var = env_vars[p], rug = T)
    setTxtProgressBar(pb, p)
  }
  
  close(pb)
  
  resp_curves <- ggarrange(plotlist = plotlist)
  resp_curv_plot <- 
    annotate_figure(resp_curves, 
                    top = text_grob(paste0("Response curves for ", name), 
                                    face = "bold", size = 14))
  
  print.noquote("Performing spatial prediction...")
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
  
  print.noquote(paste0("Generating output objects for ", name,"."))
  
  name_friendly <- gsub(" ", "_", final_dir)
  name_friendly <- gsub("-", "_", name_friendly)
  
  # 1
  print.noquote(paste0("Writing ", name_friendly, "_(01)_best_mod.txt"))
  bstmdl <- capture.output(print(best_mod), type = "message")
  writeLines(bstmdl, con = file.path(final_dir, paste0(name_friendly, "_(01)_best_mod.txt")))
  
  # 2
  print.noquote(paste0("Writing ", name_friendly, "_(02)_tuning_results.csv"))
  write.csv(tuning_res, file = file.path(final_dir, paste0(name_friendly, "_(02)_tuning_results.csv")), row.names = FALSE)
  
  # 3 & 4
  print.noquote(paste0("Writing ", name_friendly, "_(03-4)_auc_tss.csv"))
  write.csv(data.frame(AUC = AUC, TSS = TSS), file = file.path(final_dir, paste0(name_friendly, "_(03-4)_auc_tss.csv")), row.names = FALSE)
  
  # 5
  print.noquote(paste0("Writing ", name_friendly, "_(05)_thresholds.csv"))
  write.csv(thresh, file = file.path(final_dir, paste0(name_friendly, "_(05)_thresholds.csv")), row.names = FALSE)
  
  # 6
  print.noquote(paste0("Writing ", name_friendly, "_(06)_thresh_me.png"))
  png(filename = file.path(final_dir, paste0(name_friendly, "_(06)_thresh_me.png")), width = 1200, height = 800)
  plot(thresh_me)
  dev.off()

  # 7
  print.noquote(paste0("Writing ", name_friendly, "_(07)_mod_eval_metrics.csv"))
  write.csv(mod, file = file.path(final_dir, paste0(name_friendly, "_(07)_mod_eval_metrics.csv")), row.names = FALSE)
  
  # 8
  print.noquote(paste0("Writing ", name_friendly, "_(08)_var_importance.csv"))
  write.csv(vi, file = file.path(final_dir, paste0(name_friendly, "_(08)_var_importance.csv")), row.names = FALSE)
  
  # 9
  print.noquote(paste0("Writing ", name_friendly, "_(09)_roc_plot.png"))
  ggsave(file.path(final_dir, paste0(name_friendly, "_(09)_roc_plot.png")), plot = roc_plot, width = 8, height = 6, dpi = 300)
  
  # 10
  print.noquote(paste0("Writing ", name_friendly, "_(10)_vi_plot.png"))
  ggsave(file.path(final_dir, paste0(name_friendly, "_(10)_vi_plot.png")), plot = vi_plot, width = 8, height = 6, dpi = 300)
  
  # 11
  print.noquote(paste0("Writing ", name_friendly, "_(11)_response_curves.png"))
  ggsave(file.path(final_dir, paste0(name_friendly, "_(11)_response_curves.png")), plot = resp_curv_plot, width = 12, height = 8, dpi = 300)
  
  # 12
  print.noquote(paste0("Writing ", name_friendly, "_(12)_mean_prediction.png"))
  ggsave(file.path(final_dir, paste0(name_friendly, "_(12)_mean_prediction.png")), plot = mean_pred_plot, width = 8, height = 6, dpi = 300)
  
  # 13
  
  # 13.14
  print.noquote(paste0("Writing ", name_friendly, "_(14)_mean.png"))
  png(filename = file.path(final_dir, paste0(name_friendly, "_(14)_mean.png")), width = 1200, height = 800)
  plot(predict_mean, main = paste0("Mean prediction: ", name))
  dev.off()
  
  # 13.15
  print.noquote(paste0("Writing ", name_friendly, "_(15)_sd.png"))
  png(filename = file.path(final_dir, paste0(name_friendly, "_(15)_sd.png")), width = 1200, height = 800)
  plot(predict_sd, main = paste0("Prediction SD: ", name))
  dev.off()
  
  # 16
  print.noquote(paste0("Writing ", name_friendly, "_(16)_TSS_binary_(TSS_me).png"))
  png(filename = file.path(final_dir, paste0(name_friendly, "_(16)_TSS_binary_(TSS_me).png")), width = 1200, height = 800)
  plot(TSS_me, main = paste0("TSS thresholded: ", name))
  dev.off()
  
  # 17
  print.noquote(paste0("Writing ", name_friendly, "_(17)_LPT_binary_(LPT_me).png"))
  png(filename = file.path(final_dir, paste0(name_friendly, "_(17)_LPT_binary_(LPT_me).png")), width = 1200, height = 800)
  plot(LPT_me, main = paste0("LPT thresholded: ", name))
  dev.off()
  
  # 18
  print.noquote(paste0("Writing ", name_friendly, "_(18)_n_cells.txt"))
  sink(file = file.path(final_dir, paste0(name_friendly, "_(18)_n_cells.txt")))
  print(n_cells)
  sink()
  
  # 19
  print.noquote(paste0("Writing ", name_friendly, "_(19)_cell_area.png"))
  #ggsave(file.path(final_dir, paste0(name_friendly, "_(19)_cell_area.png")), plot = cell_area, width = 8, height = 6, dpi = 300)
  
  png(filename = file.path(final_dir, paste0(name_friendly, "_(19)_cell_area.png")), width = 1200, height = 800)
  plot(cell_area, main = paste0("Cell area: ", name))
  dev.off()
  
  sink(file = file.path(final_dir, paste0(name_friendly, "_(19)_cell_area.txt")))
  print(cell_area)
  sink()
  
  # 20
  print.noquote(paste0("Writing ", name_friendly, "_(20)_suitable_area.png"))
  #ggsave(file.path(final_dir, paste0(name_friendly, "_(20)_suitable_area.png")), plot = suitable_area, width = 8, height = 6, dpi = 300)
  
  png(filename = file.path(final_dir, paste0(name_friendly, "_(20)_suitable_area.png")), width = 1200, height = 800)
  plot(suitable_area, main = paste0("Cell area: ", name))
  dev.off()
  
  sink(file = file.path(final_dir, paste0(name_friendly, "_(20)_suitable_area.txt")))
  print(suitable_area)
  sink()
  
  # 21
  print.noquote(paste0("Writing ", name_friendly, "_(21)_total_area.txt"))
  sink(file = file.path(final_dir, paste0(name_friendly, "_(21)_total_area.txt")))
  print(total_area)
  sink()
  
  # 22
  print.noquote(paste0("Writing ", name_friendly, "_(22)_opt_mod_results.csv"))
  write.csv(opt_mod@results, file = file.path(final_dir, paste0(name_friendly, "_(22)_opt_mod_results.csv")), row.names = FALSE)
  
  # 23
  print.noquote(paste0("Writing ", name_friendly, "_(23)_hyper.txt"))
  sink(file = file.path(final_dir, paste0(name_friendly, "_(23)_hyper.txt")))
  print(hyper)
  sink()
  
  # 24
  print.noquote(paste0("Writing ", name_friendly, "_(24)_model.txt"))
  mdl <- capture.output(print(model), type = "message")
  writeLines(mdl, con = file.path(final_dir, paste0(name_friendly, "_(24)_model.txt")))
  
  # 25
  print.noquote(paste0("Writing ", name_friendly, "_(25)_trained_model.txt"))
  tmdl <- capture.output(print(trained_model), type = "message")
  writeLines(tmdl, con = file.path(final_dir, paste0(name_friendly, "_(25)_trained_model.txt")))
  
  # store into an object
  final_friendly <- gsub(" ", "_", final_dir)
  final_friendly <- gsub("-", "_", final_friendly)
  
  print.noquote(paste0("Storing output as an R object (", final_friendly, ")."))
  
  assign(x = final_friendly, 
         value = list(best_mod = best_mod, # 1
                      tuning_res = tuning_res, # 2
                      AUC = AUC, # 3
                      TSS = TSS, # 4
                      thresh = thresh, # 5
                      thresh_me = thresh_me, # 6
                      mod_eval_metrics = mod, # 7
                      var_importance = vi, # 8
                      roc_plot = roc_plot, # 9
                      vi_plot = vi_plot, # 10
                      response_curv = resp_curv_plot, # 11
                      mean_predict_plot = mean_pred_plot, # 12
                      predict = predict, # 13
                      predict_mean = predict_mean, # 14
                      predict_sd = predict_sd, # 15
                      TSS_me = TSS_me, # 16
                      LPT_me = LPT_me, # 17
                      n_cells = n_cells, # 18
                      cell_area = cell_area, # 19 - empty
                      suitable_area = suitable_area, # 20 - empty
                      total_area = total_area, # 21
                      opt_mod = opt_mod, # 22
                      hyper = hyper, # 23
                      model = model, # 24 -
                      folds = folds, # 25 -
                      trained_model = trained_model), # 26 -
         envir = .GlobalEnv)
  
  print.noquote(paste0("Output stored in `", final_friendly, "`. Access output with ", final_friendly, "$"))
  print.noquote("Complete.")

  sink()
}




