# sdm - 02 optimise model to arrive at best model

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

### read in model input ###

# predictors
env_in <- terra::rast("./model_input/predictors/predictor_rasterstack.tif")
# re-assign names to stack layers
names(env_in) <- c(
  "Bathymetry",
  "Mean.SST",
  "Distance.to.reef",
  "Scaled.non.passage.probability",
  "Ancestry.coefficient"
)

# occurrence data
leaf_in <- sf::st_read("./model_input/occurrence/leaf_occurrence.shp") %>% as_Spatial()
short_in <- sf::st_read("./model_input/occurrence/short_occurrence.shp") %>% as_Spatial()

# background layer
leaf_background <- sf::st_read("./model_input/background/leaf_background.shp") %>% as_Spatial()
short_background <- sf::st_read("./model_input/background/short_background.shp") %>% as_Spatial()

##############################################################################################################

## modelling step
#
## put in a loop
#datasets <- list(
#  short = list(p = short_in, a = short_background),
#  leaf = list(p = leaf_in, a = leaf_background)
#)
#
#models <- list()
#trained_models <- list()
#
#for (name in names(datasets)) {
#  p_coords <- coordinates(datasets[[name]]$p)
#  a_coords <- coordinates(datasets[[name]]$a)
#  
#  # Prepare SWD object
#  s <- SDMtune::prepareSWD(
#    species = "species",
#    p = p_coords,
#    a = a_coords,
#    env = env_in
#  )
#  
#  # Generate folds
#  folds <- randomFolds(data = s, k = 5, only_presence = TRUE)
#  
#  # Train the model
#  trained_model <- SDMtune::train(method = "Maxent", data = s, folds = folds)
#  
#  # Store the results
#  models[[name]] <- s
#  trained_models[[name]] <- trained_model
#}
#
## output list
#models$
#trained_models$
#
## Assign hyperparameters to tune (replicate ENMevaluate() function parameters)
#hyper <- list(reg = seq(0.5, 5, 0.5),
#              fc = c("lqh","lph","lqph"))
#
## Optimise trained_models$short using a genetic algorithm (quicker than ENMevaluate())
#short_opt_mod <- 
#  optimizeModel(
#    model = trained_models$short,
#    hypers = hyper,
#    metric = "auc"
#)
#
## Write optimised model (short_opt_mod)
#base::saveRDS(short_opt_mod, "./model_out/short_nosed_opt_mod.rds")
#
## Optimise trained_models$leaf using a genetic algorithm (quicker than ENMevaluate())
#leaf_opt_mod <- 
#  optimizeModel(
#    model = trained_models$leaf,
#    hypers = hyper,
#    metric = "auc"
#)
#
## Write optimised model (leaf_opt_mod)
#base::saveRDS(leaf_opt_mod, "./model_out/leaf_scaled_opt_mod.rds")

#  ! RDS file already written - commented out to avoid accidentally running.

##############################################################################################################

# extract species-specific opt_mod output

model_files <- c(
  "short_nosed_opt_mod.rds",
  "leaf_scaled_opt_mod.rds"
)

results_list <- list()

for (file_name in model_files) {
  # Read the model object
  model_path <- file.path("./model_out", file_name)
  model_obj <- base::readRDS(model_path)
  
  # Extract the results slot
  model_files_without_ext <- sub("\\.rds$", "", file_name)

  tuning_res <- model_obj@results
  #write.csv(tuning_res, paste0("./model_out/", model_files_without_ext,"_tuning_results.csv"))

  best_mod <- model_obj@models[[which.max(model_obj@results$test_AUC)]]
  bstmdl <- capture.output(print(best_mod), type = "message")
  #writeLines(bstmdl, con = file.path(paste0("./model_out/", model_files_without_ext,"_best_model.txt")))

  roc_plot <- plotROC(best_mod@models[[1]])
  #ggsave(file.path(paste0("./model_out/", model_files_without_ext,"_roc_plot.png")), plot = roc_plot, width = 8, height = 6, dpi = 300)

  AUC <- auc(best_mod)
  #write.csv(AUC, paste0("./model_out/", model_files_without_ext,"_AUC.csv"))

  TSS <- tss(best_mod)
  #write.csv(AUC, paste0("./model_out/", model_files_without_ext,"_TSS.csv"))

  thresh <- thresholds(best_mod@models[[1]], type = "cloglog")
  #write.csv(thresh, paste0("./model_out/", model_files_without_ext,"_thresholds.csv"))

  mod <- 
    model_obj@results %>%
    slice(which.max(model_obj@results$test_AUC)) %>% 
    as_tibble() %>%
    mutate(
      TSS = thresh[3,2],
      LPT = thresh[1,2]
    )
  #write.csv(mod, paste0("./model_out/", model_files_without_ext,"_model_evaluation_metrics.csv"))
  
  vi <- maxentVarImp(best_mod)
  #write.csv(vi, paste0("./model_out/", model_files_without_ext,"_variable_importance.csv"))

  vi_plot <- vi %>%
    ggplot(aes(x = reorder(Variable, Percent_contribution), y = Percent_contribution)) +
    geom_bar(stat = "identity") +
    labs(y = "Variable contribution (%)", x = "", title = paste0("Variable importance (", model_files_without_ext, ")")) +
    coord_flip() +
    theme_bw()
  #ggsave(file.path(paste0("./model_out/", model_files_without_ext,"_variable_importance_plot.png")), plot = vi_plot, width = 8, height = 6, dpi = 300)

  env_vars <- names(env_in)
  
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
                    top = text_grob(paste0("Response curves for ", model_files_without_ext), 
                                    face = "bold", size = 14))
  #ggsave(file.path(paste0("./model_out/", model_files_without_ext,"_response_curve_plot.png")), plot = resp_curv_plot, width = 8, height = 6, dpi = 300)
 
  
  # Store in list with name
  results_list[[file_name]] <- 
    list(
      tuning_results = tuning_res, 
      best_model = best_mod,
      roc_plot = roc_plot,
      AUC = AUC,
      TSS = TSS,
      thresholds = thresh,
      model_evaluation_metrics = mod,
      variable_importance = vi,
      variable_importance_plot = vi_plot,
      response_curve_plot = resp_curv_plot
  )
}

results_list$short_nosed_opt_mod.rds$variable_importance_plot
results_list$leaf_scaled_opt_mod.rds$