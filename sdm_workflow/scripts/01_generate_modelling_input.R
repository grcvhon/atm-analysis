# sdm - 01 generate modelling input

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




##############################################################################################################
#### 1. Input model extent ###################################################################################
##############################################################################################################

# Load nwshelf shapefile
nw_shelf <- st_read("./nw-shelf/NWShelf.shp", quiet = TRUE) %>% 
  st_transform(4326)





##############################################################################################################
#### 2. Input occurrence data ################################################################################
##############################################################################################################

# Occurrence data sources:
# - `trawled_seasnakes.xlsx` (2009, 2014 - 2021)
# - Atlas of Living Australia records (downloaded 21 January 2025)
# - ATM_DPIRD Fisheries surveys (June 2024 - December 2025)
# - One spreadsheet (KLS lab catalogue records)
#
# Coordinate values: 
# - point coordinates 
# - mid latitude and longitude values (average of start/end trawl coordinates)

#occurrence <- read_csv("./occurrence-data/ATM_master-occurrence-dataset.csv")
#
## [function] filter by species
#filter_by_species <- function(species) {
#  occurrence %>%
#    mutate(
#      lat = as.numeric(MidLat),
#      lat = if_else(is.na(lat), as.numeric(Latitude), lat),
#      long = as.numeric(MidLong),
#      long = if_else(is.na(long), as.numeric(Longitude), long)
#    ) %>%
#    filter(Species == species,
#           !is.na(lat),
#           !is.na(long)) %>%
#    clean_names()
#}
#
## filtered df by species
#short_occ <- filter_by_species("apraefrontalis")
#leaf_occ <- filter_by_species("foliosquama")
#
## generate sf from df
#short_sf <- 
#  short_occ %>% 
#  st_as_sf(coords = c("long", "lat"), 
#           crs = 4326) %>% 
#  distinct(.keep_all = T)
#
#leaf_sf <- 
#  leaf_occ %>% 
#  st_as_sf(coords = c("long", "lat"), 
#           crs = 4326) %>% 
#  distinct(.keep_all = T)
#
## map distribution and colour by species
#ggplot() +
#  geom_sf(data = nw_shelf, fill = NA) +
#  geom_sf(data = short_sf, col = "orange", cex = 3, pch = 16) +
#  geom_sf(data = leaf_sf, col = "maroon", cex = 3, pch = 16) + 
#  annotation_scale(mapping = aes(location = "br")) +
#  theme_bw()
#
## write occurrence points to file
#sf::st_write(obj = leaf_sf, dsn = "./model_input/occurrence/leaf_occurrence.shp")
#sf::st_write(obj = short_sf, dsn = "./model_input/occurrence/short_occurrence.shp")

#  ! written to file - commented out to avoid accidentally running.



##############################################################################################################
#### 3. Bias layer ###########################################################################################
##############################################################################################################

# Generate bias layer from point occurrences 
# using `2020-05-20_SnakeOcc.csv` -- * updated to 2024-12-01_SnakeOcc(02).csv

## read bias point occurrences
#bias_pts <- read_csv("./bias-layer/2024-12-01_SnakeOcc(02).csv")
#bias_pts$longitude <- as.numeric(bias_pts$longitude)
#bias_pts <- na.omit(bias_pts)
#
## generate sf object from bias_pts
#bias_sf <- bias_pts %>% 
#  st_as_sf(coords = c("longitude", "latitude"),
#           crs = 4326) %>% 
#  distinct(.keep_all = T)
#
## visualise
#ggplot() +
#  geom_sf(data = bias_sf)
#
## vectorise bias_sf
#bias_vect <- vect(bias_sf)
#
## mask bias points with nwshelf boundary
#nw_bias_sf <- mask(bias_vect, mask = vect(nw_shelf)) %>% 
#  st_as_sf()
#
## visualise
#ggplot() +
#  geom_sf(data = nw_shelf, fill = NA) +
#  geom_sf(data = nw_bias_sf, aes(col = species))
#
## generate point pattern object
#bias_ppp <- nw_bias_sf %>% 
#  st_transform(crs = 3577) %>% 
#  st_as_sf() %>% 
#  as.ppp()
#
## calculate probability density of points
#bias_prob <- bias_ppp %>% 
#  density(., sigma = 0.05) %>% 
#  rast()
#
#plot(bias_prob)
#crs(bias_prob) <- "EPSG:3577"
#
#bias_prob <- bias_prob %>% 
#  project(., y= "EPSG:4326") %>% 
#  mask(mask = vect(nw_shelf)) %>% 
#  resample(x = ., y = rast("./predictor-variables/bathymetry.asc")) %>% 
#  terra::scale()
#
#bias_prob[values(bias_prob) < 0] <- NA
#
#nx_bias <- minmax(bias_prob)
#bias_prob <- (bias_prob - nx_bias[1,]) / (nx_bias[2,] - nx_bias[1,])
#
## visualise bias layer
#ggplot() +
#  geom_spatraster(data = bias_prob) +
#  scale_fill_viridis_c(na.value = "transparent") +
#  geom_sf(data = nw_shelf, fill = NA) + 
#  annotation_scale(mapping = aes(location = "br")) +
#  theme_bw()
#
## write bias layer to file
#terra::writeRaster(bias_prob, 
#  filename = "./model_input/bias/bias_layer.asc")

#  ! written to file - commented out to avoid accidentally running.



##############################################################################################################
#### 4. Background/Pseudoabsence layer #######################################################################
##############################################################################################################

##### Background: Leaf-scaled sea snake #####

# Generate 1,000 random background points within the bias layer raster 
# and excludes areas where the leaf-scaled sea snake is known to be present
# run once, write to file

#leaf_bgpts <- 
#  dismo::randomPoints(mask = raster(bias_prob), 
#                      n = 1000, 
#                      p = as_Spatial(leaf_sf), 
#                      prob = TRUE) %>% 
#  as_tibble() %>% 
#  st_as_sf(coords = c("x", "y"), crs = 4326)
#
## visualise
#ggplot() +
#  geom_spatraster(data = bias_prob) +
#  scale_fill_viridis_c(na.value = "transparent") +
#  geom_sf(data = nw_shelf, fill = NA) + 
#  geom_sf(data = leaf_bgpts, col = "aquamarine3", cex = 0.01) +
#  annotation_scale(mapping = aes(location = "br")) +
#  theme_bw() +
#  labs(title = "Background points for Leaf-scaled sea snake within bias layer (excl. presence points)") +
#  theme(plot.title = element_text(size = 11))
#
#
## Generate 1,000 random points within the nw_shelf
## note: different output every time as seed is not set
#leaf_bgpts_ext <-
#  spsample(x = as_Spatial(nw_shelf), n = 1000, type = "random") %>%
#  st_as_sf()
#
## visualise
#ggplot() +
#  geom_spatraster(data = bias_prob) +
#  scale_fill_viridis_c(na.value = "transparent") +
#  geom_sf(data = nw_shelf, fill = NA) + 
#  geom_sf(data = leaf_bgpts_ext, col = "aquamarine3", cex = 0.8) +
#  annotation_scale(mapping = aes(location = "br")) +
#  theme_bw() +
#  labs(title = "Background points for Leaf-scaled sea snake across model extent") +
#  theme(plot.title = element_text(size = 11))
#
## combine two background points layer
#leaf_bgpts_comb <- rbind(leaf_bgpts, leaf_bgpts_ext)
#
## visualise
#ggplot() +
#  geom_spatraster(data = bias_prob) +
#  scale_fill_viridis_c(na.value = "transparent") +
#  geom_sf(data = nw_shelf, fill = NA) + 
#  geom_sf(data = leaf_bgpts_comb, col = "aquamarine3", cex = 0.8) +
#  annotation_scale(mapping = aes(location = "br")) +
#  theme_bw() +
#  labs(title = "Combined background points for Leaf-scaled sea snake") +
#  theme(plot.title = element_text(size = 11))
#
#
#
###### Background: Short-nosed sea snake #####
#
## Generate 1,000 random background points within the bias layer raster 
## and excludes areas where the leaf-scaled sea snake is known to be present
## run once, write to file
#
#short_bgpts <- 
#  dismo::randomPoints(mask = raster(bias_prob), 
#                      n = 1000, 
#                      p = as_Spatial(short_sf), 
#                      prob = TRUE) %>% 
#  as_tibble() %>% 
#  st_as_sf(coords = c("x", "y"), crs = 4326)
#
## visualise
#ggplot() +
#  geom_spatraster(data = bias_prob) +
#  scale_fill_viridis_c(na.value = "transparent") +
#  geom_sf(data = nw_shelf, fill = NA) + 
#  geom_sf(data = short_bgpts, col = "aquamarine3", cex = 0.01) +
#  annotation_scale(mapping = aes(location = "br")) +
#  theme_bw() +
#  labs(title = "Background points for Short-nosed sea snake within bias layer (excl. presence points)") +
#  theme(plot.title = element_text(size = 11))
#
#
## Generate 1,000 random points within the nw_shelf
## note: different output every time as seed is not set
#short_bgpts_ext <-
#  spsample(x = as_Spatial(nw_shelf), n = 1000, type = "random") %>%
#  st_as_sf()
#
## visualise
#ggplot() +
#  geom_spatraster(data = bias_prob) +
#  scale_fill_viridis_c(na.value = "transparent") +
#  geom_sf(data = nw_shelf, fill = NA) + 
#  geom_sf(data = short_bgpts_ext, col = "aquamarine3", cex = 0.8) +
#  annotation_scale(mapping = aes(location = "br")) +
#  theme_bw() +
#  labs(title = "Background points for Short-nosed sea snake across model extent") +
#  theme(plot.title = element_text(size = 11))
#
## combine two background points layer
#short_bgpts_comb <- rbind(short_bgpts, short_bgpts_ext)
#
## visualise
#ggplot() +
#  geom_spatraster(data = bias_prob) +
#  scale_fill_viridis_c(na.value = "transparent") +
#  geom_sf(data = nw_shelf, fill = NA) + 
#  geom_sf(data = short_bgpts_comb, col = "aquamarine3", cex = 0.8) +
#  annotation_scale(mapping = aes(location = "br")) +
#  theme_bw() +
#  labs(title = "Combined background points for Short-nosed sea snake") +
#  theme(plot.title = element_text(size = 11))
#
#
## write leaf-scaled background points to file
#sf::st_write(obj = leaf_bgpts_comb, dsn = "./model_input/background/leaf_background.shp")
#
## write short-nosed background points to file
#sf::st_write(obj = short_bgpts_comb, dsn = "./model_input/background/short_background.shp")

#  ! written to file - commented out to avoid accidentally running.



##############################################################################################################
#### 5. Input predictor rasters ##############################################################################
##############################################################################################################

#sal_mean <- terra::rast("./predictor-variables/sal_mean.asc")
#names(sal_mean) <- "Mean salinity"
#
#bathy_layer <- terra::rast("./predictor-variables/bathymetry.asc")
#names(bathy_layer) <- "Bathymetry"
#
#sst_mean <- terra::rast("./predictor-variables/sst_mean.asc")
#names(sst_mean) <- "Mean SST"
#
#disttoland <- terra::rast("./predictor-variables/DistToLand.asc")
#names(disttoland) <- "Distance to land"
#
#disttoreef <- terra::rast("./predictor-variables/DistToReef.asc")
#names(disttoreef) <- "Distance to reef"
#
#disttofw <- terra::rast("./predictor-variables/DistToFW.asc")
#names(disttofw) <- "Distance to freshwater"
#
#scld_nonpass <- terra::rast("../passage_layer/output/scaled_nonpass.asc")
#names(scld_nonpass) <- "Scaled non-passage probability"
#
## scale is -1 to 1
#clust <- terra::rast("../genetic_layer/laevis/output/k1_layer.asc")
#names(clust) <- "Ancestry coefficient"
## Get min and max values
#min_clust <- min(values(clust), na.rm = TRUE)
#max_clust <- max(values(clust), na.rm = TRUE)
## Normalize to range -1 to 1
#values(clust) <- 2 * (values(clust) - min_clust) / (max_clust - min_clust) - 1
#plot(clust)
#
## generate rasterstack
#env_init <- stack(
#  raster(sal_mean),
#  raster(bathy_layer),
#  raster(sst_mean),
#  raster(disttoland),
#  raster(disttoreef),
#  raster(disttofw),
#  raster(scld_nonpass),
#  raster(clust)
#)
#
## Find correlated variables
#
## Test for multicollinearity
#env_values <- values(env_init)
#env_corr <- cor(env_values, method = "pearson", use = "complete.obs")
#env_corr
#
## List variable names to remove
#rm_vars <- findCorrelation(env_corr, cutoff = 0.7, names = T)
#rm_vars
#
## Run again to get column number
#rm_vars <- findCorrelation(env_corr, cutoff = 0.7)
#rm_vars
#
## List environmental variables
#env_pass <- colnames(env_corr[, -rm_vars])
#env_pass
#
## subset predictors which passed multicollinearity test from the initial stack of predictors
#env_use <- subset(env_init, env_pass)
## Bathymetry, Mean.SST, Distance.to.reef, Ancestry.coefficient.K1, Scaled.non.passage.probability
#
## write bias layer to file
#terra::writeRaster(env_use, 
#  filename = "./model_input/predictors/predictor_rasterstack.tif")

#  ! written to file - commented out to avoid accidentally running.