### MaxEnt modelling framework ###
### 2025-01-16 ###

### Install packages ---------------------------------------
# install.packages("sf")
# install.packages("leaflet")
# install.packages("readr")
# install.packages("janitor")
# install.packages("dplyr")
# install.packages("mapview")
# install.packages("spatstat")
# install.packages("tidyverse")
# install.packages("raster")
# install.packages("dismo")
# install.packages("lubridate")
# install.packages("SDMtune")
# install.packages("readxl")
# install.packages("stars")
# install.packages("lwgeom")
# install.packages("maptools", repos = "http://R-Forge.R-project.org")
# install.packages("ggspatial")
# install.packages("prettymapr")

### Load package list ------------------------------------------
packages <- c("sf","leaflet","readr","janitor","dplyr",
              "mapview","spatstat","tidyverse","raster",
              "dismo","lubridate","SDMtune","readxl",
              "terra","stars","lwgeom","maptools","ggspatial",
              "prettymapr")
invisible(lapply(packages, library, character.only = TRUE))

### Input model extent -----------------------------------------

# Load Northwest shelf bounding area
nw_shelf <- st_read("data/shapefiles/nw-shelf/NWShelf.shp")
mapview(nw_shelf)

### Input environmental predictors -----------------------------

# Load initial environmental predictors
env_init <- stack("data/predictor-variables/sal_mean.asc",
                  "data/predictor-variables/sal_amp.asc",
                  "data/predictor-variables/bathymetry.asc",
                  "data/predictor-variables/sst_mean.asc",
                  "data/predictor-variables/sst_amp.asc",
                  "data/predictor-variables/chlor_mean.asc",
                  "data/predictor-variables/DistToLand.asc",
                  "data/predictor-variables/DistToReef.asc",
                  "data/predictor-variables/DistToFW.asc")
env_init

### Input occurrence data --------------------------------------

# Load occurrence data
# csv contains combined occurrence data for 
# Aipysurus apraefrontalis and A. foliosquama
rawdat <- read_csv("data/spreadsheets/ATM_2023_0715-running-master.csv")

# filter for apraefrontalis
aprae <- rawdat %>%
  mutate(
    lat = as.numeric(MidLat),
    lat = if_else(is.na(lat), as.numeric(Latitude), lat),
    long = as.numeric(MidLong),
    long = if_else(is.na(long), as.numeric(Longitude), long),) %>%
  filter(Species == "apraefrontalis",
         !is.na(lat),
         !is.na(long)) %>% 
  clean_names()

# filter for foliosquama
folio <- rawdat %>%
  mutate(
    lat = as.numeric(MidLat),
    lat = if_else(is.na(lat), as.numeric(Latitude), lat),
    long = as.numeric(MidLong),
    long = if_else(is.na(long), as.numeric(Longitude), long),) %>%
  filter(Species == "foliosquama",
         !is.na(lat),
         !is.na(long)) %>% 
  clean_names()

# species-specific occurrence map + bounding area
mapview(aprae, xcol = "long", ycol = "lat", crs = 4326) + nw_shelf  
mapview(folio, xcol = "long", ycol = "lat", crs = 4326) + nw_shelf

# convert to simple feature (sf) 
aprae_sf <- aprae %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  distinct(.keep_all = T)

folio_sf <- folio %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  distinct(.keep_all = T)

### Create bias layer (for point observations) -----------------
#### apraefrontalis ####
aprae_ppp <- aprae_sf %>%
  st_transform(crs = 3577) %>% 
  as.ppp()

ggplot() + 
  annotation_map_tile("cartolight", zoom = 4) +
  layer_spatial(aprae_sf, colour = "brown") +
  layer_spatial(nw_shelf, fill = NA, col = "black", lwd = 0.3) +
  annotation_scale(width_hint = 0.2, location = "br") +
  theme_void()

# Calculate Gaussian density distribution of points
aprae_pts_bias <- aprae_ppp %>% 
  density(., sigma = 0.05) %>% 
  raster()

crs(aprae_pts_bias) <- CRS("+init=epsg:3577")

aprae_pts_bias <- aprae_pts_bias %>% 
  projectRaster(., crs = CRS("+init=epsg:4326")) %>% 
  mask(mask = nw_shelf) %>% 
  resample(x = ., y = env_init[[1]]) # env_init are initial subset of environmental predictors - Line 40

values(aprae_pts_bias) <- values(aprae_pts_bias) + min(values(aprae_pts_bias), na.rm = T)

aprae_bias_layer <- aprae_pts_bias
aprae_bias_layer[values(aprae_bias_layer) < 0] <- NA

mapview(aprae_bias_layer)

#### foliosquama ####
folio_ppp <- folio_sf %>%
  st_transform(crs = 3577) %>% 
  as.ppp()

ggplot() + 
  annotation_map_tile("cartolight", zoom = 4) +
  layer_spatial(folio_sf, colour = "brown") +
  layer_spatial(nw_shelf, fill = NA, col = "black", lwd = 0.3) +
  annotation_scale(width_hint = 0.2, location = "br") +
  theme_void()

# Calculate Gaussian density distribution of points
folio_pts_bias <- folio_ppp %>% 
  density(., sigma = 0.05) %>% 
  raster()

crs(folio_pts_bias) <- CRS("+init=epsg:3577")

folio_pts_bias <- folio_pts_bias %>% 
  projectRaster(., crs = CRS("+init=epsg:4326")) %>% 
  mask(mask = nw_shelf) %>% 
  resample(x = ., y = env_init[[1]]) # env_init are initial subset of environmental predictors - Line 40

values(folio_pts_bias) <- values(folio_pts_bias) + min(values(folio_pts_bias), na.rm = T)

folio_bias_layer <- folio_pts_bias
folio_bias_layer[values(folio_bias_layer) < 0] <- NA

mapview(folio_bias_layer)

mapview(folio, xcol = "long", ycol = "lat", crs = 4326) + folio_bias_layer

### Create bias layer (for transect lines) ---------------------

sample <- folio %>% 
  select(c(start_long,start_lat,end_long,end_lat)) %>% 
  filter(apply(folio[, c("start_long", "start_lat", "end_long", "end_lat")], 1, function(x) all(!is.na(x))))

# create an empty list to store transect line geometries
trn_lines <- list()

for(i in 1:nrow(sample)) {
  start_gps <- c(sample$start_long[i], sample$start_lat[i])
  end_gps <- c(sample$end_long[i], sample$end_lat[i])
  
  trn_gps <- rbind(start_gps, end_gps)
  trn_line <- st_sfc(st_linestring(trn_gps), crs = 4326)
  trn_lines[[i]] <- trn_line
}

trn_sf <- st_sf(geometry = do.call(c, trn_lines))
mapview(trn_sf)
