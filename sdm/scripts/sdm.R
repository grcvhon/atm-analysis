# Species distribution modelling

# setwd: sdm directory
setwd("C:/Users/a1235304/Dropbox/Short-nosed and Leaf-scaled sea snake TSSC/atm-analysis/sdm")

# install/load packages
packages <- c("sf","leaflet","readr","janitor","dplyr",
              "mapview","spatstat","tidyverse","raster",
              "dismo","lubridate","SDMtune","readxl",
              "terra","stars","lwgeom","maptools","ggspatial",
              "prettymapr","tidyterra")

# load package list all at once
invisible(lapply(packages, library, character.only = TRUE))

#### 1. Input model extent ####

# Load nwshelf shapefile
nw_shelf <- st_read("./nw-shelf/NWShelf.shp", quiet = TRUE) %>% 
  st_transform(4326)

# preview nwshelf
plot(nw_shelf)

#### 2. Input environmental rasters ####

# Load single environmental raster
bathymetry <- raster("./predictor-variables/bathymetry.asc")

# Load initial set of environmental rasters
env_init <- stack("./predictor-variables/sal_mean.asc",
                  "./predictor-variables/sal_amp.asc",
                  "./predictor-variables/bathymetry.asc",
                  "./predictor-variables/sst_mean.asc",
                  "./predictor-variables/sst_amp.asc",
                  "./predictor-variables/chlor_mean.asc",
                  "./predictor-variables/DistToLand.asc",
                  "./predictor-variables/DistToReef.asc",
                  "./predictor-variables/DistToFW.asc")

#### 3. Input occurrence data ####

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



#### 4. Bias layer ####

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
  geom_sf(data = nw_shelf, fill = NA)
