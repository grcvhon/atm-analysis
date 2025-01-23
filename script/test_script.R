### MaxEnt modelling framework ###
### 2025-01-16 ###

### i. Install packages ----------------------------------------------------------
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

### ii. Load package list ---------------------------------------------------------
packages <- c("sf","leaflet","readr","janitor","dplyr",
              "mapview","spatstat","tidyverse","raster",
              "dismo","lubridate","SDMtune","readxl",
              "terra","stars","lwgeom","maptools","ggspatial",
              "prettymapr")
invisible(lapply(packages, library, character.only = TRUE))

### 1. Input model extent -----------------------------------------------------

# Load Northwest shelf bounding area
nw_shelf <- st_read("data/shapefiles/nw-shelf/NWShelf.shp")
mapview(nw_shelf)

### 2. Input environmental predictors -----------------------------------------

# Load initial set of environmental rasters
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

bathymetry <- raster("data/predictor-variables/bathymetry.asc")
mapview(bathymetry) + nw_shelf

### 3. Input occurrence data --------------------------------------------------

# Load occurrence data csv contains combined occurrence data for 
# Aipysurus apraefrontalis and A. foliosquama
rawdat <- read_csv("data/spreadsheets/ATM_2023_0715-running-master.csv")

# function to filter by species
filter_species <- function(species_name) {
  rawdat %>%
    mutate(
      lat = as.numeric(MidLat),
      lat = if_else(is.na(lat), as.numeric(Latitude), lat),
      long = as.numeric(MidLong),
      long = if_else(is.na(long), as.numeric(Longitude), long)
    ) %>%
    filter(Species == species_name,
           !is.na(lat),
           !is.na(long)) %>%
    clean_names()
}

# create data frame filtered for each species
aprae <- filter_species("apraefrontalis")
folio <- filter_species("foliosquama")

# function to convert dataframe to simple feature (sf) 
convert_2_sf <- function(df_name) {
  df_name %>% 
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
    distinct(.keep_all = T)
}

# convert dataframe to simple feature (sf)
aprae_sf <- convert_2_sf(aprae)
folio_sf <- convert_2_sf(folio)
seasnake_sf <- bind_rows(aprae_sf,folio_sf) # combine aprae and folio sf occurrence points

### 3.1 Input transect data
# shapefiles generated with "Generate and save shapefile 
# using transect start and end coordinates.R"
aprae_transect <- st_read("data/shapefiles/apraefrontalis_transect.shp")
folio_transect <- st_read("data/shapefiles/foliosquama_transect.shp")
seasnake_trn <- st_read("data/shapefiles/seasnake_transect.shp") # combined aprae and folio transect lines

### 4. Create Bias Layer ------------------------------------------------------
# convert points and transect lines into a point pattern object (ppp or psp)

# using combined apraefrontalis and foliosquama occurrence points
seasnake_occ_ppp <- seasnake_sf %>% 
  st_transform(crs = 3577) %>% 
  as_Spatial() %>% 
  maptools::as.ppp.SpatialPointsDataFrame(.)

# using combined apraefrontalis and foliosquama transect lines
seasnake_trn_psp <- seasnake_trn %>% 
  st_transform(crs = 3577) %>% 
  as_Spatial() %>% 
  maptools::as.psp.SpatialLinesDataFrame(.)

# calculate Gaussian density distribution of points and transects
seasnake_pts_bias <- seasnake_occ_ppp %>% 
  density(., sigma = 0.05) %>% 
  raster()
crs(seasnake_pts_bias) <- CRS(SRS_string = "EPSG:3577")
seasnake_pts_bias <- seasnake_pts_bias %>% 
  projectRaster(., crs = CRS(SRS_string = "EPSG:4326")) %>% 
  mask(mask = nw_shelf) %>% 
  resample(x = ., y = bathymetry)
values(seasnake_pts_bias) <- values(seasnake_pts_bias) + min(values(seasnake_pts_bias), na.rm = T)

seasnake_trn_bias <- seasnake_trn_psp %>% 
  density(., sigma = 0.05) %>% 
  raster()
crs(seasnake_trn_bias) <- CRS(SRS_string = "EPSG:3577")
seasnake_trn_bias <- seasnake_trn_bias %>% 
  projectRaster(., crs = CRS(SRS_string = "EPSG:4326")) %>% 
  mask(mask = nw_shelf) %>% 
  resample(x = ., y = bathymetry)
values(seasnake_trn_bias) <- values(seasnake_trn_bias) + min(values(seasnake_trn_bias), na.rm = T)

mapview(seasnake_pts_bias)
mapview(seasnake_trn_bias)

seasnake_bias_layer <- seasnake_pts_bias + seasnake_trn_bias
seasnake_bias_layer[values(seasnake_bias_layer) < 0] <- NA
mapview(seasnake_bias_layer)


# species-specific occurrence map + bounding area
mapview(aprae, xcol = "long", ycol = "lat", crs = 4326) + nw_shelf  
mapview(folio, xcol = "long", ycol = "lat", crs = 4326) + nw_shelf
