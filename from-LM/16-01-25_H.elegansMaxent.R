## MaxEnt modelling framework ## 

# Call libraries 
library(sf)
library(leaflet)
library(readr)
library(janitor)
library(dplyr)
library(mapview)
library(spatstat)
library(tidyverse)
library(raster)
library(dismo)

library(tidyverse)
library(lubridate)
library(raster)
library(SDMtune)
library(spatstat)

## 1. Input model extent  - area around australia
eez <- st_read("Australia_EEZ/Australia_EEZ.shp")

#Crop eez shapefile to just Australia 
eez_crop <- st_crop(
  eez, 
  # x= long, y = lat - use google maps to get coordinates
  xmin = 106, xmax = 129, 
  ymin = -38, ymax = -10
)

leaflet() %>%
  addTiles() %>%
  addPolygons(data = eez_crop)

## 2. Input Occurrence data
H.elegans <- read_csv("2025-01-15_H.elegansoccur.csv") %>%
  mutate(
    lat = as.numeric(lat),
    long = as.numeric(long))%>%
  clean_names()

# 2.1 Filtering out missing coordinates 
H.elegans_clean <- filter(
  H.elegans,
  !is.na(lat),!is.na(long)
)

# 2.2. Filter out suspect points
H.elegans_clean2 <- H.elegans_clean %>% filter(suspect == 0)

# Identify duplicates within dataset
duplicates <- duplicated(H.elegans_clean2[, c("lat", "long")]) | duplicated(H.elegans_clean2[, c("lat", "long")], fromLast = TRUE)
view(duplicates)

#Colours for different occurrence datasets usede
dataset_colours <- c(
  "ALA" = "#d53e4f",
  "Dampier MP 2024" = "#f46d43",
  "Masterdataset" = "#fee08b",
  "RoebuckMP 2024-05" = "#abdda4",
  "Shark Bay 2024-06" = "#66c2a5",
  "Vinaydataset_2020" = "#3288bd"
)

mapview(
  H.elegans_clean2, 
  xcol = "long", ycol = "lat",
  crs = 4326,
  grid = FALSE, 
  zcol = "dataset", 
  col.regions = dataset_colours
)

# 2.3 Keep points witin eez shapefile
occ_sf <- H.elegans_clean2 %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  distinct(.keep_all = T) %>%
  st_crop(eez, )

mapview(occ_sf, xcol = "long", ycol = "lat", crs = 4326, grid = FALSE, 
        zcol = "dataset", col.regions = dataset_colours)
#Doesnt do much
  
## 4. Create bias layer 
# convert points and transect lines into a point pattern object (ppp or psp)

#This code doesnt work
occ_ppp <- 
  occ_sf %>% 
  st_transform(crs = 3577) %>% #Check whether the CRS needs to change 
  as.ppp() %>%


ggplot() + 
  annotation_map_tile('cartolight', zoom = 4) +
  layer_spatial(occ_sf, color = "red") +
  layer_spatial(eez, fill = NA, col = "black", lwd = 0.5) +
  annotation_scale(width_hint = 0.2, location = "br") +
  theme_void()

# 5. Calculate gausian density distribution of points 
pts_bias <-
  occ_ppp %>% 
  density(., sigma = 0.05) %>% 
  raster()
crs(pts_bias) <-  CRS("+init=epsg:3577")
pts_bias <-
  pts_bias %>% 
  projectRaster(., crs = CRS("+init=epsg:4326")) %>% 
  mask(mask = eez) %>% 
  resample(x = ., y = env[[1]]) # Need to input environmental rasters for this to work
values(pts_bias) <- values(pts_bias) + min(values(pts_bias), na.rm = T)

bias_layer <- pts_bias 
bias_layer[values(bias_layer) < 0] <- NA

plot(bias_layer)

## 6. Pseudo-absence generation -----------------------------------------------------------------------

bias.pts.bias <- 
  dismo::randomPoints(mask = bias_layer, n = 1000, p = as_Spatial(occ_sf), prob=TRUE) %>% 
  as_tibble() %>% 
  st_as_sf(coords=c("x","y"), crs=4326)

bias.pts.bg <-
  spsample(x = as_Spatial(eez), n = 500, type = "random") %>%
  st_as_sf() 

bias.pts <- rbind(bias.pts.bias, bias.pts.bg)

# Quick plot to check input data
mapview(bias_layer, na.color = "transparent", legend = F, homebutton = F, layer.name = "Bias Layer") +
  mapview(bias.pts, col.regions = "red", alpha = 0, legend = F, homebutton = F, layer.name = "Pseudo-absences") +
  mapview(occ_sf, col.regions = "white", alpha = 0, legend = F, homebutton = F, layer.name = "Occurrences")


#Ignore for now, will fix up later

# Generate points for just western Australia 

pts_bias <-
  occ_ppp %>% 
  density(., sigma = 0.05) %>% 
  raster()
crs(pts_bias) <-  CRS("+init=epsg:3577")
pts_bias <-
  pts_bias %>% 
  projectRaster(., crs = CRS("+init=epsg:4326")) %>% 
  mask(mask = eez_crop) %>% 
  resample(x = ., y = env[[1]])
values(pts_bias) <- values(pts_bias) + min(values(pts_bias), na.rm = T)

bias_layer <- pts_bias 
bias_layer[values(bias_layer) < 0] <- NA

plot(bias_layer)

#Pseudo-absence generation -----------------------------------------------------------------------

bias.pts.bias <- 
  dismo::randomPoints(mask = bias_layer, n = 1000, p = as_Spatial(occ_sf), prob=TRUE) %>% 
  as_tibble() %>% 
  st_as_sf(coords=c("x","y"), crs=4326)

bias.pts.bg <-
  spsample(x = as_Spatial(eez_crop), n = 500, type = "random") %>%
  st_as_sf() 

bias.pts <- rbind(bias.pts.bias, bias.pts.bg)

# Quick plot to check input data
mapview(bias_layer, na.color = "transparent", legend = F, homebutton = F, layer.name = "Bias Layer") +
  mapview(bias.pts, col.regions = "red", alpha = 0, legend = F, homebutton = F, layer.name = "Pseudo-absences") +
  mapview(occ_sf, col.regions = "white", alpha = 0, legend = F, homebutton = F, layer.name = "Occurrences")

# Input environmental layers:
# install.packages("devtools")
install.packages("devtools")
devtools::install_github("bio-oracle/biooracler")

library(biooracler)
list_layers()
