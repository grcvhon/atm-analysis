# # Ecological flow
# from: https://www.alexbaecher.com/post/connectivity-script/#lets-get-started
# modified 29 July 2025

library(gdistance)
library(tidyverse)
library(rgeos)
library(elevatr)
library(ggplot2)
library(tigris)
library(spocc)
library(raster)
library(viridis)
library(ggthemes)
library(tidyverse)
library(biooracler)

# BioOracle layers are at the spatial resolution of 0.05 x 0.05 decimal degrees and decadal temporal resolution

# list available layers
layers <- list_layers()
swd <- list_layers("SeaWaterDirection")
sws <- list_layers("SeaWaterSpeed")

# check layer info
swd_id <- "swd_baseline_2000_2019_depthsurf"
sws_id <- "sws_baseline_2000_2019_depthsurf"

# time based on info_layer
oracle_time <- c("2000-01-01T00:00:00Z", "2010-01-01T00:00:00Z")
# extent based on "nw_shelf" shape i.e., ext(nw_shelf)
oracle_lat <- c(-26.7363816780876, -9.69138797000501)
oracle_lon <- c(111.544995584098, 130.342214992267)

constraints <- list(oracle_time, oracle_lat, oracle_lon)
names(constraints) <- c("time", "latitude", "longitude")

swd_var <- "swd_mean" # sea water direction - mean
swd_layer <- download_layers(swd_id, swd_var, constraints)
swd_layer_raster <- raster::raster(swd_layer$swd_mean_1)

sws_var <- "sws_mean" # sea water speed - mean
sws_layer <- download_layers(sws_id, sws_var, constraints)
sws_layer_raster <- raster::raster(sws_layer$sws_mean_1)

# Load original raster
swd_lay_rasload <- swd_layer_raster
sws_lay_rasload <- sws_layer_raster
# Create a clean copy by writing and re-reading it (strips NetCDF metadata)
swd_lay_ras <- raster(swd_lay_rasload)
sws_lay_ras <- raster(sws_lay_rasload)

# write to a temporary GeoTIFF, then re-load it
tempfile <- tempfile(fileext = ".tif")
writeRaster(swd_lay_rasload, tempfile, formate = "GTiff", overwrite = TRUE)
# Read it back — now it's stripped of zvar, z-value, band info
swd_layer_raster <- raster(tempfile)
# re-name `names` attribute to "bearing"
names(swd_layer_raster) <- "bearing"
# Print
swd_layer_raster

# write to a temporary GeoTIFF, then re-load it
tempfile <- tempfile(fileext = ".tif")
writeRaster(sws_lay_rasload, tempfile, formate = "GTiff", overwrite = TRUE)
# Read it back — now it's stripped of zvar, z-value, band info
sws_layer_raster <- raster(tempfile)
# re-name `names` attribute to "bearing"
names(sws_layer_raster) <- "speed"
# Print
sws_layer_raster

# Sample coordinates
# - from all nw laevis samples, take only coordinates of the first 3 samples from listed pops
samp_coords <- laevis_nw %>% 
  filter(pop %in% c("Ashmore", "Broome", "Exmouth_Gulf", "Pilbara", "Shark_Bay")) %>% 
  group_by(pop) %>% slice(1:3) %>% ungroup()

# keep only longitude and latitude columns
samp_coords <- as.data.frame(samp_coords[,c(4,5)])

# bearing

# generate transition layer
swd_layer_tr <- transition(swd_layer_raster, transitionFunction = mean, directions = 8) %>% 
  geoCorrection(type = "c", multpl = F)

ggplot(as.data.frame(swd_layer_raster, xy=T)) + 
  geom_raster(aes(x=x, y=y, fill = bearing)) + 
  geom_point(data=as.data.frame(samp_coords), aes(x=longitude, y=latitude), size=2, col="red") +
  scale_fill_continuous(na.value=NA) + 
  theme_map() + theme(legend.position = "right")

laevis_combn <- 
  combn(nrow(samp_coords),2) %>%
  t() %>%
  as.matrix()

bearing_passages <- list()                                                     # Create a list to store the passage probability rasters in
system.time(                                                           # Keep track of how long this takes
  for (i in 1:nrow(laevis_combn)) {           
    locations <- SpatialPoints(rbind(samp_coords[laevis_combn[i,1],1:2],     # create origin points
                                     samp_coords[laevis_combn[i,2],1:2]),   # create destination (or goal) points, to traverse
                               proj4string = CRS("+init=epsg:4326"))
    bearing_passages[[i]] <- passage(swd_layer_tr,                                   # run the passage function 
                             origin=locations[1],                 # set orgin point
                             goal=locations[2],                   # set goal point
                             theta = 0.00001)                             # set theta (tuning parameter, see notes below)
    print(paste((i/nrow(laevis_combn))*100, "% complete"))
  }
)

bearing_passages <- stack(bearing_passages)                                            # create a raster stack of all the passage probabilities
bearing_passages_overlay <- sum(bearing_passages)/nrow(laevis_combn)                       # calculate average

colors <- c("grey50", viridis_pal(option="inferno", begin = 0.3, end = 1)(20))
ggplot(as.data.frame(bearing_passages_overlay, xy=T)) + 
  geom_raster(aes(x=x,y=y,fill=layer)) +
  scale_fill_gradientn(colors = colors, na.value = NA) + 
  #geom_point(data=as.data.frame(samp_coords), aes(x=longitude, y=latitude), size=1, col="red") +
  theme_map() +
  theme(legend.position = "right")

# speed

# generate transition layer
sws_layer_tr <- transition(sws_layer_raster, transitionFunction = mean, directions = 8) %>% 
  geoCorrection(type = "c", multpl = F)

ggplot(as.data.frame(sws_layer_raster, xy=T)) + 
  geom_raster(aes(x=x, y=y, fill = speed)) + 
  geom_point(data=as.data.frame(samp_coords), aes(x=longitude, y=latitude), size=2, col="red") +
  scale_fill_continuous(na.value=NA) + 
  theme_map() + theme(legend.position = "right")

laevis_combn <- 
  combn(nrow(samp_coords),2) %>%
  t() %>%
  as.matrix()

speed_passages <- list()                                                     # Create a list to store the passage probability rasters in
system.time(                                                           # Keep track of how long this takes
  for (i in 1:nrow(laevis_combn)) {           
    locations <- SpatialPoints(rbind(samp_coords[laevis_combn[i,1],1:2],     # create origin points
                                     samp_coords[laevis_combn[i,2],1:2]),   # create destination (or goal) points, to traverse
                               proj4string = CRS("+init=epsg:4326"))
    speed_passages[[i]] <- passage(sws_layer_tr,                                   # run the passage function 
                             origin=locations[1],                 # set orgin point
                             goal=locations[2],                   # set goal point
                             theta = 0.00001)                             # set theta (tuning parameter, see notes below)
    print(paste((i/nrow(laevis_combn))*100, "% complete"))
  }
)

speed_passages <- stack(speed_passages)                                            # create a raster stack of all the passage probabilities
speed_passages_overlay <- sum(speed_passages)/nrow(laevis_combn)                       # calculate average

colors <- c("grey50", viridis_pal(option="inferno", begin = 0.3, end = 1)(20))
ggplot(as.data.frame(speed_passages_overlay, xy=T)) + 
  geom_raster(aes(x=x,y=y,fill=layer)) +
  scale_fill_gradientn(colors = colors, na.value = NA) + 
  geom_point(data=as.data.frame(samp_coords), aes(x=longitude, y=latitude), size=2, col="red") +
  theme_map() +
  theme(legend.position = "right")
