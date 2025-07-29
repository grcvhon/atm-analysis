# Ecological flow
# from: https://www.alexbaecher.com/post/connectivity-script/#lets-get-started
# modified 28 July 2025

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

states <- states()

se <- states %>%
  subset(REGION == "3") 

TN_NC <- se %>%     # Subsetting the data to Tennessee and North Carolina
  subset(NAME %in% c("Tennessee", "North Carolina"))

Pj <- occ(query = "Plethodon jordani",                 # JRCS scientific name
          from = "gbif",                               # limiting query to *the first* 1000 records
          limit=1000,                                  # limiting query to *the first* 1000 records
          has_coords = T)                              # limiting those 1000 records to those that have geo-referenced data

Pj_sp <- Pj$gbif$data$Plethodon_jordani %>%            # Grabbing the Darwin-core data from the spocc object
  dplyr::select(longitude,                             # Keep locations and year, discard the rest
                latitude,
                year) %>%                   
  dplyr::filter(year > 2000) %>%                       # Filter records to only those after year 2000
  filter(!duplicated(round(longitude, 2),              # Remove duplicate records using rounded decimals (this removes points very near to one-another)       
                     round(latitude, 2)) == TRUE) %>%  # >> See notes below about ^^
  dplyr::mutate(lon = scale(longitude),                # Remove points far outside the cluster of occurrences
                lat = scale(latitude)) %>%             # >> See notes below about ^^
  dplyr::filter(!abs(lon)>2) %>%
  dplyr::filter(!abs(lat)>2) %>%
  dplyr::select(longitude,
                latitude) 

Pj_sp <- SpatialPoints(coords = Pj_sp, proj4string = CRS("+init=epsg:4269"))


Pj_coords <- Pj_sp@coords                                                  
Pj_chull <- chull(Pj_sp@coords)                           # Creating convex hull

Pj_chull_ends <- Pj_sp@coords[c(Pj_chull, Pj_chull[1]),]  # generate the end points of polygon. 
Pj_poly <- SpatialPolygons(list(Polygons(list(Polygon(Pj_chull_ends)), ID=1)), proj4string = CRS("+init=epsg:4269"))       # convert coords to SpatialPolygons

Pj_poly_buff <- gBuffer(Pj_poly, width = 0.05, byid=T)

ggplot() + 
  #geom_polygon(data=TN_NC, aes(x=long, y=lat), col="grey40", fill="light blue") + 
  geom_polygon(data = Pj_poly_buff, aes(x = long, y = lat, group = group), col="grey40", fill="pink") +
  geom_point(data = as.data.frame(Pj_sp@coords), aes(x = longitude, y=latitude), size = 0.01) + 
  coord_quickmap() + 
  theme_map()


elevation <- get_elev_raster(sf::st_as_sf(Pj_poly_buff), z = 8)         # This will find a DEM tile nearest to our polygon

elv <- elevation %>% crop(Pj_poly_buff) %>% mask(Pj_poly_buff)

asp <- terrain(elv, opt = "aspect", neighbors = 8)

ggplot(as.data.frame(asp, xy=T)) + 
  geom_raster(aes(x=x, y=y, fill=aspect)) + 
  scale_fill_continuous(na.value=NA) + 
  theme_map() + 
  theme(legend.position = "right")

#set.seed(6)                                               # To make your results match mine
Pj_sample <- Pj_coords[sample(nrow(Pj_coords), 20),]       # Take 5 random locations

ggplot(as.data.frame(asp, xy=T)) + geom_raster(aes(x=x, y=y, fill=aspect)) + 
  geom_point(data=as.data.frame(Pj_sample), aes(x=longitude, y=latitude), size=2, col="white") +
  scale_fill_continuous(na.value=NA) + theme_map()

Pj_combn <- combn(nrow(Pj_sample),2) %>%
  t() %>%
  as.matrix()

asp_tr <- transition(asp, transitionFunction = mean, 4) %>%
  geoCorrection(type="c",multpl=F)

passages <- list()                                                     # Create a list to store the passage probability rasters in
system.time(                                                           # Keep track of how long this takes
  for (i in 1:nrow(Pj_combn)) {           
    locations <- SpatialPoints(rbind(Pj_sample[Pj_combn[i,1],1:2],     # create origin points
                                     Pj_sample[Pj_combn[i,2],1:2]),   # create destination (or goal) points, to traverse
                               proj4string = CRS("+init=epsg:4269"))
    passages[[i]] <- passage(asp_tr,                                   # run the passage function 
                             origin=locations[1],                 # set orgin point
                             goal=locations[2],                   # set goal point
                             theta = 0.00001)                             # set theta (tuning parameter, see notes below)
    print(paste((i/nrow(Pj_combn))*100, "% complete"))
  }
)

passages <- stack(passages)                                            # create a raster stack of all the passage probabilities
passages_overlay <- sum(passages)/nrow(Pj_combn)                       # calculate average

colors <- c("grey60", viridis_pal(option="plasma", begin = 0.3, end = 1)(20))
ggplot(as.data.frame(passages_overlay, xy=T)) + 
  geom_raster(aes(x=x,y=y,fill=layer)) +
  scale_fill_gradientn(colors = colors, na.value = NA) + 
  theme_map() +
  theme(legend.position = "right")

## Warning: Removed 18361 rows containing missing values (geom_raster).

###############################################################################

# load package
library(biooracler)

# BioOracle layers are at the spatial resolution of 0.05 x 0.05 decimal degrees and decadal temporal resolution

# list available layers
layers <- list_layers()
swd <- list_layers("SeaWaterDirection")

# check layer info
swd_id <- "swd_baseline_2000_2019_depthsurf"

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

# Load original raster
swd_lay_rasload <- swd_layer_raster 
# Create a clean copy by writing and re-reading it (strips NetCDF metadata)
swd_lay_ras <- raster(swd_lay_rasload) 
# write to a temporary GeoTIFF, then re-load it
tempfile <- tempfile(fileext = ".tif")
writeRaster(swd_lay_rasload, tempfile, formate = "GTiff", overwrite = TRUE)
# Read it back — now it's stripped of zvar, z-value, band info
swd_layer_raster <- raster(tempfile)
# re-name `names` attribute to "bearing"
names(swd_layer_raster) <- "bearing"
# Print
swd_layer_raster

# Sample coordinates
# - from all nw laevis samples, take only coordinates of the first 3 samples from listed pops
samp_coords <- laevis_nw %>% 
  filter(pop %in% c("Ashmore", "Broome", "Pilbara", "Exmouth_Gulf", "Shark_Bay")) %>% 
  group_by(pop) %>% slice(1:3) %>% ungroup()

samp_coords <- as.data.frame(samp_coords[,c(4,5)])

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

passages <- list()                                                     # Create a list to store the passage probability rasters in
system.time(                                                           # Keep track of how long this takes
  for (i in 1:nrow(laevis_combn)) {           
    locations <- SpatialPoints(rbind(samp_coords[laevis_combn[i,1],1:2],     # create origin points
                                     samp_coords[laevis_combn[i,2],1:2]),   # create destination (or goal) points, to traverse
                               proj4string = CRS("+init=epsg:4326"))
    passages[[i]] <- passage(swd_layer_tr,                                   # run the passage function 
                             origin=locations[1],                 # set orgin point
                             goal=locations[2],                   # set goal point
                             theta = 0.00001)                             # set theta (tuning parameter, see notes below)
    print(paste((i/nrow(laevis_combn))*100, "% complete"))
  }
)

passages <- stack(passages)                                            # create a raster stack of all the passage probabilities
passages_overlay <- sum(passages)/nrow(laevis_combn)                       # calculate average

colors <- c("grey50", viridis_pal(option="inferno", begin = 0.3, end = 1)(20))
ggplot(as.data.frame(passages_overlay, xy=T)) + 
  geom_raster(aes(x=x,y=y,fill=layer)) +
  scale_fill_gradientn(colors = colors, na.value = NA) + 
  #geom_point(data=as.data.frame(samp_coords), aes(x=longitude, y=latitude), size=1, col="red") +
  theme_map() +
  theme(legend.position = "right")
