# Generating mean passage probability across the northwest shelf
# based on ocean current bearing (direction)
# using spatially balanced points

timestamp()
print(quote = FALSE, " ")

library(gdistance)
library(tidyverse)
library(ggplot2)
library(raster)
library(viridis)
library(ggthemes)
library(biooracler)
library(sf)
library(terra)

print(quote=FALSE, "Packages loaded!")


### ***        Bearing (direction)        *** ###
### *** lines for speed are commented out *** ###

# BioOracle layers are at the spatial resolution of 0.05 x 0.05 decimal degrees and 
# decadal temporal resolution

# list available layers
#layers <- list_layers()
#swd <- list_layers("SeaWaterDirection")
#sws <- list_layers("SeaWaterSpeed")

# check layer info
#swd_id <- "swd_baseline_2000_2019_depthsurf"
#sws_id <- "sws_baseline_2000_2019_depthsurf"

# time based on info_layer
#oracle_time <- c("2000-01-01T00:00:00Z", "2010-01-01T00:00:00Z")
# extent based on "nw_shelf" shape i.e., ext(nw_shelf)
#oracle_lat <- c(-26.7363816780876, -9.69138797000501)
#oracle_lon <- c(111.544995584098, 130.342214992267)

#constraints <- list(oracle_time, oracle_lat, oracle_lon)
#names(constraints) <- c("time", "latitude", "longitude")

#swd_var <- "swd_mean" # sea water direction - mean
#swd_layer <- download_layers(swd_id, swd_var, constraints)
#swd_layer_raster <- raster::raster(swd_layer$swd_mean_1)

#sws_var <- "sws_mean" # sea water speed - mean
#sws_layer <- download_layers(sws_id, sws_var, constraints)
#sws_layer_raster <- raster::raster(sws_layer$sws_mean_1)

# Load original raster
#swd_lay_rasload <- swd_layer_raster
#sws_lay_rasload <- sws_layer_raster
# Create a clean copy by writing and re-reading it (strips NetCDF metadata)
#swd_lay_ras <- raster(swd_lay_rasload)
#sws_lay_ras <- raster(sws_lay_rasload)

# write to a temporary GeoTIFF, then re-load it
#tempfile <- tempfile(fileext = ".tif")
#writeRaster(swd_lay_rasload, tempfile, formate = "GTiff", overwrite = TRUE)
# Read it back — now it's stripped of zvar, z-value, band info
#swd_layer_raster <- raster(tempfile)
# re-name `names` attribute to "bearing"
#names(swd_layer_raster) <- "bearing"
# Print
#swd_layer_raster
#writeRaster(swd_layer_raster, filename = "swd_layer_raster.tif")

# write to a temporary GeoTIFF, then re-load it
#tempfile <- tempfile(fileext = ".tif")
#writeRaster(sws_lay_rasload, tempfile, formate = "GTiff", overwrite = TRUE)
# Read it back — now it's stripped of zvar, z-value, band info
#sws_layer_raster <- raster(tempfile)
# re-name `names` attribute to "bearing"
#names(sws_layer_raster) <- "speed"
# Print
#sws_layer_raster
#writeRaster(sws_layer_raster, filename = "sws_layer_raster.tif")

### To circumvent cURL timeout:
### *** The raster files have been downloaded and now exist in rasterfiles/ dir *** ###

# Load the tif file (class will be SpatRaster)
swd_layer_raster <- rast("/hpcfs/users/a1235304/atm_passage/rasterfiles/swd_layer_raster.tif")
#sws_layer_raster <- rast("/hpcfs/users/a1235304/atm_passage/rasterfiles/sws_layer_raster.tif")

# Convert to RasterLayer
swd_layer_raster <- raster(swd_layer_raster)
#sws_layer_raster <- raster(sws_layer_raster)



### *** crop bearing and speed rasters based on NW shelf shapefile extent *** ###

# load shapefile
nw_shelf <- st_read("/hpcfs/users/a1235304/atm_passage/shapefiles/nw_shelf/NWShelf.shp", quiet = TRUE) %>% st_transform(4326)

# crop bearing layer
#swd_layer_raster <- crop(swd_layer_raster, extent(nw_shelf))
#swd_layer_raster <- mask(swd_layer_raster, nw_shelf)
#plot(swd_layer_raster)

# crop speed layer
#sws_layer_raster <- crop(sws_layer_raster, extent(nw_shelf))
#sws_layer_raster <- mask(sws_layer_raster, nw_shelf)
#plot(sws_layer_raster)

library(sf)
nw_shelf_utm <- st_transform(nw_shelf, crs = 32750)

library(dssduoa)
nw_shelf_region <- make.region(shape = nw_shelf_utm)

#ecoflow_points <- st_as_sfc(st_bbox(nw_shelf_utm)) #|>
#st_make_grid(cellsize = 1000) |>
#st_as_sf()

#ecoflow_points_inPoly <- st_intersection(ecoflow_points, nw_shelf_utm)

library(spsurvey)
seed <- 100
set.seed(seed)
n_base <- 100
ecoflow_pts <- grts(nw_shelf_utm, n_base = n_base)

#plot(ecoflow_pts)

ecoflow_pts_sf <- st_as_sf(ecoflow_pts$sites_base)

#library(ggplot2)
#ggplot() + geom_sf(data = nw_shelf_utm) + geom_sf(data = ecoflow_pts_sf)

ecoflow_pts_lon <- ecoflow_pts$sites_base$lon_WGS84
ecoflow_pts_lat <- ecoflow_pts$sites_base$lat_WGS84
ecoflow_pts_coords <- cbind(ecoflow_pts_lon,ecoflow_pts_lat)
colnames(ecoflow_pts_coords) <- c("longitude","latitude")
df_ecoflow_pts_coords <- as.data.frame(ecoflow_pts_coords)

# introduce manually selected points in EG and SB
manual_pts <- 
  data.frame(longitude = c(114.37139,114.31644,114.29309,113.69914,113.25505,113.45073),
             latitude = c(-21.86184,-22.08711,-22.31837,-26.38297,-25.73358,-25.13991))

# append manually selected points to sbs generated pts
df_ecoflow_pts_coords <- rbind(df_ecoflow_pts_coords, manual_pts)

# check for duplicate points
dim(df_ecoflow_pts_coords) == dim(unique(df_ecoflow_pts_coords)) # no duplicate points

# generate origin/destination combinations across the 1000 points
ecoflow_pts_comb <- 
  combn(nrow(df_ecoflow_pts_coords),2) %>%
  t() %>%
  as.matrix()
#ecoflow_pts_comb

# bearing
swd_layer_tr <- transition(swd_layer_raster, transitionFunction = mean, directions = 8) %>% 
  geoCorrection(type = "c", multpl = F)

ecoflow_bearing_passages <- list()                                                     

system.time( # Keep track of how long this takes
  for (i in 1:nrow(ecoflow_pts_comb)) {           
    locations <- SpatialPoints(rbind(df_ecoflow_pts_coords[ecoflow_pts_comb[i,1],1:2],   # create origin points
                                     df_ecoflow_pts_coords[ecoflow_pts_comb[i,2],1:2]),  # create destination (or goal) points, to traverse
                               proj4string = CRS("+init=epsg:4326"))
    ecoflow_bearing_passages[[i]] <- passage(swd_layer_tr,                        # run the passage function 
                                             origin=locations[1],                 # set orgin point
                                             goal=locations[2],                   # set goal point
                                             theta = 0.00001)                     # set theta (tuning parameter, see notes below)
    print(paste((i/nrow(ecoflow_pts_comb))*100, "% complete"))
  }
)

ecoflow_bearing_passages <- stack(ecoflow_bearing_passages) # create a raster stack of all the passage probabilities
ecoflow_bearing_passages_overlay <- sum(ecoflow_bearing_passages)/nrow(ecoflow_pts_comb) # calculate average



### *** save output *** ###

dir <- paste0("/hpcfs/users/a1235304/atm_passage/output/sbs_bearing_seed",seed,"_",n_base,"pts_",format(Sys.time(),"%Hh%Mm%Ss"),"/")
dir.create(dir, recursive = TRUE)

# as csv
library(terra)
ovr <- rast(ecoflow_bearing_passages_overlay)
df_ovr <- as.data.frame(ovr, xy = TRUE)

library(utils)
write.csv(df_ovr, file = paste0(dir, "sbs_bearing_seed",seed,"_",n_base,"pts_",format(Sys.time(),"%Hh%Mm%Ss"),".csv"))
          
# as pdf

pdf(file = paste0(dir, "sbs_bearing_seed",seed,"_",n_base,"pts_",format(Sys.time(),"%Hh%Mm%Ss"),".pdf"), height = 8.5, width = 11)

colors <- c("grey50", viridis_pal(option="inferno", begin = 0.3, end = 1)(20))
ggplot(as.data.frame(ecoflow_bearing_passages_overlay, xy=T)) +
  geom_raster(aes(x=x,y=y,fill=layer)) +
  scale_fill_gradientn(colors = colors, na.value = NA) +
  #geom_point(data=as.data.frame(samp_coords), aes(x=longitude, y=latitude), size=1, col="red") +
  theme_map() +
  theme(legend.position = "right")

dev.off()

print(quote = FALSE, paste0("Output saved in ", dir,"."))
print(quote = FALSE, " ")

sessionInfo()
print(quote=FALSE, " ")


timestamp()
