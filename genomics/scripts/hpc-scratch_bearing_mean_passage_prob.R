# Generating mean passage probability across the northwest shelf
# based on ocean bearing (direction)

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

sessionInfo()






### *** BioOracle layers: Bearing (direction) and Speed *** ###

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
writeRaster(swd_layer_raster, filename = "swd_layer_raster.tif")

# write to a temporary GeoTIFF, then re-load it
tempfile <- tempfile(fileext = ".tif")
writeRaster(sws_lay_rasload, tempfile, formate = "GTiff", overwrite = TRUE)
# Read it back — now it's stripped of zvar, z-value, band info
sws_layer_raster <- raster(tempfile)
# re-name `names` attribute to "bearing"
names(sws_layer_raster) <- "speed"
# Print
sws_layer_raster
writeRaster(sws_layer_raster, filename = "sws_layer_raster.tif")

speed_rast <- rast("./genomics/sws_layer_raster.tif")





### *** crop bearing and speed rasters based on NW shelf shapefile extent *** ###

# load shapefile
nw_shelf <- st_read("/hpcfs/users/a1235304/atm_passage/shapefiles/nw_shelf/NWShelf.shp", quiet = TRUE) %>% st_transform(4326)

# crop bearing layer
swd_layer_raster <- crop(swd_layer_raster, extent(nw_shelf))
swd_layer_raster <- mask(swd_layer_raster, nw_shelf)
plot(swd_layer_raster)

# crop speed layer
sws_layer_raster <- crop(sws_layer_raster, extent(nw_shelf))
sws_layer_raster <- mask(sws_layer_raster, nw_shelf)
plot(sws_layer_raster)






### *** Generate regular points across the extent *** ###

# regular

set.seed(1) # generates exactly 1000 points
size <- 50
reg_points <- st_sample(nw_shelf, 
                        size = size, 
                        type = "regular", 
                        replace = FALSE
                        # exact = TRUE,
                        ) # even the `exact` arg does not give exact number of points
plot(reg_points, pch = 19, cex = 0.1)
df_reg_points <- sfheaders::sfc_to_df(reg_points)
df_reg_points <- df_reg_points[,-c(1,2)]
colnames(df_reg_points) <- c("longitude", "latitude")
print(dim(df_reg_points)) # check number of rows and columns
print(dim(unique(df_reg_points))) # check unique number of rows and columns
print(dim(unique(df_reg_points)) == dim(df_reg_points)) # should be TRUE TRUE

# generate origin/destination combinations across the 1000 points
reg_pts_combn <- 
  combn(nrow(df_reg_points),2) %>%
  t() %>%
  as.matrix()
reg_pts_combn






### *** Generate transition layer for bearing *** ###

# bearing
swd_layer_tr <- transition(swd_layer_raster, transitionFunction = mean, directions = 8) %>% 
  geoCorrection(type = "c", multpl = F)

#ggplot(as.data.frame(swd_layer_raster, xy=T)) + 
#  geom_raster(aes(x=x, y=y, fill = bearing)) + 
#  geom_point(data=as.data.frame(samp_coords), aes(x=longitude, y=latitude), size=2, col="red") +
#  scale_fill_continuous(na.value=NA) + 
#  theme_map() + theme(legend.position = "right")






### *** Generate passage probability layer based on bearing *** ###

# Create a list to store the passage probability rasters in
reg_bearing_passages <- list()                                                     

system.time( # Keep track of how long this takes
  for (i in 1:nrow(reg_pts_combn)) {           
    locations <- SpatialPoints(rbind(df_reg_points[reg_pts_combn[i,1],1:2],   # create origin points
                                     df_reg_points[reg_pts_combn[i,2],1:2]),  # create destination (or goal) points, to traverse
                               proj4string = CRS("+init=epsg:4326"))
    reg_bearing_passages[[i]] <- passage(swd_layer_tr,                        # run the passage function 
                                         origin=locations[1],                 # set orgin point
                                         goal=locations[2],                   # set goal point
                                         theta = 0.00001)                     # set theta (tuning parameter, see notes below)
    print(paste((i/nrow(reg_pts_combn))*100, "% complete"))
  }
)

reg_bearing_passages <- stack(reg_bearing_passages) # create a raster stack of all the passage probabilities
reg_bearing_passages_overlay <- sum(reg_bearing_passages)/nrow(reg_pts_combn) # calculate average






### *** save output *** ###

dir <- paste0("/hpcfs/users/a1235304/atm_passage/output/regular_", size, "pts/")
dir.create(dir)

# as csv
library(terra)
ovr <- rast(random_bearing_passages_overlay)
df_ovr <- as.data.frame(ovr, xy = TRUE)

library(utils)
write.csv(df_ovr, "/hpcfs/users/a1235304/atm_passage/output/regular_", size ,"pts.csv")

# as pdf
pdf(paste0("/hpcfs/users/a1235304/atm_passage/output/regular_", size ,"pts.pdf"), height = 8.5, width = 11)
    
colors <- c("grey50", viridis_pal(option="inferno", begin = 0.3, end = 1)(20))
ggplot(as.data.frame(reg_bearing_passages_overlay, xy=T)) + 
  geom_raster(aes(x=x,y=y,fill=layer)) +
  scale_fill_gradientn(colors = colors, na.value = NA) + 
  #geom_point(data=as.data.frame(samp_coords), aes(x=longitude, y=latitude), size=1, col="red") +
  theme_map() +
  theme(legend.position = "right")

dev.off()

print(quote = FALSE, paste0("Output saved in ", dir,"."))
print(quote = FALSE, " ")

timestamp()

