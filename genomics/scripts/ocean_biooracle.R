### try BioOracle instead...
library(biooracler)

# BioOracle layers are at the spatial resolution of 0.05 x 0.05 decimal degrees and decadal temporal resolution

# list available layers
layers <- list_layers()
swd <- list_layers("SeaWaterDirection")
#View(swd)
# want swd_baseline_2000_2019_depthsurf | Bio-Oracle SeaWaterDirection [depthSurf]Baseline 2000-2019.

sws <- list_layers("SeaWaterSpeed")
#View(sws)
# want sws_baseline_2000_2019_depthsurf | Bio-Oracle SeaWaterSpeed [depthSurf]Baseline 2000-2019

# check layer info
swd_id <- "swd_baseline_2000_2019_depthsurf"
#info_layer(swd_id)

sws_id <- "sws_baseline_2000_2019_depthsurf"
#info_layer(sws_id)

# set constraints before downloading data

# time based on info_layer
oracle_time <- c("2000-01-01T00:00:00Z", "2010-01-01T00:00:00Z")
# extent based on "nw_shelf" shape i.e., ext(nw_shelf)
oracle_lat <- c(-26.7363816780876, -9.69138797000501)
oracle_lon <- c(111.544995584098, 130.342214992267)

constraints <- list(oracle_time, oracle_lat, oracle_lon)
names(constraints) <- c("time", "latitude", "longitude")

# define wanted variables to download
swd_var <- "swd_mean" # sea water direction - mean
sws_var <- "sws_mean" # sea water speed - mean

# perform download
swd_layer <- download_layers(swd_id, swd_var, constraints)
sws_layer <- download_layers(sws_id, sws_var, constraints)

### side quest: marmap adaptation
  #swd_rast <- raster::raster(swd_layer)
  #plot(swd_rast)
  #swd_bathy <- marmap::as.bathy(swd_rast)
  #plot(swd_bathy, image = TRUE)
  
  #library(raster)
  #library(terra)
  #library(marmap)
  swd_mask <- raster::crop(swd_layer, nw_shelf)
  swd_mask <- terra::mask(swd_mask, mask = terra::vect(nw_shelf))
  swd_mask <- raster::raster(swd_mask)
  plot(swd_mask)
  swd_mask_bathy <- marmap::as.bathy(swd_mask)
  summary(swd_mask_bathy)
  swd_trans <- 
  swd_lcdist <- marmap::lc.dist(swd_trans, sub_laevis, res = "path")
  lapply(nw_lcdist200,lines,col="dodgerblue",lwd=3,lty=1)
  
  #library(rWind)
  # working with swd_layer downloaded from biooracle
  swd_layer_coords <- terra::crds(swd_layer$swd_mean_1,  df = TRUE) # get coords from swd_mean_1 layer and into a df 
  swd_layer_dr1_df <- as.data.frame(swd_layer$swd_mean_1) # get mean direction  1
  sws_layer_sp1_df <- as.data.frame(sws_layer$sws_mean_1) # get mean speed 1
  ocean_uv <- ds2uv(ocean_c$dir, ocean_c$speed) # convert direction and speed to eastward and northward vectors
  ocean_c <- cbind(swd_layer_coords,swd_layer_dr1_df,sws_layer_sp1_df) # combine columns
  colnames(ocean_c) <- c("lon","lat","dir","speed") # rename headers
  class(ocean_c) <- c("rWind", "data.frame") # reclass
  ggplot(ocean_c, aes(x = lon, y = lat, color = dir)) + geom_point()
  
  ocean_c_ras <- wind2raster(ocean_c)
  ocean_conduct <- flow.dispersion(ocean_c_ras)
  GtoS <- shortestPath(ocean_conduct, 
                       goal = c(sub_laevis[2,1], sub_laevis[2,2]),
                       origin = c(sub_laevis[1,1], sub_laevis[1,2]),
                       output = "SpatialLines")
  
  StoG <- shortestPath(ocean_conduct, 
                       origin = c(sub_laevis[2,1], sub_laevis[2,2]),
                       goal = c(sub_laevis[1,1], sub_laevis[1,2]),
                       output = "SpatialLines")

  
  #example
  data(wind.data)
  wd_ras <- wind2raster(wind.data)
  conductance <- flow.dispersion(wd_ras)
  AtoB <- shortestPath(conductance, c(-5.5,37), c(-5.5,35), output = "SpatialLines")
  BtoA <- shortestPath(conductance, c(-5.5,35), c(-5.5,37), output = "SpatialLines")
  
# recall: nw_shelf - object containing shapefile of extent
nw_swd <- raster::crop(swd_layer, nw_shelf)
nw_swd <- terra::mask(nw_swd, mask = terra::vect(nw_shelf))
mapview::mapview(nw_swd, na.color = NA)

nw_sws <- raster::crop(sws_layer, nw_shelf)
nw_sws <- terra::mask(nw_sws, mask = terra::vect(nw_shelf))
mapview::mapview(nw_sws, na.color = NA)


###

library(marmap)
library(raster)
library(gdistance)

# get bathymetric data for Australia
nw_bathy <- getNOAA.bathy(lon1=110.544995584098, lon2=118.342214992267, 
                           lat1=-18.69138797000501, lat2=-29.7363816780876,
                           resolution=1, keep=TRUE)

# subset laevis sites
sub_laevis <- laevis_nw[c(122,170),c(4,5)] # Exmouth and Shark Bay
sub_laevis

nw_trans <- trans.mat(nw_bathy)
nw_trans200 <- trans.mat(nw_bathy, min.depth = 0, max.depth = -200)

nw_lcdist <- lc.dist(nw_trans,sub_laevis,res="path")
nw_lcdist200 <- lc.dist(nw_trans200,sub_laevis,res="path")

plot(nw_bathy, image=TRUE)
lapply(nw_lcdist200,lines,col="orange",lwd=3,lty=1)
points(sub_laevis, pch = 19, col = "red", cex = 1.5)

###









