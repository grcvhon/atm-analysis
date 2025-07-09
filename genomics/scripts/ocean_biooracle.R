### try BioOracle instead...
library(biooracler)

# BioOracle layers are at the spatial resolution of 0.05 x 0.05 decimal degrees and decadal temporal resolution

# list available layers
layers <- list_layers()
swd <- list_layers("SeaWaterDirection")
View(swd)
# want swd_baseline_2000_2019_depthsurf | Bio-Oracle SeaWaterDirection [depthSurf]Baseline 2000-2019.

sws <- list_layers("SeaWaterSpeed")
View(sws)
# want sws_baseline_2000_2019_depthsurf | Bio-Oracle SeaWaterSpeed [depthSurf]Baseline 2000-2019

# check layer info
swd_id <- "swd_baseline_2000_2019_depthsurf"
info_layer(swd_id)

sws_id <- "sws_baseline_2000_2019_depthsurf"
info_layer(sws_id)

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

# recall: nw_shelf - object containing shapefile of extent
nw_swd <- raster::crop(swd_layer, nw_shelf)
nw_swd <- terra::mask(nw_swd, mask = terra::vect(nw_shelf))
mapview::mapview(nw_swd, na.color = NA)

nw_sws <- raster::crop(sws_layer, nw_shelf)
nw_sws <- terra::mask(nw_sws, mask = terra::vect(nw_shelf))
mapview::mapview(nw_sws, na.color = NA)