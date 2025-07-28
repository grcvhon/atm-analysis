# load package
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

#devtools::install_github("jabiologo/rWind")
library(rWind)
#install.packages("gdistance")
library(gdistance)

#### FOR REFERENCE - RUN IF NEEDED ####
# this is the example from rWind vignette (https://github.com/jabiologo/rWind/tree/master)
data(wind.data)
wd_ras <- wind2raster(wind.data)
conductance <- flow.dispersion(wd_ras)
AtoB <- shortestPath(conductance, c(-5.5,37), c(-5.5,35), output = "SpatialLines")
BtoA <- shortestPath(conductance, c(-5.5,35), c(-5.5,37), output = "SpatialLines")
#######################################

# in the following lines, I'm trying to copy the format of `data(wind.data)` from the `rWind` package

# working with swd_layer downloaded from biooracle
swd_layer_coords <- terra::crds(swd_layer$swd_mean_1,  df = TRUE) # get coords from swd_mean_1 layer and into a df 
swd_layer_dr1_df <- as.data.frame(swd_layer$swd_mean_1) # get mean direction  1
sws_layer_sp1_df <- as.data.frame(sws_layer$sws_mean_1) # get mean speed 1
ocean_c <- cbind(swd_layer_coords,swd_layer_dr1_df,sws_layer_sp1_df) # combine columns
colnames(ocean_c) <- c("lon","lat","dir","speed") # rename headers
class(ocean_c) <- c("rWind", "data.frame") # reclass

# plot to see if map looks as expected
library(ggplot2)
ggplot(ocean_c, aes(lon, lat, color = dir)) + 
  geom_point() + 
  scale_color_continuous(type = "viridis")

# now on to generating paths with rWind and gdistance (`shortestPath`)

# attempt to follow sequence of functions to generate shortest paths accounting for current bearing/direction

# using ocean_c which should be similar structure as the example `wind.data`
ocean_c_ras <- rWind::wind2raster(ocean_c) # note warning here...
ocean_conduct <- flow.dispersion(ocean_c_ras)

# provide origin and goal coordinates for path generation
sub_coords <- data.frame(longitude = c(114.2825, 114.0446),
                         latitude = c(-22.11952, -26.39628))
sub_coords # first row ExmouthG, second row SharkB

GtoS <- shortestPath(ocean_conduct, 
                     goal = c(sub_coords[2,1], sub_coords[2,2]),
                     origin = c(sub_coords[1,1], sub_coords[1,2]),
                     output = "SpatialLines")

# note order of origin and goal args swapped
StoG <- shortestPath(ocean_conduct, 
                     origin = c(sub_coords[2,1], sub_coords[2,2]),
                     goal = c(sub_coords[1,1], sub_coords[1,2]),
                     output = "SpatialLines")

plot(swd_layer$swd_mean_1)
lines(GtoS, col = "orange", lwd = 3) 
lines(StoG, col = "orangered", lwd = 3)

# NOTES/COMMENTS
# 
# The lines generated so far the least cost path accounting for current direction (orange is Exmouth to Shark Bay and then SB to EG in red).
# Although, the obvious issue/challenge is that paths were generated over land, although I think there shouldn't be values in the whitespace.
#
# If we can get around this issue, I think we can do the following:
# - get coordinates of the dots that form the line 
# - then based on those coordinates, compute the total distance (i.e., connecting the coordinate points to form a line and measure that line)
# 
# This computed distance would represent how far an individual would travel from SB to EG if it were to follow the path that has 
# the least resistance due to current direction. I think, then, we can use this computed distance and regress it with genetic distance 
# to answer the question: does current direction (here represented by the information contained in the computed distance) impose resistance 
# that can influence population genetic connectivity?
# 
# another challenge though is that we won't have maps like those produced in circuitscape since the rWind approach only finds one solution 
# and therefore does not explore/display other potential pathways (i.e., not best are still possible pathways) across the extent.
#
# Maybe useful code: the rWind package provided the calculations for cost in the function `flow.dispersion` and the cost function itself `cost.FMGS`