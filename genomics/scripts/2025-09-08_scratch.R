### script below depends on `resistance_gdistance_passage.R`

# bearing

# generate transition layer
swd_layer_tr <- transition(swd_layer_raster, transitionFunction = mean, directions = 8) %>% 
  geoCorrection(type = "c", multpl = F)

ggplot(as.data.frame(swd_layer_raster, xy=T)) + 
  geom_raster(aes(x=x, y=y, fill = bearing)) + 
  geom_point(data=as.data.frame(samp_coords), aes(x=longitude, y=latitude), size=2, col="red") +
  scale_fill_continuous(na.value=NA) + 
  theme_map() + theme(legend.position = "right")


# random points across extent
library(sf)
set.seed(123)
rand_points <- st_sample(nw_shelf, size = 50)
df_rand_points <- sfheaders::sfc_to_df(rand_points)
df_rand_points <- df_rand_points[,-c(1,2)]
colnames(df_rand_points) <- c("longitude", "latitude")
df_rand_points

# generate combinations across 100 random points
rand_pts_combn <- 
  combn(nrow(df_rand_points),2) %>%
  t() %>%
  as.matrix()


random_bearing_passages <- list()                                                     # Create a list to store the passage probability rasters in

system.time(                                                           # Keep track of how long this takes
  for (i in 1:nrow(rand_pts_combn)) {           
    locations <- SpatialPoints(rbind(df_rand_points[rand_pts_combn[i,1],1:2],     # create origin points
                                     df_rand_points[rand_pts_combn[i,2],1:2]),   # create destination (or goal) points, to traverse
                               proj4string = CRS("+init=epsg:4326"))
    random_bearing_passages[[i]] <- passage(swd_layer_tr,                                   # run the passage function 
                                     origin=locations[1],                 # set orgin point
                                     goal=locations[2],                   # set goal point
                                     theta = 0.00001)                             # set theta (tuning parameter, see notes below)
    print(paste((i/nrow(rand_pts_combn))*100, "% complete"))
  }
)

random_bearing_passages <- stack(bearing_passages)                                            # create a raster stack of all the passage probabilities
random_bearing_passages_overlay <- sum(random_bearing_passages)/nrow(rand_pts_combn)                       # calculate average

# save output as csv
library(terra)
ovr <- rast(random_bearing_passages_overlay)
df_ovr <- as.data.frame(ovr, xy = TRUE)

library(utils)
write.csv(df_ovr, "location")

colors <- c("grey50", viridis_pal(option="inferno", begin = 0.3, end = 1)(20))
ggplot(as.data.frame(random_bearing_passages_overlay, xy=T)) + 
  geom_raster(aes(x=x,y=y,fill=layer)) +
  scale_fill_gradientn(colors = colors, na.value = NA) + 
  #geom_point(data=as.data.frame(samp_coords), aes(x=longitude, y=latitude), size=1, col="red") +
  theme_map() +
  theme(legend.position = "right")


# regular
reg_rand_points <- st_sample(nw_shelf, size = 1000, type = "regular", replace = FALSE)
plot(reg_rand_points, pch = 19, cex = 0.1)
df_reg_rand_points <- sfheaders::sfc_to_df(reg_rand_points)
df_reg_rand_points <- df_reg_rand_points[,-c(1,2)]
colnames(df_reg_rand_points) <- c("longitude", "latitude")
print(dim(df_reg_rand_points)) # check number of rows and columns
print(dim(unique(df_reg_rand_points))) # check unique number of rows and columns
print(dim(unique(df_reg_rand_points)) == dim(df_reg_rand_points)) # should be TRUE TRUE

# hexagonal
hex_grid <- st_make_grid(nw_shelf[1,], what = "corners", square = F, n = 50)
plot(hex_grid, cex = 0.3)
hex_grid

df_hex_grid <- sfheaders::sfc_to_df(hex_grid)
df_hex_grid <- df_hex_grid[,-c(1,2)]
colnames(df_hex_grid) <- c("longitude", "latitude")
print(dim(df_hex_grid)) # check number of rows and columns
print(dim(unique(df_hex_grid))) # check unique number of rows and columns
print(dim(unique(df_hex_grid)) == dim(df_hex_grid)) # should be TRUE TRUE

# Comes up with duplicate sets of coordinates

# spatially balanced
# dssduoa - adopt from other projects

# required nw_shelf
plot(nw_shelf)

library(sf)
nw_shelf_utm <- st_transform(nw_shelf, crs = 32750)

library(dssduoa)
nw_shelf_region <- make.region(shape = nw_shelf_utm)

#ecoflow_points <- st_as_sfc(st_bbox(nw_shelf_utm)) #|>
  #st_make_grid(cellsize = 1000) |>
  #st_as_sf()

#ecoflow_points_inPoly <- st_intersection(ecoflow_points, nw_shelf_utm)

library(spsurvey)
set.seed(100)
ecoflow_pts <- grts(nw_shelf_utm, n_base = 100)

plot(ecoflow_pts)

ecoflow_pts_sf <- st_as_sf(ecoflow_pts$sites_base)

library(ggplot2)
ggplot() + geom_sf(data = nw_shelf_utm) + geom_sf(data = ecoflow_pts_sf)

ecoflow_pts_lon <- ecoflow_pts$sites_base$lon_WGS84
ecoflow_pts_lat <- ecoflow_pts$sites_base$lat_WGS84
ecoflow_pts_coords <- cbind(ecoflow_pts_lon,ecoflow_pts_lat)
colnames(ecoflow_pts_coords) <- c("longitude","latitude")
df_ecoflow_pts_coords <- as.data.frame(ecoflow_pts_coords)
dim(df_ecoflow_pts_coords) == dim(unique(df_ecoflow_pts_coords)) # no duplicate points

# generate origin/destination combinations across the 1000 points
ecoflow_pts_comb <- 
  combn(nrow(df_ecoflow_pts_coords),2) %>%
  t() %>%
  as.matrix()
ecoflow_pts_comb

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




