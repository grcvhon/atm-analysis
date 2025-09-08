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
reg_rand_points

# hexagonal
hex_grid <- st_make_grid(nw_shelf[1,], what = "corners", square = F, n = 1000)
plot(hex_grid, cex = 0.3)
hex_grid



