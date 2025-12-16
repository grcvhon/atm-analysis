# Species distribution modelling
This directory contains code and input data for developing an updated species distribution model (SDM) for the Short-nosed (<i>Aipysurus apraefrontalis</i>) and Leaf-scaled sea snakes (<i>A. foliosquama</i>). In developing these species-specific SDMs, we will incorporate the following layers as predictor variables: (1) environmental and habitat layers, (2) [genetic layer](https://github.com/grcvhon/atm-analysis/tree/master/genetic_layer), (3) [passage probability layer](https://github.com/grcvhon/atm-analysis/tree/master/passage_layer).<br>
<br>
The code below, at the moment, displays the SDM workflow so far.

### 1) Input model extent
```r
# Load nwshelf shapefile
nw_shelf <- st_read("./nw-shelf/NWShelf.shp", quiet = TRUE) %>% 
  st_transform(4326)
```
### 2) Input environmental rasters
```r
# Load single environmental raster
bathymetry <- raster("./predictor-variables/bathymetry.asc")

# Load initial set of environmental rasters
env_init <- stack("./predictor-variables/sal_mean.asc",
                  "./predictor-variables/sal_amp.asc",
                  "./predictor-variables/bathymetry.asc",
                  "./predictor-variables/sst_mean.asc",
                  "./predictor-variables/sst_amp.asc",
                  "./predictor-variables/chlor_mean.asc",
                  "./predictor-variables/DistToLand.asc",
                  "./predictor-variables/DistToReef.asc",
                  "./predictor-variables/DistToFW.asc")
```

### 3) Input occurrence data
The occurrence data were compiled into a master dataset ([`ATM_master-occurrence-dataset.csv`](https://github.com/grcvhon/atm-analysis/tree/master/sdm/occurrence-data/ATM_master-occurrence-dataset.csv)). These records of Short-nosed and Leaf-scaled sea snakes were sourced from the following:<br>
1. `trawled_seasnakes.xlsx`(2009, 2014 - 2021)<br>
2. Atlas of Living Australia (downloaded on 21 January 2025)
3. DPIRD Fisheries surveys in Exmouth Gulf and Shark Bay, Western Australia (June 2024 - December 2025)
4. One spreadsheet (KLS Lab catalogue records)

Occurrence data are point coordinates. When captured via trawling (with start and end coordinates), the mid latitude and mid longitude values were used as proxy for point coordinates; and were calculated as the average latitude and longitude values of the start and end coordinates.
```r
occurrence <- read_csv("./occurrence-data/ATM_master-occurrence-dataset.csv")

# [function] filter by species
filter_by_species <- function(species) {
  occurrence %>%
    mutate(
      lat = as.numeric(MidLat),
      lat = if_else(is.na(lat), as.numeric(Latitude), lat),
      long = as.numeric(MidLong),
      long = if_else(is.na(long), as.numeric(Longitude), long)
      ) %>%
    filter(Species == species,
           !is.na(lat),
           !is.na(long)) %>%
    clean_names()
}

# filtered df by species
short_occ <- filter_by_species("apraefrontalis")
leaf_occ <- filter_by_species("foliosquama")

# generate sf from df
short_sf <- 
  short_occ %>% 
  st_as_sf(coords = c("long", "lat"), 
           crs = 4326) %>% 
  distinct(.keep_all = T)

leaf_sf <- 
  leaf_occ %>% 
  st_as_sf(coords = c("long", "lat"), 
           crs = 4326) %>% 
  distinct(.keep_all = T)

# map distribution and colour by species
ggplot() +
  geom_sf(data = nw_shelf, fill = NA) +
  geom_sf(data = short_sf, col = "orange", cex = 3, pch = 16) +
  geom_sf(data = leaf_sf, col = "maroon", cex = 3, pch = 16) + 
  annotation_scale(mapping = aes(location = "br")) +
  theme_bw()
```
<p align = center>
<img src="https://raw.githubusercontent.com/grcvhon/atm-analysis/master/sdm/plots/plot_occurrence-data.png", width = 75%, height = 75%>
<div align = "center">
Occurrence dataset plotted within northwest shelf extent. Short-nosed sea snake (yellow), Leaf-scaled sea snake (maroon).
</div>
</p>

### 4) Bias layer
We will generate a bias layer using point occurrences from `2020-05-20_SnakeOcc.csv`
```r
# read bias point occurrences
bias_pts <- read_csv("./bias-layer/2020-05-20_SnakeOcc.csv")

# generate sf object from bias_pts
bias_sf <- bias_pts %>% 
  st_as_sf(coords = c("longitude", "latitude"),
           crs = 4326) %>% 
  distinct(.keep_all = T)

# visualise
ggplot() +
  geom_sf(data = bias_sf)

# vectorise bias_sf
bias_vect <- vect(bias_sf)

# mask bias points with nwshelf boundary
nw_bias_sf <- mask(bias_vect, mask = vect(nw_shelf)) %>% 
  st_as_sf()

# visualise
ggplot() +
  geom_sf(data = nw_shelf, fill = NA) +
  geom_sf(data = nw_bias_sf, aes(col = species))

# generate point pattern object
bias_ppp <- nw_bias_sf %>% 
  st_transform(crs = 3577) %>% 
  st_as_sf() %>% 
  as.ppp()

# calculate probability density of points
bias_prob <- bias_ppp %>% 
  density(., sigma = 0.05) %>% 
  rast()

plot(bias_prob)
crs(bias_prob) <- "EPSG:3577"

bias_prob <- bias_prob %>% 
  project(., y= "EPSG:4326") %>% 
  mask(mask = vect(nw_shelf)) %>% 
  resample(x = ., y = rast(bathymetry)) %>% 
  terra::scale()
  
bias_prob[values(bias_prob) < 0] <- NA

nx_bias <- minmax(bias_prob)
bias_prob <- (bias_prob - nx_bias[1,]) / (nx_bias[2,] - nx_bias[1,])

# visualise bias layer
ggplot() +
  geom_spatraster(data = bias_prob) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = nw_shelf, fill = NA) + 
  annotation_scale(mapping = aes(location = "br")) +
  theme_bw()
```
<p align = center>
<img src="https://raw.githubusercontent.com/grcvhon/atm-analysis/master/sdm/plots/plot_bias-probability.png", width = 75%, height = 75%>
<div align = "center">
Bias layer with probability density (sigma = 0.05).
</div>
</p>

### 5) Pseudoabsence layer