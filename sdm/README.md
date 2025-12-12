# Species distribution modelling
This directory contains code and input data for developing an updated species distribution model (SDM) for the Short-nosed (<i>Aipysurus apraefrontalis</i>) and Leaf-scaled sea snakes (<i>A. foliosquama</i>). In developing these species-specific SDMs, we will incorporate the following layers as predictor variables: (1) environmental and habitat layers, (2) [genetic layer](https://github.com/grcvhon/atm-analysis/tree/master/genetic_layer), (3) [passage probability layer](https://github.com/grcvhon/atm-analysis/tree/master/passage_layer).<br>
<br>
The code below, at the moment, displays the SDM workflow so far.

#### 1) Input model extent
```r
# Load nwshelf shapefile
nw_shelf <- st_read("./nw-shelf/NWShelf.shp", quiet = TRUE) %>% 
  st_transform(4326)
```

#### 2) Input occurrence data
The occurrence data which we will use were sourced from:<br>
1. `trawled_seasnakes.xlsx`(2009, 2014-2021)<br>
2. [Atlas of Living Australia](https://ala.org.au) (downloaded on 21 January 2025)
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
Plot of occurrence dataset. Short-nosed sea snake (yellow), Leaf-scaled sea snake (maroon).
</div>
</p>