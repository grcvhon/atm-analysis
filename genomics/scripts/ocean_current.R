# Global Ocean Surface Currents - Monthly Mean
# ESRI https://www.arcgis.com/home/item.html?id=b02f417ebbed4dc69edefd848dc69715

library(raster)
library(terra)
library(mapview)

nw_current <- rast("./genomics/ocean-cur/nw_1000.tif")
View(nw_current)
mapview(nw_current, na.color = NA)