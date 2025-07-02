# Global Ocean Surface Currents - Monthly Mean
# ESRI https://www.arcgis.com/home/item.html?id=b02f417ebbed4dc69edefd848dc69715

# *** insert how tif was exported from ArcGIS ***

library(raster)
library(terra)
library(mapview)

# *** have to find out which time slice ***
nw_current <- rast("./genomics/ocean-cur/nw_1000.tif")

# check layers
plot(nw_current$Band_1) # assuming Band_1 is U (eastward velocity)
plot(nw_current$Band_2) # assuming Band_2 is V (northward velocity)

# assign U and V
u <- nw_current$Band_1
v <- nw_current$Band_2

# formula for magnitude
# sqrt(u^2 + v^2)
nw_current_dir <- sqrt(u^2 + v^2)
plot(nw_current_dir)

# formula for bearing
# atan2(v,u)*180/pi
nw_current_ber <- atan2(v,u)*180/pi
plot(nw_current_ber)

# plot side-by-side
par(mfrow=c(2,1))
plot(nw_current_dir)
title("Magnitude", line = 3)
plot(nw_current_ber)
title("Bearing", line = 3)
