# Global Ocean Surface Currents - Monthly Mean
# ESRI https://www.arcgis.com/home/item.html?id=b02f417ebbed4dc69edefd848dc69715

# *** insert how tif was exported from ArcGIS ***
# *** insert how monthly means were subset from multidimensional data ***

library(raster)
library(terra)
library(mapview)

# *** have to find out which time slice ***
nw_current <- rast("./genomics/ocean-cur/2001-03-01.tif") # tif manually subset in ArcGIS

# check layers
plot(nw_current$Band_1) # assuming Band_1 is U (eastward velocity)
plot(nw_current$Band_2) # assuming Band_2 is V (northward velocity)

# assign U and V
u <- nw_current$Band_1
v <- nw_current$Band_2

# formula for magnitude
# sqrt(u^2 + v^2)
nw_current_mag <- sqrt(u^2 + v^2)
plot(nw_current_mag)

# formula for bearing
# atan2(v,u)*180/pi
nw_current_ber <- atan2(v,u)*180/pi
plot(nw_current_ber)

# one page plot
par(mfrow=c(2,1))
plot(nw_current_mag)
title("Magnitude", line = 3)
plot(nw_current_ber)
title("Bearing", line = 3)

# only keep within nw_shelf shape extent
nw_current_mag <- terra::project(x = nw_current_mag, y = "EPSG:4326")
plot(nw_current_mag)
nw_current_mag_ext <- raster::crop(nw_current_mag, nw_shelf)
nw_current_mag_ext <- terra::mask(nw_current_mag_ext, mask = terra::vect(nw_shelf))
mapview(nw_current_mag_ext, na.color = NA) + mapview(nw_shelf, alpha.region = 0)

nw_current_ber <- terra::project(x = nw_current_ber, y = "EPSG:4326")
plot(nw_current_ber)
nw_current_ber_ext <- raster::crop(nw_current_ber, nw_shelf)
nw_current_ber_ext <- terra::mask(nw_current_ber_ext, mask = terra::vect(nw_shelf))
mapview(nw_current_ber_ext, na.color = NA) + mapview(nw_shelf, alpha.region = 0)
