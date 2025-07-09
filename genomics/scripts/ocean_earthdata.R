# 1) sign up for an Earthdata profile
#
# 2) set up cmdline to recognise Earthdata credentials
# > cd ~
# > touch .netrc
# > echo "machine urs.earthdata.nasa.gov login toughturf password +0u6hturF9801" > .netrc
# > chmod 0600 .netrc
#
# 3) create cookies for efficiency; makes credentials persist
# > cd ~
# > touch .urs_cookies
#
# 4) download data - modify date range to download
# Ocean Surface Current Analyses Real-time (OSCAR) Surface Currents - Final 0.25 Degree (Version 2.0)
# url: https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_FINAL_V2.0#
# > podaac-data-downloader -c OSCAR_L4_OC_FINAL_V2.0 -d ./ --start-date 1993-01-01T00:00:00Z --end-date 1993-01-08T00:00:00Z -e .nc

# load packages
library(ncdf4)
library(terra)

# open netcdf file
osc <- nc_open("./genomics/ocean-cur/podaac_data/oscar_currents_final_19930108.nc")

# access "u" (eastward) and "v" (northward) variables
u <- ncvar_get(osc, "u")
v <- ncvar_get(osc, "v")

# rasterise
z <- rast(u)
n <- rast(v)

# flip north-side up
z <- flip(z)
n <- flip(n)

# provide extent
z <- terra::rotate(z)
ext(z) <- c(-180, 180, -89.875,89.875) 
n <- terra::rotate(n)
ext(n) <- c(-180, 180, -89.875,89.875)

# res(z) and res(n) is 0.25 0.25 with these numbers; based on dataset, spatial coverage: -180, 180, -89.75, 89.75

# projection
crs(z) <- "EPSG:4326"
crs(n) <- "EPSG:4326"

# eastward
nw_z <- raster::crop(z, nw_shelf)
nw_z <- terra::mask(nw_z, mask = terra::vect(nw_shelf))
mapview::mapview(nw_z, na.color = NA)

# northward
nw_n <- raster::crop(n, nw_shelf)
nw_n <- terra::mask(nw_n, mask = terra::vect(nw_shelf))
mapview::mapview(nw_n, na.color = NA)

# formula for magnitude
# sqrt(u^2 + v^2)
nw_zn_mag <- sqrt(nw_z^2 + nw_n^2)
plot(nw_zn_mag)

# formula for bearing
# atan2(v,u)*180/pi
nw_zn_ber <- atan2(nw_z,nw_n)*180/pi
plot(nw_zn_ber)