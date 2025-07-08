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
osc <- nc_open("./oscar_currents_final_19930108.nc")

# access "u" variable
u <- ncvar_get(osc, "u")

# rasterise
z <- rast(u)

# provide extent
ext(z) <- c(-180, 180, -90, 90)

# flip north-side up
z <- flip(z)

# projection
crs(z) <- "EPSG:4326"

# plot
plot(z)
mapview::mapview(z)

# The orientation of the map (Pacific Ocean in the centre so east and west do not line up). Check projection