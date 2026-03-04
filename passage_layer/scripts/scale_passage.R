# scale passage probability values

# set working directory
setwd("C:/Users/a1235304/Dropbox/Short-nosed and Leaf-scaled sea snake TSSC/atm-analysis/passage_layer/")

# set palette
col_palette <- heat.colors(100)
palette(col_palette)

# scaling function
raster01 <- function(r){
  
  # get the min max values
  minmax_r = range(values(r), na.rm=TRUE) 
  
  # rescale 
  return( (r-minmax_r[1]) / (diff(minmax_r)))
}

# import passage layer asc file
passage_asc <- terra::rast(x = "output/sbs_bearing_seed100_100pts_03h49m31s/passage_layer.asc")

# scaled as passage probability
scaled_passage <- raster01(passage_asc)
plot(scaled_passage)

# save as ascii raster
terra::writeRaster(
  x = scaled_passage,
  filename = "./output/scaled_passage.asc")

# scaled as non-passage probability
scaled_nonpass <- raster01(1-passage_asc) # 1 - passage probability values
plot(scaled_nonpass)

# save as ascii raster
terra::writeRaster(
  x = scaled_nonpass,
  filename = "./output/scaled_nonpass.asc")

