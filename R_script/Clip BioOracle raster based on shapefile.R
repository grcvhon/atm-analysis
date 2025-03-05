### Clip Bio Oracle raster based on a shapefile ###

eez <- st_read("misc/from-LM/australia-eez/Australia_EEZ.shp")
eez_crop <- st_crop(eez, 
                    xmin = 106, xmax = 129, 
                    ymin = -38, ymax = -10)
mapview(eez_crop)

bathymetry <- 
  raster("data/predictor-variables/biooracle_bathy/bathymetry_mean_Layer.tif")
mapview(bathymetry)

class(eez_crop)
class(bathymetry)

eez_crop <- 
  # project eez_crop with same projection as raster i.e., bathymetry
  st_transform(eez_crop, crs(bathymetry)) 

bathy_crop <- 
  # crops only based on the x-min/max and y-min/max
  crop(bathymetry, extent(eez_crop)) 

bathy_mask <- 
  # only area within the polygon is kept
  mask(bathy_crop, eez_crop) 

# compare
mapview(bathy_crop)
mapview(bathy_mask)