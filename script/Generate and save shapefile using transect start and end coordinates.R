## Generate and save shapefile using ##
## start and end coordinates ##
### of transects/trawls ###

transect_to_shp <- function(species) {
  sp <- species %>% 
    # select the columns of start and end lats and longs from data table
    select(c(start_long,start_lat,end_long,end_lat)) %>% 
    # make sure that only non-NAs are included
    filter(apply(species[, c("start_long", "start_lat", "end_long", "end_lat")], 1, function(x) all(!is.na(x))))
  
  # create an empty list to store transect line geometries
  trn_lines <- list()
  
  for(i in 1:nrow(sp)) {
    start_gps <- c(sp$start_long[i], sp$start_lat[i])
    end_gps <- c(sp$end_long[i], sp$end_lat[i])
  
    trn_gps <- rbind(start_gps, end_gps)
    trn_line <- st_sfc(st_linestring(trn_gps), crs = 4326)
    trn_lines[[i]] <- trn_line
  }
  
  trn_sf <- st_sf(geometry = do.call(c, trn_lines))
}

aprae_trn <- transect_to_shp(aprae)
mapview(aprae_trn)
st_write(aprae_trn, "./data/shapefiles/apraefrontalis_transect.shp")

folio_trn <- transect_to_shp(folio)
mapview(folio_trn)
st_write(folio_trn, "./data/shapefiles/foliosquama_transect.shp")