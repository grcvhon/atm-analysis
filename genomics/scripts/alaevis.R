


# TESS
# STRUCTURE-like analyses
laevis_tess <- 
  algatr::tess_ktest(gen = laevis_dosage, 
                     coords = laevis_nw_proj, 
                     Kvals = 1:7, 
                     ploidy = 2, 
                     K_selection = "auto")
# ^note that K changes every iteration ...

# Get TESS object and best K from results (i.e., laevis_tess)
laevis_tessobj <- laevis_tess$tess3_obj
laevis_bestK <- laevis_tess[["K"]]

# Get matrix of ancestry coefficients
laevis_qmat <- 
  tess3r::qmatrix(tess3 = laevis_tessobj, 
                  K = laevis_bestK)

# Prepare for kriging
mspec_sh_krig <- (mspec_sh_stack[[1]]) # environmental PC1, nw shape
y_krig_sh_raster <- raster::projectRaster(mspec_sh_krig, crs = "epsg:3112")

x <- sf::st_as_sf(laevis_nw_coords, coords = c("x", "y"), crs = 4326)
x_proj <- sf::st_transform(x, crs = 3112) # EPSG:3112; GDA94 Geoscience Australia projection 

z_krig_sh_admix <- 
  algatr::tess_krig(qmat = laevis_qmat, 
                    coords = x_proj, 
                    grid = y_krig_sh_raster)

z_map_sh_admix <- 
  algatr::tess_ggplot(z_krig_sh_admix,
                      plot_method = "maxQ",
                      plot_axes = TRUE,
                      coords = x_proj)
