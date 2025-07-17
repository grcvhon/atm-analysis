### try BioOracle instead...
library(biooracler)

# BioOracle layers are at the spatial resolution of 0.05 x 0.05 decimal degrees and decadal temporal resolution

# list available layers
layers <- list_layers()
swd <- list_layers("SeaWaterDirection")
#View(swd)
# want swd_baseline_2000_2019_depthsurf | Bio-Oracle SeaWaterDirection [depthSurf]Baseline 2000-2019.

sws <- list_layers("SeaWaterSpeed")
#View(sws)
# want sws_baseline_2000_2019_depthsurf | Bio-Oracle SeaWaterSpeed [depthSurf]Baseline 2000-2019

# check layer info
swd_id <- "swd_baseline_2000_2019_depthsurf"
#info_layer(swd_id)

sws_id <- "sws_baseline_2000_2019_depthsurf"
#info_layer(sws_id)

# set constraints before downloading data

# time based on info_layer
oracle_time <- c("2000-01-01T00:00:00Z", "2010-01-01T00:00:00Z")
# extent based on "nw_shelf" shape i.e., ext(nw_shelf)
oracle_lat <- c(-26.7363816780876, -9.69138797000501)
oracle_lon <- c(111.544995584098, 130.342214992267)

constraints <- list(oracle_time, oracle_lat, oracle_lon)
names(constraints) <- c("time", "latitude", "longitude")

# define wanted variables to download
swd_var <- "swd_mean" # sea water direction - mean
sws_var <- "sws_mean" # sea water speed - mean

# perform download
swd_layer <- download_layers(swd_id, swd_var, constraints)
sws_layer <- download_layers(sws_id, sws_var, constraints)

### side quest: marmap adaptation
  #swd_rast <- raster::raster(swd_layer)
  #plot(swd_rast)
  #swd_bathy <- marmap::as.bathy(swd_rast)
  #plot(swd_bathy, image = TRUE)
  
  #library(raster)
  #library(terra)
  #library(marmap)
  swd_mask <- raster::crop(swd_layer, nw_shelf)
  swd_mask <- terra::mask(swd_mask, mask = terra::vect(nw_shelf))
  swd_mask <- raster::raster(swd_mask)
  plot(swd_mask)
  swd_mask_bathy <- marmap::as.bathy(swd_mask)
  summary(swd_mask_bathy)
  swd_trans <- 
  swd_lcdist <- marmap::lc.dist(swd_trans, sub_laevis, res = "path")
  lapply(nw_lcdist200,lines,col="dodgerblue",lwd=3,lty=1)
  
  #devtools::install_github("jabiologo/rWind")
  library(rWind)
  # in the following lines, I'm trying to copy the format of `data(wind.data)` from the `rWind` package
  
  # working with swd_layer downloaded from biooracle
  swd_layer_coords <- terra::crds(swd_layer$swd_mean_1,  df = TRUE) # get coords from swd_mean_1 layer and into a df 
  swd_layer_dr1_df <- as.data.frame(swd_layer$swd_mean_1) # get mean direction  1
  sws_layer_sp1_df <- as.data.frame(sws_layer$sws_mean_1) # get mean speed 1
  ocean_c <- cbind(swd_layer_coords,swd_layer_dr1_df,sws_layer_sp1_df) # combine columns
  colnames(ocean_c) <- c("lon","lat","dir","speed") # rename headers
  class(ocean_c) <- c("rWind", "data.frame") # reclass
  # plot to see if map looks as expected
  library(ggplot2)
  ggplot(ocean_c, aes(lon, lat, color = dir)) + 
    geom_point() + 
    scale_color_continuous(type = "viridis")
  
  # following rWind workflow (example from vignette below)
  ocean_c_ras <- rWind::wind2raster(ocean_c)
  ocean_conduct <- flow.dispersion(ocean_c_ras)
  
  # provide origin and goal coordinates for path generation
  sub_coords <- data.frame(longitude = c(114.2825, 114.0446),
                           latitude = c(-22.11952, -26.39628))
  sub_coords # first row ExmouthG, second row SharkB
  
  GtoS <- shortestPath(ocean_conduct, 
                       goal = c(sub_coords[2,1], sub_coords[2,2]),
                       origin = c(sub_coords[1,1], sub_coords[1,2]),
                       output = "SpatialLines")
  
  # note order of origin and goal args swapped
  StoG <- shortestPath(ocean_conduct, 
                       origin = c(sub_coords[2,1], sub_coords[2,2]),
                       goal = c(sub_coords[1,1], sub_coords[1,2]),
                       output = "SpatialLines")
  
  plot(swd_layer$swd_mean_1)
  lines(GtoS, col = "orange", lwd = 3) 
  lines(StoG, col = "orangered", lwd = 3)
  #
  # issue/challenge: line crosses the land...
  #
  # if fixed: can get coords of the dots that build the line
  #   `GtoS@lines[[1]]@Lines[[1]]@coords`
  #   `StoG@lines[[1]]@Lines[[1]]@coords`
  # 
  # can compute distance by connecting the dots. Not sure if readily doable in the package or R
  # this computed distance represents the distance that follows the least cost based on current direction
  # which I think is analogous to the resistance distance which can be regressed against genetic distance
  # to answer question, does current direction impose resistance to population genetic connectivity.
  # 
  # we may not have good maps to visualise this though, although circuitscape has the capacity to do this but
  # we run into the challenge of classifying current bearing values with their computed resistance value
  # based on their straightline bearing, which is not valid as some straightlines pass through land features.
  
  # The lines represent the least cost path accounting for current direction (orange is Exmouth to Shark Bay and then SB to EG in red)
  #
  # the issue/challenge here is it is making path into the land, although I think there shouldn't be values in the whitespace
  #
  # if we can get around this issue, I think we can do the following:
  # get coordinates of the dots that form the line 
  # based on those coordinates, compute the total distance
  # 
  # This computed distance would represent how far an individual from SB to EG would travel if it were to follow the path of least resistance due to current direction.
  # 
  # Then I think we can use this computed distance and regress it with genetic distance to answer the question: does current direction (here represented by the computed distance) impose resistance to the population genetic connectivity?
  # 
  # another challenge though is that we won't have good maps like those produced in circuitscape, since this approach only finds one solution (i.e., one line) therefore no exploration of other potential pathways.
  
  # some calculations in the function `flow.dispersion`
  
  # NOTES/COMMENTS
  # 
  # The lines generated so far the least cost path accounting for current direction (orange is Exmouth to Shark Bay and then SB to EG in red).
  # Although, the obvious issue/challenge is that paths were generated over land, although I think there shouldn't be values in the whitespace.
  #
  # If we can get around this issue, I think we can do the following:
  # - get coordinates of the dots that form the line 
  # - then based on those coordinates, compute the total distance (i.e., connecting the coordinate points to form a line and measure that line)
  # 
  # This computed distance would represent how far an individual would travel from SB to EG if it were to follow the path that has 
  # the least resistance due to current direction. I think, then, we can use this computed distance and regress it with genetic distance 
  # to answer the question: does current direction (here represented by the information contained in the computed distance) impose resistance 
  # that can influence population genetic connectivity?
  # 
  # another challenge though is that we won't have maps like those produced in circuitscape since the rWind approach 
  # only finds one solution and therefore no exploration of other potential pathways (while not best are still possible pathways) across the extent.
  #
  # Maybe useful code: the rWind package provided the calculations for cost in the function `flow.dispersion` and the cost function itself `cost.FMGS`
  
  #example
  data(wind.data)
  wd_ras <- wind2raster(wind.data)
  conductance <- flow.dispersion(wd_ras)
  AtoB <- shortestPath(conductance, c(-5.5,37), c(-5.5,35), output = "SpatialLines")
  BtoA <- shortestPath(conductance, c(-5.5,35), c(-5.5,37), output = "SpatialLines")
  
# recall: nw_shelf - object containing shapefile of extent
nw_swd <- raster::crop(swd_layer, nw_shelf)
nw_swd <- terra::mask(nw_swd, mask = terra::vect(nw_shelf))
mapview::mapview(nw_swd, na.color = NA)

nw_sws <- raster::crop(sws_layer, nw_shelf)
nw_sws <- terra::mask(nw_sws, mask = terra::vect(nw_shelf))
mapview::mapview(nw_sws, na.color = NA)


###

library(marmap)
library(raster)
library(gdistance)

# get bathymetric data for Australia
nw_bathy <- getNOAA.bathy(lon1=110.544995584098, lon2=118.342214992267, 
                           lat1=-18.69138797000501, lat2=-29.7363816780876,
                           resolution=1, keep=TRUE)

# subset laevis sites
sub_laevis <- laevis_nw[c(122,170),c(4,5)] # Exmouth and Shark Bay
sub_laevis

nw_trans <- trans.mat(nw_bathy)
nw_trans200 <- trans.mat(nw_bathy, min.depth = 0, max.depth = -200)

nw_lcdist <- lc.dist(nw_trans,sub_laevis,res="path")
nw_lcdist200 <- lc.dist(nw_trans200,sub_laevis,res="path")

plot(nw_bathy, image=TRUE)
lapply(nw_lcdist200,lines,col="orange",lwd=3,lty=1)
points(sub_laevis, pch = 19, col = "red", cex = 1.5)

###


# Cost computation following Muñoz et al., 2004; Felicísimo et al., 2008
#' @rdname flow.dispersion
#' @export
cost.FMGS <- function(wind.direction, wind.speed, target, type = "passive") {
  dif <- (abs(wind.direction - target))
  # If dif > 180 and is not NA
  dif[dif > 180 & !is.na(dif)] <- 360 - dif[dif > 180 & !is.na(dif)] # Modified from the original function
  if (type == "passive") {
    # In "passive" type, if dif > 90, movement is not allowed
    dif[dif >= 90 & !is.na(dif)] <- Inf # check
    # For sea currents, dif could be NA is there are lands around
    dif[is.na(dif)] <- Inf # Modified from the original function
    # Here we apply the formula in Felicísimo et al. 2008
    dif[dif < 90] <- 2 * dif[dif < 90]
    dif[dif == 0] <- 0.1
  }
  else {
    print("Only passive movements are currently allowed")
    # For sea currents, dif could be NA is there are lands around
    #dif[is.na(dif)] <- Inf
    # For "active" type movements against flow are allowed, so simply we
    # multiply by 2 the dif, following Felicísimo et al., 2008
    #dif[!is.na(dif)] <- 2 * dif[!is.na(dif)]
    #dif[dif == 0] <- 0.1
  }
  
  wind.speed[is.na(wind.speed)] <- 0
  
  dif / wind.speed
}

#' Compute flow-based cost or conductance
#'
#' flow.dispersion_int computes movement conductance through a flow either, sea
#' or wind currents. It implements the formula described in Felícisimo et al.
#' 2008:
#'
#' Cost=(1/Speed)*(HorizontalFactor)
#'
#' being HorizontalFactor a "function that incrementally penalized angular
#' deviations from the wind direction" (Felicísimo et al. 2008). Only _passive_
#' movements are currently allowed.
#'
#' @param stack RasterStack object with layers obtained from wind2raster
#' function ("rWind" package) with direction and speed flow values.
#' @param fun A function to compute the cost to move between cells. The default
#' is \code{cost.FMGS} from Felicísimo et al. (2008), see details.
#' @param output This argument allows to select different kinds of output. "raw"
#' mode creates a matrix (class "dgCMatrix") with transition costs between all
#' cells in the raster. "transitionLayer" creates a TransitionLayer object with
#' conductance values to be used with "gdistance" package.
#' @param ... Further arguments passed to or from other methods.
#' @return In "transitionLayer" output, the function returns conductance values
#' (1/cost)to move between all cells in a raster having into account flow speed
#' and direction obtained from wind.fit function("rWind" package). As wind or
#' sea currents implies directionality, flow.dispersion produces an anisotropic
#' conductance matrix (asymmetric). Conductance values are used later to built a
#' TransitionLayer object from "gdistance" package.
#'
#' In "raw" output, flow.dispersion creates a sparse Matrix with cost values.
#' @note Note that for large data sets, it could take a while. For large study
#' areas is strongly advised perform the analysis in a remote computer or a
#' cluster.
#' @author Javier Fernández-López; Klaus Schliep; Yurena Arjona
#' @seealso \code{\link{wind.dl}}, \code{\link{wind2raster}}
#' @references
#'
#' Felicísimo, Á. M., Muñoz, J., & González-Solis, J. (2008). Ocean surface
#' winds drive dynamics of transoceanic aerial movements. PLoS One, 3(8),
#' e2928.
#'
#' Jacob van Etten (2017). R Package gdistance: Distances and Routes on
#' Geographical Grids. Journal of Statistical Software, 76(13), 1-21.
#' doi:10.18637/jss.v076.i13
#' @keywords ~anisotropy ~conductance
#' @examples
#'
#' data(wind.data)
#' wind <- wind2raster(wind.data)
#' Conductance <- flow.dispersion(wind, type = "passive")
#' \dontrun{
#' require(gdistance)
#' transitionMatrix(Conductance)
#' image(transitionMatrix(Conductance))
#' }
#' @importClassesFrom raster RasterLayer
#' @importFrom raster ncell
#' @importMethodsFrom raster as.matrix
#' @importFrom Matrix sparseMatrix
#' @importFrom gdistance transition transitionMatrix<-
#' @keywords internal
flow.dispersion_int <- function(stack, fun = cost.FMGS, output = "transitionLayer",
                                ...) {
  output <- match.arg(output, c("raw", "transitionLayer"))
  
  DL <- as.matrix(stack$direction)
  SL <- as.matrix(stack$speed)
  M <- matrix(as.integer(1:ncell(stack$direction)),
              nrow = nrow(stack$direction), byrow = TRUE
  )
  nr <- nrow(M)
  nc <- ncol(M)
  
  ###################################################################
  
  directions <- c(315, 0, 45, 270, 90, 225, 180, 135)
  
  ###################################################################
  
  # Go Nortwest
  
  north.west.from <- as.vector(M[-1, -1])
  north.west.to <- as.vector(M[-nr, -nc])
  north.west.cost <- fun(DL[-1, -1], SL[-1, -1], directions[1], ...)
  
  ###################################################################
  
  # Go North
  
  north.from <- as.vector(M[-1, ])
  north.to <- as.vector(M[-nr, ])
  north.cost <- as.vector(fun(DL[-1, ], SL[-1, ], directions[2], ...))
  
  ###################################################################
  
  # Go Norteast
  
  north.east.from <- as.vector(M[-1, -nc])
  north.east.to <- as.vector(M[-nr, -1])
  north.east.cost <- as.vector(fun(DL[-1, -nc], SL[-1, -nc], directions[3], ...))
  
  ###################################################################
  
  # Go West
  
  west.from <- as.vector(M[, -1])
  west.to <- as.vector(M[, -nc])
  west.cost <- as.vector(fun(DL[, -1], SL[, -1], directions[4], ...))
  
  ###################################################################
  
  # Go East
  
  east.from <- as.vector(M[, -nc])
  east.to <- as.vector(M[, -1])
  east.cost <- as.vector(fun(DL[, -nc], SL[, -nc], directions[5], ...))
  
  ###################################################################
  
  # Go Southwest
  
  south.west.from <- as.vector(M[-nr, -1])
  south.west.to <- as.vector(M[-1, -nc])
  south.west.cost <- as.vector(fun(DL[-nr, -1], SL[-nr, -1], directions[6], ...))
  
  ###################################################################
  
  # Go South
  
  south.from <- as.vector(M[-nr, ])
  south.to <- as.vector(M[-1, ])
  south.cost <- as.vector(fun(DL[-nr, ], SL[-nr, ], directions[7], ...))
  
  ###################################################################
  
  # Go Southeast
  
  south.east.from <- as.vector(M[-nr, -nc])
  south.east.to <- as.vector(M[-1, -1])
  south.east.cost <- as.vector(fun(DL[-nr, -nc], SL[-nr, -nc], directions[8], ...))
  
  ###################################################################
  
  ii <- c(north.west.from, north.from, north.east.from, west.from, east.from, south.west.from, south.from, south.east.from)
  jj <- c(north.west.to, north.to, north.east.to, west.to, east.to, south.west.to, south.to, south.east.to)
  xx <- c(north.west.cost, north.cost, north.east.cost, west.cost, east.cost, south.west.cost, south.cost, south.east.cost)
  
  tl <- sparseMatrix(i = ii, j = jj, x = xx)
  if (output == "raw") {
    return(tl)
  }
  if (output == "transitionLayer") {
    tmp <- transition(stack$direction, transitionFunction = function(x) 0, directions = 8)
    transitionMatrix(tmp) <- sparseMatrix(i = ii, j = jj, x = 1 / xx)
    return(tmp)
  }
  return(NULL)
}
  



