# gdistance vignette notes

# conductance = 1/resistance
# In gdistance, conductance rather than resistance values are
# expected in the transition matrix.

display <- function(raster){
  # Get cell coordinates
  coords <- xyFromCell(raster, cell = 1:ncell(r))
  # Get cell values
  vals <- getValues(raster)
  # Add text labels to the plot
  text(coords, labels = vals, cex = 1) # Adjust cex for text size  
}

# theoretical ----
library(gdistance) # load package

set.seed(123)


# Raster* class
r <- raster(ncol = 3, nrow = 3) # create simple 3x3 raster
r[] <- 1:ncell(r) # assign numbers 1 to 9
r
plot(r)
display(r)

# Transition* class
# Row 1 and column 1 in the transition matrix correspond to cell 1 in the original raster, 
# row 2 and column 2 to cell 2, and so on.

r[] <- 1 # use the raster created and set all its values to unit
plot(r)
display(r)
tr1 <- transition(r, transitionFunction = mean, directions = 8)
tr1 # dsCMatrix = symmetric transition
plot(raster(tr1))
display(raster(tr1))

# create an asymmetric matrix, first create a non-commutative distance function `ncdf`
# then use `ncdf` as argument in the transition function.
# set symm argument in transition as FALSE as well
r[] <- runif(9) # assign numbers between 0 and 1 (default args) to each cell (random every time)
plot(r)
display(r)

ncf <- function(x) max(x) - x[1] + x[2]
tr2 <- transition(r, transitionFunction = ncf, directions = 4, symm = FALSE)
tr2 # dgCMatrix = asymmetric transition
tr2[1:9,1:9]

plot(raster(tr2))
display(raster(tr2))

tr3 <- tr1 * tr2
tr3 <- tr1 + tr2
tr3 <- tr1 * 3
tr3 <- sqrt(tr1)

tr3[cbind(1:9, 1:9)] <- tr2[cbind(1:9, 1:9)]
tr3[1:9, 1:9] <- tr2[1:9, 1:9]
tr3[1:5, 1:5]

tr1C <- geoCorrection(tr1, type = "c")
tr2C <- geoCorrection(tr2, type = "c")

r3 <- raster(ncol = 18, nrow = 9)
r3 <- setValues(r3, runif(18 * 9) + 5)

tr3 <- transition(r3, mean, 4)
tr3C <- geoCorrection(tr3, type = "c", multpl = FALSE, scl = TRUE)
tr3R <- geoCorrection(tr3, type = "r", multpl = FALSE, scl = TRUE)

CorrMatrix <- geoCorrection(tr3, type = "r", multpl = TRUE, scl = TRUE)
tr3R <- tr3 * CorrMatrix

sP <- cbind(c(-100, -100, 100), c(50, -50, 50))
plot(sP)

costDistance(tr3C, sP)
commuteDistance(tr3R, sP)
rSPDistance(tr3R, sP, sP, theta = 1e-12, totalNet = "total")

origin <- SpatialPoints(cbind(0, 0))
rSPraster <- passage(tr3C, origin, sP[1, ], theta = 3)
plot(rSPraster)


# Hiking example ----

# read in altitude data
r <- raster(system.file("external/maungawhau.grd", package = "gdistance"))
plot(r)
r

altDiff <- function(x) {x[2] - x[1]}

hd <- transition(r, altDiff, 8, symm = FALSE)
slope <- geoCorrection(hd)

adj <- adjacent(r, cells = 1:ncell(r), pairs = TRUE, directions = 8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))

conductance <- geoCorrection(speed)
plot(raster(conductance))

A <- c(2667670, 6479000)
B <- c(2667800, 6479400)

AtoB <- shortestPath(conductance, A, B, output = "SpatialLines")
BtoA <- shortestPath(conductance, B, A, output = "SpatialLines")
plot(r, xlab = "x coordinate (m)", ylab = "y coordinate (m)", legend.lab = "Altitude (masl)")
lines(AtoB, col = "red", lwd = 2)
lines(BtoA, col = "blue")
text(A[1] - 10, A[2] - 10, "A")
text(B[1] + 10, B[2] + 10, "B")

plot(passage)

passage <- passage(conductance,A,B,theta = 0.01)


library(tidyverse)
library(viridis)
library(ggthemes)
colors <- c("grey50", viridis_pal(option="inferno", begin = 0.3, end = 1)(20))
ggplot(as.data.frame(passage, xy=T)) + 
  geom_raster(aes(x=x,y=y,fill=layer)) +
  scale_fill_gradientn(colors = colors, na.value = NA) + 
  #geom_point(data=as.data.frame(samp_coords), aes(x=longitude, y=latitude), size=1, col="red") +
  theme_map() +
  theme(legend.position = "right")
