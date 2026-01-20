# Generate conductance surface using passage probability 
# based on ocean surface current bearing

# *** For use in Circuitscape ***

# set working directory
setwd("C:/Users/a1235304/Dropbox/Short-nosed and Leaf-scaled sea snake TSSC/atm-analysis/passage_layer/")

###############################
# ***    Passage layer      ***
###############################

# import passage layer asc file
passage_asc <- terra::rast(x = "output/sbs_bearing_seed100_100pts_03h49m31s/passage_layer.asc")
plot(passage_asc)

values(passage_asc)[is.nan(values(passage_asc))] <- NA

# write new asc file
terra::writeRaster(
  x = passage_asc,
  filename = "./circuitscape/cs_input/passage_layer_9999.asc",
  NAflag = -9999, overwrite = TRUE
  )

###############################
# *** Leaf-scaled sea snake ***
###############################

# take cols of interest
leaf_uniq <- (unique(leaf_occ[, c(7,8,9,12,13)]))
leaf_uniq <- as.data.frame(leaf_uniq)

# create columns with long and lat values
leaf_vals <- t(apply(leaf_uniq[, 2:5], 1, function(x) x[!is.na(x)]))
leaf_uniq <- cbind(leaf_uniq, leaf_vals)
leaf_uniq <- leaf_uniq[,c(1,7,6)]
write_tsv(leaf_uniq, file = "./circuitscape/cs_input/leaf_uniq.txt")

# Change localities/regions to focal node numbers
# *** Take note that not NAs are Ashmore Reef 
# *** Edit in Notepad++

# after editing outside R for locality info of NAs, read in txt file
leaf_uniq_edit <- read_tsv("./circuitscape/cs_input/leaf_uniq_edit.txt", col_names = F)

# average coordinates of samples
mean_SB <- colMeans(as.matrix(leaf_uniq_edit[leaf_uniq_edit$X1 == "Shark Bay", 2:3]))
mean_EG <- colMeans(as.matrix(leaf_uniq_edit[leaf_uniq_edit$X1 == "Exmouth Gulf", 2:3]))
mean_PB <- colMeans(as.matrix(leaf_uniq_edit[leaf_uniq_edit$X1 == "Pilbara", 2:3]))
mean_AR <- colMeans(as.matrix(leaf_uniq_edit[leaf_uniq_edit$X1 == "Ashmore Reef", 2:3]))

# prepare focal node file
leaf_mean_loc <- rbind(mean_SB,mean_EG,mean_PB,mean_AR)
leaf_mean_loc <- cbind(X1 = c(1,2,3,4), leaf_mean_loc)
leaf_mean_loc <- as.data.frame(leaf_mean_loc)
write_tsv(leaf_mean_loc, file = "leaf_mean_loc.txt", col_names = F)

# read tsv
leaf_mean_pts <- read_tsv("./circuitscape/cs_input/leaf_mean_loc.txt", col_names = F)
leaf_mean_pts <- (as.data.frame(leaf_mean_pts[,c(2,3)]))
colnames(leaf_mean_pts) <- c("x","y")
leaf_mean_pts <- st_as_sf(leaf_mean_pts, coords = c("x", "y"), crs = 4326)

# The input for Circuitscape are as follows:
# *** 1) raster (conductance): passage_layer_9999.asc
# *** 2) focal node: leaf_mean_loc.txt (1 = SB, 2 = EG, 3 = PB, 4 = AR)

# In Circuitscape: Under options, compute log-transform current map

# plot results
leaf_main <- rast("./circuitscape/cs_output/leaf_cs_main/leaf_cs_main_cum_curmap.asc")
ggplot() +
  geom_spatraster(data = leaf_main) +
  scale_fill_continuous(
    palette = c("black","grey20","grey40","grey60","cyan"), 
    na.value = "transparent") +
  geom_sf(data = leaf_mean_pts, col = "cyan") +
  theme_bw()

###############################
# *** Short-nosed sea snake ***
###############################

# take cols of interest
short_uniq <- (unique(short_occ[, c(7,8,9,12,13)]))
short_uniq <- as.data.frame(short_uniq)

# create columns with long and lat values
short_vals <- t(apply(short_uniq[, 2:5], 1, function(x) x[!is.na(x)]))
short_uniq <- cbind(short_uniq, short_vals)
short_uniq <- short_uniq[,c(1,7,6)]
write_tsv(short_uniq, file = "./circuitscape/cs_input/short_uniq.txt")

# Change localities/regions to focal node numbers
# *** Take note that not NAs are Ashmore Reef 
# *** Edit in Notepad++

# after editing outside R for locality info of NAs, read in txt file
short_uniq_edit <- read_tsv("./circuitscape/cs_input/short_uniq_edit.txt", col_names = F)

# average coordinates of samples
ave_SB <- colMeans(as.matrix(short_uniq_edit[short_uniq_edit$X1 == "Shark Bay", 2:3]))
ave_EG <- colMeans(as.matrix(short_uniq_edit[short_uniq_edit$X1 == "Exmouth Gulf", 2:3]))
ave_PB <- colMeans(as.matrix(short_uniq_edit[short_uniq_edit$X1 == "Pilbara", 2:3]))
ave_BM <- colMeans(as.matrix(short_uniq_edit[short_uniq_edit$X1 == "Broome", 2:3]))
ave_AR <- colMeans(as.matrix(short_uniq_edit[short_uniq_edit$X1 == "Ashmore Reef", 2:3]))
ave_AB <- colMeans(as.matrix(short_uniq_edit[short_uniq_edit$X1 == "Anson Bay", 2:3]))

# prepare focal node file
short_mean_loc <- rbind(ave_SB,ave_EG,ave_PB,ave_BM,ave_AR,ave_AB)

short_mean_loc <- cbind(X1 = c(1,2,3,4,5,6), short_mean_loc)
short_mean_loc <- as.data.frame(short_mean_loc)
write_tsv(short_mean_loc, file = "./circuitscape/cs_input/short_mean_loc.txt", col_names = F)

# read tsv
short_mean_pts <- read_tsv("./circuitscape/cs_input/short_mean_loc.txt", col_names = F)
short_mean_pts <- (as.data.frame(short_mean_pts[,c(2,3)]))
colnames(short_mean_pts) <- c("x","y")
short_mean_pts <- st_as_sf(short_mean_pts, coords = c("x", "y"), crs = 4326)

# plot results
short_main <- rast("./circuitscape/cs_output/short_cs_main/short_cs_main_cum_curmap.asc")
ggplot() +
  geom_spatraster(data = short_main) +
  scale_fill_continuous(
    palette = c("black","grey20","grey40","grey60","cyan"), 
    na.value = "transparent") +
  geom_sf(data = short_mean_pts, col = "cyan") +
  theme_bw()



