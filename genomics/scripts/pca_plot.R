# Function to produce PCA plot after running VCF in plink2
pca.plot <- function(eigenvec,eigenval,species){
  
  require(tidyverse)
  
  # read in eigenvector file
  vec <- read_table(eigenvec, col_names = TRUE)
  
  # read in eigenvalues file
  val <- scan(eigenval)
  
  # remove nuisance column
  vec <- vec[,-1]
  names(vec)[1] <- "ind"
  names(vec)[2:ncol(vec)] <- paste0("PC", 1:(ncol(vec)-1))
  vec
  
  # supply sample locations
  if(species == "foliosquama") {
    Location <- rep(NA, length(vec$ind))
    Location[grep("fo", vec$ind)] <- "Ashmore Reef"
    Location[grep("401", vec$ind)] <- "Shark Bay"
    Location[grep("1014", vec$ind)] <- "Pilbara"
    Location[grep("1202", vec$ind)] <- "Pilbara"
  }

  if(species == "apraefrontalis") {
    Location <- rep(NA, length(vec$ind))
    Location[grep("Aaprae", vec$ind)] <- "Ashmore Reef"
    Location[grep("SS", vec$ind)] <- "Pilbara"
    Location[grep("KLS", vec$ind)] <- "Exmouth Gulf"
  }
  
  vec <- as.tibble(data.frame(vec, Location))
  
  pve <- data.frame(PC = 1:ncol(vec %>% select(where(is.numeric))), pve = val/sum(val)*100)
  
  # generate PCA plot
  ggplot(vec, aes(PC1, PC2, col = Location)) + 
    geom_point(size = 3) +
    coord_fixed(xlim = NULL, ylim = NULL, ratio = 1, expand = TRUE, clip = "on") +
    theme_light() +
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
  
}
