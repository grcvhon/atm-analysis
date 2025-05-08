library(tidyverse)

var_based_stats <- function(lqual,
                            ldepth.mean,
                            lmiss,
                            idepth,
                            imiss){
  # variant quality
  var_qual <- read_delim(lqual, delim = "\t",
                         col_names = c("chr", "pos", "qual"), skip = 1)
  a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
  
  # variant mean depth
  var_depth <- read_delim(ldepth.mean, delim = "\t",
                          col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
  b <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
  
  # variant missingness
  var_miss <- read_delim(lmiss, delim = "\t",
                         col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  c <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
  
  # mean depth ind
  ind_depth <- read_delim(idepth, delim = "\t",
                          col_names = c("ind", "nsites", "depth"), skip = 1)
  
  d <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
  
  # missing ind
  ind_miss  <- read_delim(imiss, delim = "\t",
                          col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
  e <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
  
  #a + theme_light()
  b + theme_light()
  #c + theme_light()
  #d + theme_light()
  #e + theme_light()
  
}

var_based_stats(lqual = "./genomics/vcftools/AFO-reference.lqual",
                ldepth.mean = "./genomics/vcftools/AFO-reference.ldepth.mean",
                lmiss = "./genomics/vcftools/AFO-reference.lmiss",
                idepth = "./genomics/vcftools/AFO-reference.idepth",
                imiss = "./genomics/vcftools/AFO-denovo.imiss")
