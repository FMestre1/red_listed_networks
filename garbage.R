
#################################################
# process in parallel
library(doParallel) 
#library(spatialEco)
detectCores()
cl <- makeCluster(4, type='PSOCK') # number of cores adjusted to the total number (maybe not use all!)
registerDoParallel(cl)
#################### COMPARE #################### 

?raster.modified.ttest

#IVI
ivi_compare <- raster.modified.ttest(
  ivi_nt_spatial_raster, 
  ivi_t_spatial_raster, 
  sample = "random",
  p = 0.1
)

save(ivi_compare, file = "ivi_compare.RData")

#CENTRALITY
centrality_compare <- raster.modified.ttest(
  centrality_nt_spatial_raster, 
  centrality_t_spatial_raster, 
  sample = "random",
  p = 0.1
)

save(centrality_compare, file = "centrality_compare.RData")

#IN-DEGREE
indegree_compare <- raster.modified.ttest(
  indegree_nt_spatial_raster, 
  indegree_t_spatial_raster, 
  sample = "random",
  p = 0.1
)

save(indegree_compare, file = "indegree_compare.RData")


#OUT-DEGREE
outdegree_compare <- raster.modified.ttest(
  outdegree_nt_spatial_raster, 
  outdegree_t_spatial_raster, 
  sample = "random",
  p = 0.1
)

save(outdegree_compare, file = "outdegree_compare.RData")

#CLOSENNESS
closeness_compare <- raster.modified.ttest(
  closeness_nt_spatial_raster, 
  closeness_t_spatial_raster, 
  sample = "random",
  p = 0.1
)

save(closeness_compare, file = "closeness_compare.RData")

#TROPHIC LEVEL
tl_compare <- raster.modified.ttest(
  tl_nt_spatial_raster, 
  tl_t_spatial_raster, 
  sample = "random",
  p = 0.1
)

save(tl_compare, file = "tl_compare.RData")


################################################################################

#Use smaller area to test function

portugal <- terra::vect("C:/Users/mestr/Documents/0. Artigos/IUCN_networks/shapefiles/portugal.shp")
#terra::crs(portugal)
portugal_p <- terra::project(portugal, ivi_nt_spatial_raster)
#terra::crs(portugal_p)
#plot(portugal)
#terra::crs(ivi_nt_spatial_raster)

#plot(ivi_nt_spatial_raster)
#plot(portugal_p, add=TRUE)

ivi_nt_spatial_raster_PT <- terra::crop(ivi_nt_spatial_raster, portugal_p)
#plot(ivi_nt_spatial_raster_PT)
ivi_t_spatial_raster_PT <- terra::crop(ivi_t_spatial_raster, portugal_p)
#plot(ivi_t_spatial_raster_PT)

#ivi_compare <- raster.modified.ttest(
#  ivi_nt_spatial_raster_PT, 
#  ivi_t_spatial_raster_PT, 
#  sample = "random",
#  p = 0.01
#)

####
#Compare maps
#plot(terra::compare(ivi_nt_spatial_raster_PT, ivi_t_spatial_raster_PT, oper = "<"))
#terra::compareGeom(ivi_nt_spatial_raster_PT, ivi_t_spatial_raster_PT)

#Using diffeR
#?terra::compare

#test <- diffeR::differenceMR(ivi_nt_spatial_raster_PT, ivi_t_spatial_raster_PT, eval = "original")
#test <- diffeR::differenceMR(ivi_nt_spatial_raster, ivi_t_spatial_raster)
#test

# process in parallel
#library(doParallel) 
#detectCores()
#cl <- makeCluster(4, type='PSOCK') # number of cores adjusted to the total number (maybe not use all!)
#registerDoParallel(cl)

#test_ivi <- spatialEco::raster.change(ivi_nt_spatial_raster_PT, ivi_t_spatial_raster_PT, stat = "t.test")


#################

#cor.test(as.vector(ivi_nt_spatial_raster_PT), as.vector(ivi_t_spatial_raster_PT))

#################



library(ade4)

tl_nt_spatial_raster
tl_t_spatial_raster

# Reduce resolution by a factor of 2 (or more if needed)
tl_nt_spatial_raster_resampled <- aggregate(tl_nt_spatial_raster, fact = 2)
tl_t_spatial_raster_resampled <- aggregate(tl_t_spatial_raster, fact = 2)

tl_ssi_compare_resample <- SSIMmap::ssim_raster(tl_nt_spatial_raster_resampled, 
                                                tl_t_spatial_raster_resampled, 
                                                global = FALSE)

plot_tl_SSIM <- rasterVis::levelplot(tl_ssi_compare_resample[[1]], contour=FALSE, par.settings = myTheme, main= "SSIM", margin = FALSE)
plot_tl_SIM <- rasterVis::levelplot(tl_ssi_compare_resample[[2]], contour=FALSE, par.settings = myTheme, main= "SIM", margin = FALSE)
plot_tl_SIV <- rasterVis::levelplot(tl_ssi_compare_resample[[3]], contour=FALSE, par.settings = myTheme, main= "SIV", margin = FALSE)
plot_tl_SIP <- rasterVis::levelplot(tl_ssi_compare_resample[[4]], contour=FALSE, par.settings = myTheme, main= "SIP", margin = FALSE)

gridExtra::grid.arrange(plot_tl_SSIM,
                        plot_tl_SIM,
                        plot_tl_SIV,
                        plot_tl_SIP,
                        ncol=2,
                        top=grid::textGrob("Trophic Level"))

# Assuming map1 and map2 are matrices or rasters
#dist_map1 <- dist(as.vector(tl_nt_spatial_raster_resampled))
#dist_map2 <- dist(as.vector(tl_t_spatial_raster_resampled))

# Perform Mantel test
mantel_result <- mantel.rtest(dist_map1, dist_map2, nrepet = 999)

##############################################

library(spdep)

# Create a spatial neighborhood structure (based on map coordinates)
coords <- coordinates(map1)  # Assuming map1 is a raster
nb <- dnearneigh(coords, 0, 1.5)  # Define neighbors within a certain distance
lw <- nb2listw(nb, style = "W")

# Compute the actual correlation or similarity
actual_similarity <- cor(map1_vals, map2_vals)

# Perform a spatially constrained permutation test
n_permutations <- 999
similarity_scores <- numeric(n_permutations)

for (i in 1:n_permutations) {
  # Perform spatially constrained permutation
  permuted_map2_vals <- spsample(map2_vals, lw)
  
  # Compute the similarity for the permuted map
  similarity_scores[i] <- cor(map1_vals, permuted_map2_vals)
}

# Calculate the p-value
p_value <- sum(similarity_scores >= actual_similarity) / n_permutations

print(paste("Spatially constrained permutation test p-value:", p_value))

################################################################################
################################################################################

# Load required packages
library(terra)
library(spdep)
library(ggplot2)
library(gridExtra)
library(viridis)
library(doParallel)

# Load your maps
tl_nt_spatial_raster <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\transformed\\tl_nt_spatial_raster_z_raster.tif")
tl_t_spatial_raster <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\transformed\\tl_t_spatial_raster_z_raster.tif")
ivi_nt_spatial_raster <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\transformed\\ivi_nt_spatial_raster_z_raster.tif")
ivi_t_spatial_raster <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\transformed\\ivi_t_spatial_raster_z_raster.tif")
centrality_nt_spatial_raster <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\transformed\\centrality_nt_spatial_raster_z_raster.tif")
centrality_t_spatial_raster <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\transformed\\centrality_t_spatial_raster_z_raster.tif")
closeness_nt_spatial_raster <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\transformed\\closeness_nt_spatial_raster_z_raster.tif")
closeness_t_spatial_raster <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\transformed\\closeness_t_spatial_raster_z_raster.tif")
indegree_nt_spatial_raster <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\transformed\\indegree_nt_spatial_raster_z_raster.tif")
indegree_t_spatial_raster <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\transformed\\indegree_t_spatial_raster_z_raster.tif")
outdegree_nt_spatial_raster <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\transformed\\outdegree_nt_spatial_raster_z_raster.tif")
outdegree_t_spatial_raster <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\transformed\\outdegree_t_spatial_raster_z_raster.tif")

# Ensure both maps have the same extent and resolution
#map2 <- resample(map2, map1)

# Detect number of cores
ncores <- detectCores() - 1  # Leave one core free for system tasks

# Set up and register cluster
cl <- makeCluster(ncores)
registerDoParallel(cl)

#######

# Parallel computation using foreach
results <- foreach(i = 1:100) %dopar% kappa <- kappa_calc(tl_nt_spatial_raster, tl_t_spatial_raster)

# Stop cluster when done
stopCluster(cl)

## 1. Kappa Statistic
# Note: terra doesn't have a direct kappa function, so we'll use a custom implementation
kappa_calc <- function(map1, map2) {
  confusion_matrix <- table(values(map1), values(map2))
  n <- sum(confusion_matrix)
  nc <- ncol(confusion_matrix)
  diag <- diag(confusion_matrix)
  rowsums <- rowSums(confusion_matrix)
  colsums <- colSums(confusion_matrix)
  expected <- rowsums %*% t(colsums) / n
  
  kappa <- (sum(diag) - sum(expected)) / (n - sum(expected))
  return(kappa)
}

kappa <- kappa_calc(tl_nt_spatial_raster, tl_t_spatial_raster)
print(paste("Kappa Statistic:", kappa))

## 2. Moran's I Autocorrelation Index
pts <- as.points(map1)
nb <- dnearneigh(crds(pts), 0, res(map1)[1])
lw <- nb2listw(nb, style = "W")

moran_map1 <- moran.test(values(map1), lw)
moran_map2 <- moran.test(values(map2), lw)

print("Moran's I for Map 1:")
print(moran_map1)
print("Moran's I for Map 2:")
print(moran_map2)

## 3. Raster Subtraction
diff_map <- map1 - map2

# Visualize the difference
plot(diff_map, main = "Difference between Map 1 and Map 2")

## 4. Exploratory Spatial Data Analysis (ESDA)
# Local Moran's I
local_moran_map1 <- localmoran(values(map1), lw)
local_moran_map2 <- localmoran(values(map2), lw)

# Visualize Local Moran's I
local_moran1 <- rast(map1)
values(local_moran1) <- local_moran_map1[, 1]
local_moran2 <- rast(map2)
values(local_moran2) <- local_moran_map2[, 1]

par(mfrow = c(1, 2))
plot(local_moran1, main = "Local Moran's I - Map 1")
plot(local_moran2, main = "Local Moran's I - Map 2")

## 5. Statistical Comparison
# Perform a paired t-test
t_test_result <- t.test(values(map1), values(map2), paired = TRUE)
print("Paired t-test result:")
print(t_test_result)

## 6. Visualization
# Create a function to plot rasters
plot_raster <- function(raster, title) {
  df <- as.data.frame(raster, xy = TRUE)
  names(df) <- c("x", "y", "value")
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_vi
}  
