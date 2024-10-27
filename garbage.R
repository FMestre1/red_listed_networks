
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


##############################################################################################
##############################################################################################


################################################################################
#                               Using spatialEco
################################################################################

#test_ivi <- spatialEco::raster.change(ivi_nt_spatial_raster, ivi_t_spatial_raster, stat = "t.test")
#test_ivi <- spatialEco::raster.change(ivi_nt_spatial_raster, ivi_t_spatial_raster, stat = "t.test")
#test_centrality <- spatialEco::raster.change(centrality_nt_spatial_raster, centrality_t_spatial_raster, stat = "t.test")
#test_closeness <- spatialEco::raster.change(closeness_nt_spatial_raster, closeness_t_spatial_raster, stat = "t.test")
#test_indegree <- spatialEco::raster.change(indegree_nt_spatial_raster, indegree_t_spatial_raster, stat = "t.test")
#test_outdegree <- spatialEco::raster.change(outdegree_nt_spatial_raster, outdegree_t_spatial_raster, stat = "t.test")
#test_tl <- spatialEco::raster.change(tl_nt_spatial_raster, tl_t_spatial_raster, stat = "t.test")

#Save ...
#terra::writeRaster(test_ivi[[1]], filename = "ivi_raster_change_1.tif")
#terra::writeRaster(test_ivi[[2]], filename = "ivi_raster_change_2.tif")

#terra::writeRaster(test_centrality[[1]], filename = "centrality_raster_change_1.tif")
#terra::writeRaster(test_centrality[[2]], filename = "centrality_raster_change_2.tif")

#terra::writeRaster(test_closeness[[1]], filename = "closeness_raster_change_1.tif")
#terra::writeRaster(test_closeness[[2]], filename = "closeness_raster_change_2.tif")

#terra::writeRaster(test_indegree[[1]], filename = "indegree_raster_change_1.tif")
#terra::writeRaster(test_indegree[[2]], filename = "indegree_raster_change_2.tif")

#terra::writeRaster(test_outdegree[[1]], filename = "outdegree_raster_change_1.tif")
#terra::writeRaster(test_outdegree[[2]], filename = "outdegree_raster_change_2.tif")

#terra::writeRaster(test_tl[[1]], filename = "tl_raster_change_1.tif")
#terra::writeRaster(test_tl[[2]], filename = "tl_raster_change_2.tif")

################################################################################
#                                 Using SSIMmap
################################################################################

##### Previously apply a z score transformation to every raster #####

############### 1. IVI ###############

# Calculate the mean and standard deviation
#mean
ivi_nt_spatial_raster_mean <- as.numeric(global(ivi_nt_spatial_raster, mean, na.rm = TRUE))
ivi_t_spatial_raster_mean <- as.numeric(global(ivi_t_spatial_raster, mean, na.rm = TRUE))
#sd
ivi_nt_spatial_raster_sd <- as.numeric(global(ivi_nt_spatial_raster, sd, na.rm = TRUE))
ivi_t_spatial_raster_sd <- as.numeric(global(ivi_t_spatial_raster, sd, na.rm = TRUE))

print(ivi_nt_spatial_raster_mean)
print(ivi_t_spatial_raster_mean)
print(ivi_nt_spatial_raster_sd)
print(ivi_t_spatial_raster_sd)

# Apply Z-score transformation
ivi_nt_spatial_raster_z_raster <- (ivi_nt_spatial_raster - ivi_nt_spatial_raster_mean)/ivi_nt_spatial_raster_sd
ivi_t_spatial_raster_z_raster <- (ivi_t_spatial_raster - ivi_t_spatial_raster_mean)/ivi_t_spatial_raster_sd

# Plot the Z-score transformed raster
#plot(ivi_nt_spatial_raster_z_raster, main = "IVI NT - Z-score Transformed")
#plot(ivi_t_spatial_raster_z_raster, main = "IVI T - Z-score Transformed")

# Save the transformed rasters
writeRaster(ivi_nt_spatial_raster_z_raster, "ivi_nt_spatial_raster_z_raster.tif", overwrite = TRUE)
writeRaster(ivi_t_spatial_raster_z_raster, "ivi_t_spatial_raster_z_raster.tif", overwrite = TRUE)

############### 2. TROPHIC LEVEL ###############

# Calculate the mean and standard deviation
#mean
tl_nt_spatial_raster_mean <- as.numeric(global(tl_nt_spatial_raster, mean, na.rm = TRUE))
tl_t_spatial_raster_mean <- as.numeric(global(tl_t_spatial_raster, mean, na.rm = TRUE))
#sd
tl_nt_spatial_raster_sd <- as.numeric(global(tl_nt_spatial_raster, sd, na.rm = TRUE))
tl_t_spatial_raster_sd <- as.numeric(global(tl_t_spatial_raster, sd, na.rm = TRUE))

print(tl_nt_spatial_raster_mean)
print(tl_t_spatial_raster_mean)
print(tl_nt_spatial_raster_sd)
print(tl_t_spatial_raster_sd)

# Apply Z-score transformation
tl_nt_spatial_raster_z_raster <- (tl_nt_spatial_raster - tl_nt_spatial_raster_mean)/tl_nt_spatial_raster_sd
tl_t_spatial_raster_z_raster <- (tl_t_spatial_raster - tl_t_spatial_raster_mean)/tl_t_spatial_raster_sd

# Plot the Z-score transformed raster
#plot(tl_nt_spatial_raster_z_raster, main = "TL NT - Z-score Transformed")
#plot(tl_t_spatial_raster_z_raster, main = "TL T - Z-score Transformed")

# Save the transformed rasters
writeRaster(tl_nt_spatial_raster_z_raster, "tl_nt_spatial_raster_z_raster.tif", overwrite = TRUE)
writeRaster(tl_t_spatial_raster_z_raster, "tl_t_spatial_raster_z_raster.tif", overwrite = TRUE)

############### 3. CENTRALITY ###############

# Calculate the mean and standard deviation
#mean
centrality_nt_spatial_raster_mean <- as.numeric(global(centrality_nt_spatial_raster, mean, na.rm = TRUE))
centrality_t_spatial_raster_mean <- as.numeric(global(centrality_t_spatial_raster, mean, na.rm = TRUE))
#sd
centrality_nt_spatial_raster_sd <- as.numeric(global(centrality_nt_spatial_raster, sd, na.rm = TRUE))
centrality_t_spatial_raster_sd <- as.numeric(global(centrality_t_spatial_raster, sd, na.rm = TRUE))

print(centrality_nt_spatial_raster_mean)
print(centrality_t_spatial_raster_mean)
print(centrality_nt_spatial_raster_sd)
print(centrality_t_spatial_raster_sd)

# Apply Z-score transformation
centrality_nt_spatial_raster_z_raster <- (centrality_nt_spatial_raster - centrality_nt_spatial_raster_mean)/centrality_nt_spatial_raster_sd
centrality_t_spatial_raster_z_raster <- (centrality_t_spatial_raster - centrality_t_spatial_raster_mean)/centrality_t_spatial_raster_sd

# Plot the Z-score transformed raster
#plot(centrality_nt_spatial_raster_z_raster, main = "Centrality NT - Z-score Transformed")
#plot(centrality_t_spatial_raster_z_raster, main = "Centrality T - Z-score Transformed")

# Save the transformed rasters
writeRaster(centrality_nt_spatial_raster_z_raster, "centrality_nt_spatial_raster_z_raster.tif", overwrite = TRUE)
writeRaster(centrality_t_spatial_raster_z_raster, "centrality_t_spatial_raster_z_raster.tif", overwrite = TRUE)

############### 4. CLOSENESS ###############

# Calculate the mean and standard deviation
#mean
closeness_nt_spatial_raster_mean <- as.numeric(global(closeness_nt_spatial_raster, mean, na.rm = TRUE))
closeness_t_spatial_raster_mean <- as.numeric(global(closeness_t_spatial_raster, mean, na.rm = TRUE))
#sd
closeness_nt_spatial_raster_sd <- as.numeric(global(closeness_nt_spatial_raster, sd, na.rm = TRUE))
closeness_t_spatial_raster_sd <- as.numeric(global(closeness_t_spatial_raster, sd, na.rm = TRUE))

print(closeness_nt_spatial_raster_mean)
print(closeness_t_spatial_raster_mean)
print(closeness_nt_spatial_raster_sd)
print(closeness_t_spatial_raster_sd)

# Apply Z-score transformation
closeness_nt_spatial_raster_z_raster <- (closeness_nt_spatial_raster - closeness_nt_spatial_raster_mean)/closeness_nt_spatial_raster_sd
closeness_t_spatial_raster_z_raster <- (closeness_t_spatial_raster - closeness_t_spatial_raster_mean)/closeness_t_spatial_raster_sd

# Plot the Z-score transformed raster
#plot(closeness_nt_spatial_raster_z_raster, main = "closeness NT - Z-score Transformed")
#plot(closeness_t_spatial_raster_z_raster, main = "closeness T - Z-score Transformed")

# Save the transformed rasters
writeRaster(closeness_nt_spatial_raster_z_raster, "closeness_nt_spatial_raster_z_raster.tif", overwrite = TRUE)
writeRaster(closeness_t_spatial_raster_z_raster, "closeness_t_spatial_raster_z_raster.tif", overwrite = TRUE)

############### 5. IN-DEGREE ###############

# Calculate the mean and standard deviation
#mean
indegree_nt_spatial_raster_mean <- as.numeric(global(indegree_nt_spatial_raster, mean, na.rm = TRUE))
indegree_t_spatial_raster_mean <- as.numeric(global(indegree_t_spatial_raster, mean, na.rm = TRUE))
#sd
indegree_nt_spatial_raster_sd <- as.numeric(global(indegree_nt_spatial_raster, sd, na.rm = TRUE))
indegree_t_spatial_raster_sd <- as.numeric(global(indegree_t_spatial_raster, sd, na.rm = TRUE))

print(indegree_nt_spatial_raster_mean)
print(indegree_t_spatial_raster_mean)
print(indegree_nt_spatial_raster_sd)
print(indegree_t_spatial_raster_sd)

# Apply Z-score transformation
indegree_nt_spatial_raster_z_raster <- (indegree_nt_spatial_raster - indegree_nt_spatial_raster_mean)/indegree_nt_spatial_raster_sd
indegree_t_spatial_raster_z_raster <- (indegree_t_spatial_raster - indegree_t_spatial_raster_mean)/indegree_t_spatial_raster_sd

# Plot the Z-score transformed raster
#plot(indegree_nt_spatial_raster_z_raster, main = "indegree NT - Z-score Transformed")
#plot(indegree_t_spatial_raster_z_raster, main = "indegree T - Z-score Transformed")

# Save the transformed rasters
writeRaster(indegree_nt_spatial_raster_z_raster, "indegree_nt_spatial_raster_z_raster.tif", overwrite = TRUE)
writeRaster(indegree_t_spatial_raster_z_raster, "indegree_t_spatial_raster_z_raster.tif", overwrite = TRUE)

############### 6. OUT-DEGREE ###############

# Calculate the mean and standard deviation
#mean
outdegree_nt_spatial_raster_mean <- as.numeric(global(outdegree_nt_spatial_raster, mean, na.rm = TRUE))
outdegree_t_spatial_raster_mean <- as.numeric(global(outdegree_t_spatial_raster, mean, na.rm = TRUE))
#sd
outdegree_nt_spatial_raster_sd <- as.numeric(global(outdegree_nt_spatial_raster, sd, na.rm = TRUE))
outdegree_t_spatial_raster_sd <- as.numeric(global(outdegree_t_spatial_raster, sd, na.rm = TRUE))

print(outdegree_nt_spatial_raster_mean)
print(outdegree_t_spatial_raster_mean)
print(outdegree_nt_spatial_raster_sd)
print(outdegree_t_spatial_raster_sd)

# Apply Z-score transformation
outdegree_nt_spatial_raster_z_raster <- (outdegree_nt_spatial_raster - outdegree_nt_spatial_raster_mean)/outdegree_nt_spatial_raster_sd
outdegree_t_spatial_raster_z_raster <- (outdegree_t_spatial_raster - outdegree_t_spatial_raster_mean)/outdegree_t_spatial_raster_sd

# Plot the Z-score transformed raster
#plot(outdegree_nt_spatial_raster_z_raster, main = "outdegree NT - Z-score Transformed")
#plot(outdegree_t_spatial_raster_z_raster, main = "outdegree T - Z-score Transformed")

# Save the transformed rasters
writeRaster(outdegree_nt_spatial_raster_z_raster, "outdegree_nt_spatial_raster_z_raster.tif", overwrite = TRUE)
writeRaster(outdegree_t_spatial_raster_z_raster, "outdegree_t_spatial_raster_z_raster.tif", overwrite = TRUE)

############### 7. COMPARE MAPS ###############

#Load required rasters (if these have been deleted from the environment)
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

#Run & Save...

ivi_ssi_compare <- SSIMmap::ssim_raster(ivi_nt_spatial_raster, ivi_t_spatial_raster, global = FALSE)
terra::writeRaster(ivi_ssi_compare[[1]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\ivi_ssi_compare_SSIM.tif")
terra::writeRaster(ivi_ssi_compare[[2]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\ivi_ssi_compare_SIM.tif")
terra::writeRaster(ivi_ssi_compare[[3]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\ivi_ssi_compare_SIV.tif")
terra::writeRaster(ivi_ssi_compare[[4]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\ivi_ssi_compare_SIP.tif")
#
tl_ssi_compare <- SSIMmap::ssim_raster(tl_nt_spatial_raster, tl_t_spatial_raster, global = FALSE)
terra::writeRaster(tl_ssi_compare[[1]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\tl_ssi_compare_SSIM.tif")
terra::writeRaster(tl_ssi_compare[[2]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\tl_ssi_compare_SIM.tif")
terra::writeRaster(tl_ssi_compare[[3]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\tl_ssi_compare_SIV.tif")
terra::writeRaster(tl_ssi_compare[[4]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\tl_ssi_compare_SIP.tif")
#
centrality_ssi_compare <- SSIMmap::ssim_raster(centrality_nt_spatial_raster, centrality_t_spatial_raster, global = FALSE)
terra::writeRaster(centrality_ssi_compare[[1]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\centrality_ssi_compare_SSIM.tif")
terra::writeRaster(centrality_ssi_compare[[2]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\centrality_ssi_compare_SIM.tif")
terra::writeRaster(centrality_ssi_compare[[3]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\centrality_ssi_compare_SIV.tif")
terra::writeRaster(centrality_ssi_compare[[4]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\centrality_ssi_compare_SIP.tif")
#
closeness_ssi_compare <- SSIMmap::ssim_raster(closeness_nt_spatial_raster, closeness_t_spatial_raster, global = FALSE)
terra::writeRaster(closeness_ssi_compare[[1]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\closeness_ssi_compare_SSIM.tif")
terra::writeRaster(closeness_ssi_compare[[2]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\closeness_ssi_compare_SIM.tif")
terra::writeRaster(closeness_ssi_compare[[3]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\closeness_ssi_compare_SIV.tif")
terra::writeRaster(closeness_ssi_compare[[4]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\closeness_ssi_compare_SIP.tif") 
#
outdegree_ssi_compare <- SSIMmap::ssim_raster(outdegree_nt_spatial_raster, outdegree_t_spatial_raster, global = FALSE)
terra::writeRaster(outdegree_ssi_compare[[1]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\outdegree_ssi_compare_SSIM.tif")
terra::writeRaster(outdegree_ssi_compare[[2]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\outdegree_ssi_compare_SIM.tif")
terra::writeRaster(outdegree_ssi_compare[[3]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\outdegree_ssi_compare_SIV.tif")
terra::writeRaster(outdegree_ssi_compare[[4]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\outdegree_ssi_compare_SIP.tif")
#
indegree_ssi_compare <- SSIMmap::ssim_raster(indegree_nt_spatial_raster, indegree_t_spatial_raster, global = FALSE)
terra::writeRaster(indegree_ssi_compare[[1]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\indegree_ssi_compare_SSIM.tif")
terra::writeRaster(indegree_ssi_compare[[2]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\indegree_ssi_compare_SIM.tif")
terra::writeRaster(indegree_ssi_compare[[3]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\indegree_ssi_compare_SIV.tif")
terra::writeRaster(indegree_ssi_compare[[4]], filename = "D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\indegree_ssi_compare_SIP.tif")

################################################################################
#                                  Plotting
################################################################################

#https://cran.r-project.org/web/packages/SSIMmap/vignettes/Introduction_to_SSIMmap.html

ivi_ssi_compare_SSIM <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\ivi_ssi_compare_SSIM.tif")
ivi_ssi_compare_SIM <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\ivi_ssi_compare_SIM.tif")
ivi_ssi_compare_SIV <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\ivi_ssi_compare_SIV.tif")
ivi_ssi_compare_SIP <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\ivi_ssi_compare_SIP.tif")
#
centrality_ssi_compare_SSIM <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\centrality_ssi_compare_SSIM.tif")
centrality_ssi_compare_SIM <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\centrality_ssi_compare_SIM.tif")
centrality_ssi_compare_SIV <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\centrality_ssi_compare_SIV.tif")
centrality_ssi_compare_SIP <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\centrality_ssi_compare_SIP.tif")
#
outdegree_ssi_compare_SSIM <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\outdegree_ssi_compare_SSIM.tif")
outdegree_ssi_compare_SIM <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\outdegree_ssi_compare_SIM.tif")
outdegree_ssi_compare_SIV <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\outdegree_ssi_compare_SIV.tif")
outdegree_ssi_compare_SIP <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\outdegree_ssi_compare_SIP.tif")
#
indegree_ssi_compare_SSIM <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\indegree_ssi_compare_SSIM.tif")
indegree_ssi_compare_SIM <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\indegree_ssi_compare_SIM.tif")
indegree_ssi_compare_SIV <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\indegree_ssi_compare_SIV.tif")
indegree_ssi_compare_SIP <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\indegree_ssi_compare_SIP.tif")
#
closeness_ssi_compare_SSIM <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\closeness_ssi_compare_SSIM.tif")
closeness_ssi_compare_SIM <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\closeness_ssi_compare_SIM.tif")
closeness_ssi_compare_SIV <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\closeness_ssi_compare_SIV.tif")
closeness_ssi_compare_SIP <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\closeness_ssi_compare_SIP.tif")
#
tl_ssi_compare_SSIM <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\tl_ssi_compare_SSIM.tif")
tl_ssi_compare_SIM <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\tl_ssi_compare_SIM.tif")
tl_ssi_compare_SIV <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\tl_ssi_compare_SIV.tif")
tl_ssi_compare_SIP <- terra::rast("D:\\THREATENED_NON_THREATENED_ SPECIES\\output\\tl_ssi_compare_SIP.tif")

##

colours_RYB <- colorRampPalette(c("darkred" ,"lightgrey","darkgreen"))
myTheme <- rasterVis::rasterTheme(region = colours_RYB(100))

plot_ivi_SSIM <- rasterVis::levelplot(ivi_ssi_compare_SSIM, contour=FALSE, par.settings = myTheme, main= "SSIM", margin = FALSE)
plot_ivi_SIM <- rasterVis::levelplot(ivi_ssi_compare_SIM, contour=FALSE, par.settings = myTheme, main= "SIM", margin = FALSE)
plot_ivi_SIV <- rasterVis::levelplot(ivi_ssi_compare_SIV, contour=FALSE, par.settings = myTheme, main= "SIV", margin = FALSE)
plot_ivi_SIP <- rasterVis::levelplot(ivi_ssi_compare_SIP, contour=FALSE, par.settings = myTheme, main= "SIP", margin = FALSE)

gridExtra::grid.arrange(plot_ivi_SSIM,
                        plot_ivi_SIM,
                        plot_ivi_SIV,
                        plot_ivi_SIP,
                        ncol=2,
                        top=grid::textGrob("IVI"))

##

plot_centrality_SSIM <- rasterVis::levelplot(centrality_ssi_compare_SSIM, contour=FALSE, par.settings = myTheme, main= "SSIM", margin = FALSE)
plot_centrality_SIM <- rasterVis::levelplot(centrality_ssi_compare_SIM, contour=FALSE, par.settings = myTheme, main= "SIM", margin = FALSE)
plot_centrality_SIV <- rasterVis::levelplot(centrality_ssi_compare_SIV, contour=FALSE, par.settings = myTheme, main= "SIV", margin = FALSE)
plot_centrality_SIP <- rasterVis::levelplot(centrality_ssi_compare_SIP, contour=FALSE, par.settings = myTheme, main= "SIP", margin = FALSE)

gridExtra::grid.arrange(plot_centrality_SSIM,
                        plot_centrality_SIM,
                        plot_centrality_SIV,
                        plot_centrality_SIP,
                        ncol=2,
                        top=grid::textGrob("Centrality"))


##

plot_tl_SSIM <- rasterVis::levelplot(tl_ssi_compare_SSIM, contour=FALSE, par.settings = myTheme, main= "SSIM", margin = FALSE)
plot_tl_SIM <- rasterVis::levelplot(tl_ssi_compare_SIM, contour=FALSE, par.settings = myTheme, main= "SIM", margin = FALSE)
plot_tl_SIV <- rasterVis::levelplot(tl_ssi_compare_SIV, contour=FALSE, par.settings = myTheme, main= "SIV", margin = FALSE)
plot_tl_SIP <- rasterVis::levelplot(tl_ssi_compare_SIP, contour=FALSE, par.settings = myTheme, main= "SIP", margin = FALSE)

gridExtra::grid.arrange(plot_tl_SSIM,
                        plot_tl_SIM,
                        plot_tl_SIV,
                        plot_tl_SIP,
                        ncol=2,
                        top=grid::textGrob("Trophic Level"))

##

plot_outdegree_SSIM <- rasterVis::levelplot(outdegree_ssi_compare_SSIM, contour=FALSE, par.settings = myTheme, main= "SSIM", margin = FALSE)
plot_outdegree_SIM <- rasterVis::levelplot(outdegree_ssi_compare_SIM, contour=FALSE, par.settings = myTheme, main= "SIM", margin = FALSE)
plot_outdegree_SIV <- rasterVis::levelplot(outdegree_ssi_compare_SIV, contour=FALSE, par.settings = myTheme, main= "SIV", margin = FALSE)
plot_outdegree_SIP <- rasterVis::levelplot(outdegree_ssi_compare_SIP, contour=FALSE, par.settings = myTheme, main= "SIP", margin = FALSE)

gridExtra::grid.arrange(plot_outdegree_SSIM,
                        plot_outdegree_SIM,
                        plot_outdegree_SIV,
                        plot_outdegree_SIP,
                        ncol=2,
                        top=grid::textGrob("Outdegree"))

##

plot_indegree_SSIM <- rasterVis::levelplot(indegree_ssi_compare_SSIM, contour=FALSE, par.settings = myTheme, main= "SSIM", margin = FALSE)
plot_indegree_SIM <- rasterVis::levelplot(indegree_ssi_compare_SIM, contour=FALSE, par.settings = myTheme, main= "SIM", margin = FALSE)
plot_indegree_SIV <- rasterVis::levelplot(indegree_ssi_compare_SIV, contour=FALSE, par.settings = myTheme, main= "SIV", margin = FALSE)
plot_indegree_SIP <- rasterVis::levelplot(indegree_ssi_compare_SIP, contour=FALSE, par.settings = myTheme, main= "SIP", margin = FALSE)

gridExtra::grid.arrange(plot_indegree_SSIM,
                        plot_indegree_SIM,
                        plot_indegree_SIV,
                        plot_indegree_SIP,
                        ncol=2,
                        top=grid::textGrob("Indegree"))

##

plot_closeness_SSIM <- rasterVis::levelplot(closeness_ssi_compare_SSIM, contour=FALSE, par.settings = myTheme, main= "SSIM", margin = FALSE)
plot_closeness_SIM <- rasterVis::levelplot(closeness_ssi_compare_SIM, contour=FALSE, par.settings = myTheme, main= "SIM", margin = FALSE)
plot_closeness_SIV <- rasterVis::levelplot(closeness_ssi_compare_SIV, contour=FALSE, par.settings = myTheme, main= "SIV", margin = FALSE)
plot_closeness_SIP <- rasterVis::levelplot(closeness_ssi_compare_SIP, contour=FALSE, par.settings = myTheme, main= "SIP", margin = FALSE)

gridExtra::grid.arrange(plot_closeness_SSIM,
                        plot_closeness_SIM,
                        plot_closeness_SIV,
                        plot_closeness_SIP,
                        ncol=2,
                        top=grid::textGrob("Closeness"))

################################################################################
#                             raster.modified.ttest
################################################################################

#FMestre
#29-09-2024

library(spatialEco)

#Load un-transformed rasters
ivi_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_ivi_15JUL.tif")
ivi_t_spatial_raster <- terra::rast("rasters_15JUL\\t_ivi_15JUL.tif")
centrality_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_centrality_15JUL.tif")
centrality_t_spatial_raster <- terra::rast("rasters_15JUL\\t_centrality_15JUL.tif")
outdegree_t_spatial_raster <- terra::rast("rasters_15JUL\\t_outdegree_15JUL.tif")
outdegree_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_outdegree_15JUL.tif")
indegree_t_spatial_raster <- terra::rast("rasters_15JUL\\t_indegree_15JUL.tif")
indegree_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_indegree_15JUL.tif")
closeness_t_spatial_raster <- terra::rast("rasters_15JUL\\t_closeness_15JUL.tif")
closeness_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_closeness_15JUL.tif")
tl_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_tl_15JUL.tif")
tl_t_spatial_raster <- terra::rast("rasters_15JUL\\t_tl_15JUL.tif")

#Load transformed rasters
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

#Delete takes too long...
#tl_ttest <- raster.modified.ttest(tl_nt_spatial_raster, tl_t_spatial_raster, d = "auto", sample = "regular", p = 0.0001)

################################################################################
#                             PEARSON CORRELATION
################################################################################

#?terra::layerCor
cor_tl <- terra::layerCor(c(tl_nt_spatial_raster, tl_t_spatial_raster), fun = "cor")
cor_indegree <- terra::layerCor(c(indegree_nt_spatial_raster, indegree_t_spatial_raster), fun = "cor")
cor_outdegree <- terra::layerCor(c(outdegree_nt_spatial_raster, outdegree_t_spatial_raster), fun = "cor")
cor_centrality <- terra::layerCor(c(centrality_nt_spatial_raster, centrality_t_spatial_raster), fun = "cor")
cor_closeness <- terra::layerCor(c(closeness_nt_spatial_raster, closeness_t_spatial_raster), fun = "cor")
cor_ivi <- terra::layerCor(c(ivi_nt_spatial_raster, ivi_t_spatial_raster), fun = "cor")
#
cor_tl$correlation[1,2]
cor_indegree$correlation[1,2]
cor_outdegree$correlation[1,2]
cor_centrality$correlation[1,2]
cor_closeness$correlation[1,2]
cor_ivi$correlation[1,2]


################################################################################
#                                    GLMER
################################################################################

#Load packages
library(lme4)

??glmer
#lme4::glmer(indices ~ fator+ (1 | lat:long)

# Extract the values from the raster
ivi_nt_spatial_raster_values <- terra::values(ivi_nt_spatial_raster)
ivi_t_spatial_raster_values <- terra::values(ivi_t_spatial_raster)

# Get the coordinates of the raster cells
ivi_nt_spatial_raster_coords <- terra::crds(ivi_nt_spatial_raster)
ivi_t_spatial_raster_coords <- terra::crds(ivi_t_spatial_raster)

# Create a data frame with lat, long, and rich columns
ivi_nt_df <- data.frame(lat = ivi_nt_spatial_raster_coords[, 2], long = ivi_nt_spatial_raster_coords[, 1], index = ivi_nt_spatial_raster_values)
ivi_t_df <- data.frame(lat = ivi_t_spatial_raster_coords[, 2], long = ivi_t_spatial_raster_coords[, 1], index = ivi_t_spatial_raster_values)

# Print the data frame
head(ivi_nt_df)
head(ivi_t_df)


################################################################################
#                                    diffeR
################################################################################

#Load packages
library(diffeR)

#Load rasters
ivi_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_ivi_15JUL.tif")
ivi_t_spatial_raster <- terra::rast("rasters_15JUL\\t_ivi_15JUL.tif")
centrality_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_centrality_15JUL.tif")
centrality_t_spatial_raster <- terra::rast("rasters_15JUL\\t_centrality_15JUL.tif")
outdegree_t_spatial_raster <- terra::rast("rasters_15JUL\\t_outdegree_15JUL.tif")
outdegree_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_outdegree_15JUL.tif")
indegree_t_spatial_raster <- terra::rast("rasters_15JUL\\t_indegree_15JUL.tif")
indegree_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_indegree_15JUL.tif")
closeness_t_spatial_raster <- terra::rast("rasters_15JUL\\t_closeness_15JUL.tif")
closeness_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_closeness_15JUL.tif")
tl_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_tl_15JUL.tif")
tl_t_spatial_raster <- terra::rast("rasters_15JUL\\t_tl_15JUL.tif")

#Difference metrics

?diffeR::differenceMR

diff_ivi <- diffeR::differenceMR(ivi_nt_spatial_raster, ivi_t_spatial_raster, eval = "multiple", percent = TRUE)
diff_centrality <- diffeR::differenceMR(centrality_nt_spatial_raster, centrality_t_spatial_raster, eval = "original")
diff_outdegee <- diffeR::differenceMR(outdegee_nt_spatial_raster, outdegee_t_spatial_raster, eval = "original")
diff_indegree <- diffeR::differenceMR(indegree_nt_spatial_raster, indegree_t_spatial_raster, eval = "original")
diff_closeness <- diffeR::differenceMR(closeness_nt_spatial_raster, closeness_t_spatial_raster, eval = "original")
diff_tl <- diffeR::differenceMR(tl_nt_spatial_raster, tl_t_spatial_raster, eval = "original")


#Save
saveRDS(diff_ivi, "diff_ivi.rds")
saveRDS(diff_centrality, "diff_centrality.rds")
saveRDS(diff_outdegee, "diff_outdegee.rds")
saveRDS(diff_indegree, "diff_indegree.rds")
saveRDS(diff_closeness, "diff_closeness.rds")
saveRDS(diff_tl, "diff_tl.rds")

#Load
diff_ivi <- readRDS("diffeR_results/diff_ivi.rds")
diff_centrality <- readRDS("diffeR_results/diff_centrality.rds")
diff_outdegee <- readRDS("diffeR_results/diff_outdegee.rds")
diff_indegree <- readRDS("diffeR_results/diff_indegree.rds")
diff_closeness <- readRDS("diffeR_results/diff_closeness.rds")
diff_tl <- readRDS("diffeR_results/diff_tl.rds")

#Create a data frame with all the results
diff_results <- data.frame(t(diff_ivi), 
                           t(diff_centrality),
                           t(diff_outdegee),
                           t(diff_indegree),
                           t(diff_closeness),
                           t(diff_tl)
)

#Rename the table of results
names(diff_results) <- c("ivi", 
                         "centrality",  
                         "outdegee", 
                         "indegree", 
                         "closeness",  
                         "tl"
) 

#View the table
View(diff_results)

# Compare at multiple scales
diff_ivi_2 <- diffeR::differenceMR(ivi_nt_spatial_raster, ivi_t_spatial_raster, eval = "multiple", percent = TRUE)
diff_centrality_2 <- diffeR::differenceMR(centrality_nt_spatial_raster, centrality_t_spatial_raster, eval = "multiple", percent = TRUE)
diff_outdegee_2 <- diffeR::differenceMR(outdegee_nt_spatial_raster, outdegee_t_spatial_raster, eval = "multiple", percent = TRUE)
diff_indegree_2 <- diffeR::differenceMR(indegree_nt_spatial_raster, indegree_t_spatial_raster, eval = "multiple", percent = TRUE)
diff_closeness_2 <- diffeR::differenceMR(closeness_nt_spatial_raster, closeness_t_spatial_raster, eval = "multiple", percent = TRUE)
diff_tl_2 <- diffeR::differenceMR(tl_nt_spatial_raster, tl_t_spatial_raster, eval = "multiple", percent = TRUE)

#Save
saveRDS(diff_ivi_2, "diff_ivi_2.rds")
saveRDS(diff_centrality_2, "diff_centrality_2.rds")
saveRDS(diff_outdegee_2, "diff_outdegee_2.rds")
saveRDS(diff_indegree_2, "diff_indegree_2.rds")
saveRDS(diff_closeness_2, "diff_closeness_2.rds")
saveRDS(diff_tl_2, "diff_tl_2.rds")

?diffeR::overallAllocD

ivi_cross <- crosstabm(ivi_nt_spatial_raster, ivi_t_spatial_raster)
ivi_overallAllocD <- overallAllocD(ivi_cross)
saveRDS(ivi_overallAllocD, "ivi_overallAllocD.rds")
#
centrality_cross <- crosstabm(centrality_nt_spatial_raster, centrality_t_spatial_raster)
centrality_overallAllocD <- overallAllocD(centrality_cross)
saveRDS(centrality_overallAllocD, "centrality_overallAllocD.rds")
#
outdegree_cross <- crosstabm(outdegree_nt_spatial_raster, outdegree_t_spatial_raster)
outdegree_overallAllocD <- overallAllocD(outdegree_cross)
saveRDS(outdegree_overallAllocD, "outdegree_overallAllocD.rds")
#
indegree_cross <- crosstabm(indegree_nt_spatial_raster, indegree_t_spatial_raster)
indegree_overallAllocD <- overallAllocD(indegree_cross)
saveRDS(indegree_overallAllocD, "indegree_overallAllocD.rds")
#
closeness_cross <- crosstabm(closeness_nt_spatial_raster, closeness_t_spatial_raster)
closeness_overallAllocD <- overallAllocD(closeness_cross)
saveRDS(closeness_overallAllocD, "closeness_overallAllocD.rds")
#
tl_cross <- crosstabm(tl_nt_spatial_raster, tl_t_spatial_raster)
tl_overallAllocD <- overallAllocD(tl_cross)
saveRDS(tl_overallAllocD, "tl_overallAllocD.rds")


#Load
ivi_overallAllocD <- readRDS("diffeR_results/ivi_overallAllocD.rds")
centrality_overallAllocD <- readRDS("diffeR_results/centrality_overallAllocD.rds")
outdegree_overallAllocD <- readRDS("diffeR_results/outdegree_overallAllocD.rds")
indegree_overallAllocD <- readRDS("diffeR_results/indegree_overallAllocD.rds")
closeness_overallAllocD <- readRDS("diffeR_results/closeness_overallAllocD.rds")
tl_overallAllocD <- readRDS("diffeR_results/tl_overallAllocD.rds")


################################################################################
#                                   fuzzySim
################################################################################

#FMestre
#03-10-2024

#Load packages
library(fuzzySim)

#Load rasters
ivi_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_ivi_15JUL.tif")
ivi_t_spatial_raster <- terra::rast("rasters_15JUL\\t_ivi_15JUL.tif")
centrality_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_centrality_15JUL.tif")
centrality_t_spatial_raster <- terra::rast("rasters_15JUL\\t_centrality_15JUL.tif")
outdegree_t_spatial_raster <- terra::rast("rasters_15JUL\\t_outdegree_15JUL.tif")
outdegree_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_outdegree_15JUL.tif")
indegree_t_spatial_raster <- terra::rast("rasters_15JUL\\t_indegree_15JUL.tif")
indegree_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_indegree_15JUL.tif")
closeness_t_spatial_raster <- terra::rast("rasters_15JUL\\t_closeness_15JUL.tif")
closeness_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_closeness_15JUL.tif")
tl_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_tl_15JUL.tif")
tl_t_spatial_raster <- terra::rast("rasters_15JUL\\t_tl_15JUL.tif")


#reduce resolution to avoid memory issues
ivi_nt_spatial_raster_lower_res <- aggregate(ivi_nt_spatial_raster, fact = 2, fun = mean)
ivi_t_spatial_raster_lower_res <- aggregate(ivi_t_spatial_raster, fact = 2, fun = mean)
centrality_nt_spatial_raster_lower_res <- aggregate(centrality_nt_spatial_raster, fact = 2, fun = mean)
centrality_t_spatial_raster_lower_res <- aggregate(centrality_t_spatial_raster, fact = 2, fun = mean)
outdegree_nt_spatial_raster_lower_res <- aggregate(outdegree_nt_spatial_raster, fact = 2, fun = mean)
outdegree_t_spatial_raster_lower_res <- aggregate(outdegree_t_spatial_raster, fact = 2, fun = mean)
indegree_nt_spatial_raster_lower_res <- aggregate(indegree_nt_spatial_raster, fact = 2, fun = mean)
indegree_t_spatial_raster_lower_res <- aggregate(indegree_t_spatial_raster, fact = 2, fun = mean)
closeness_nt_spatial_raster_lower_res <- aggregate(closeness_nt_spatial_raster, fact = 2, fun = mean)
closeness_t_spatial_raster_lower_res <- aggregate(closeness_t_spatial_raster, fact = 2, fun = mean)
tl_nt_spatial_raster_lower_res <- aggregate(tl_nt_spatial_raster, fact = 2, fun = mean)
tl_t_spatial_raster_lower_res <- aggregate(tl_t_spatial_raster, fact = 2, fun = mean)

#Convert to vector
ivi_nt_spatial_vector <- terra::values(ivi_nt_spatial_raster_lower_res, mat=FALSE)
ivi_t_spatial_vector <- terra::values(ivi_t_spatial_raster_lower_res, mat=FALSE)
centrality_nt_spatial_vector <- terra::values(centrality_nt_spatial_raster_lower_res, mat=FALSE)
centrality_t_spatial_vector <- terra::values(centrality_t_spatial_raster_lower_res, mat=FALSE)
outdegree_nt_spatial_vector <- terra::values(outdegree_nt_spatial_raster_lower_res, mat=FALSE)
outdegree_t_spatial_vector <- terra::values(outdegree_t_spatial_raster_lower_res, mat=FALSE)
indegree_nt_spatial_vector <- terra::values(indegree_nt_spatial_raster_lower_res, mat=FALSE)
indegree_t_spatial_vector <- terra::values(indegree_t_spatial_raster_lower_res, mat=FALSE)
closeness_nt_spatial_vector <- terra::values(closeness_nt_spatial_raster_lower_res, mat=FALSE)
closeness_t_spatial_vector <- terra::values(closeness_t_spatial_raster_lower_res, mat=FALSE)
tl_nt_spatial_vector <- terra::values(tl_nt_spatial_raster_lower_res, mat=FALSE)
tl_t_spatial_vector <- terra::values(tl_t_spatial_raster_lower_res, mat=FALSE)



# Create a similarity matrix for continuous rasters
?fuzzySim::similarity

#Run comparison
?fuzzySim::modOverlap
ivi_fuzzy_compare <- fuzzySim::modOverlap(ivi_nt_spatial_vector/100, ivi_t_spatial_vector/100)
centrality_fuzzy_compare <- fuzzySim::modOverlap(centrality_nt_spatial_vector/1000, centrality_t_spatial_vector/1000)
outdegree_fuzzy_compare <- fuzzySim::modOverlap(outdegree_nt_spatial_vector/100, outdegree_t_spatial_vector/100)
indegree_fuzzy_compare <- fuzzySim::modOverlap(indegree_nt_spatial_vector/100, indegree_t_spatial_vector/100)
closeness_fuzzy_compare <- fuzzySim::modOverlap(closeness_nt_spatial_vector/100, closeness_t_spatial_vector/100)
tl_fuzzy_compare <- fuzzySim::modOverlap(tl_nt_spatial_vector/100, tl_t_spatial_vector/100)

#Save
saveRDS(ivi_fuzzy_compare, "fuzzy_ivi.rds")
saveRDS(centrality_fuzzy_compare, "fuzzy_centrality.rds")
saveRDS(outdegree_fuzzy_compare, "fuzzy_outdegee.rds")
saveRDS(indegree_fuzzy_compare, "fuzzy_indegree.rds")
saveRDS(closeness_fuzzy_compare, "fuzzy_closeness.rds")
saveRDS(tl_fuzzy_compare, "fuzzy_tl.rds")

#Load
fuzzy_centrality <- readRDS("t_nt_fuzzy_outputs/fuzzy_centrality.rds")
fuzzy_closeness <- readRDS("t_nt_fuzzy_outputs/fuzzy_closeness.rds")
fuzzy_indegree <- readRDS("t_nt_fuzzy_outputs/fuzzy_indegree.rds")
fuzzy_ivi <- readRDS("t_nt_fuzzy_outputs/fuzzy_ivi.rds")
fuzzy_outdegee <- readRDS("t_nt_fuzzy_outputs/fuzzy_outdegee.rds")
fuzzy_tl <- readRDS("t_nt_fuzzy_outputs/fuzzy_tl.rds")


################################################################################
#                                     RMSE
################################################################################

#Load rasters
ivi_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_ivi_15JUL.tif")
ivi_t_spatial_raster <- terra::rast("rasters_15JUL\\t_ivi_15JUL.tif")
centrality_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_centrality_15JUL.tif")
centrality_t_spatial_raster <- terra::rast("rasters_15JUL\\t_centrality_15JUL.tif")
outdegree_t_spatial_raster <- terra::rast("rasters_15JUL\\t_outdegree_15JUL.tif")
outdegree_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_outdegree_15JUL.tif")
indegree_t_spatial_raster <- terra::rast("rasters_15JUL\\t_indegree_15JUL.tif")
indegree_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_indegree_15JUL.tif")
closeness_t_spatial_raster <- terra::rast("rasters_15JUL\\t_closeness_15JUL.tif")
closeness_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_closeness_15JUL.tif")
tl_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_tl_15JUL.tif")
tl_t_spatial_raster <- terra::rast("rasters_15JUL\\t_tl_15JUL.tif")

#reduce resolution to avoid memory issues
ivi_nt_spatial_raster_lower_res <- aggregate(ivi_nt_spatial_raster, fact = 2, fun = mean)
ivi_t_spatial_raster_lower_res <- aggregate(ivi_t_spatial_raster, fact = 2, fun = mean)
centrality_nt_spatial_raster_lower_res <- aggregate(centrality_nt_spatial_raster, fact = 2, fun = mean)
centrality_t_spatial_raster_lower_res <- aggregate(centrality_t_spatial_raster, fact = 2, fun = mean)
outdegree_nt_spatial_raster_lower_res <- aggregate(outdegree_nt_spatial_raster, fact = 2, fun = mean)
outdegree_t_spatial_raster_lower_res <- aggregate(outdegree_t_spatial_raster, fact = 2, fun = mean)
indegree_nt_spatial_raster_lower_res <- aggregate(indegree_nt_spatial_raster, fact = 2, fun = mean)
indegree_t_spatial_raster_lower_res <- aggregate(indegree_t_spatial_raster, fact = 2, fun = mean)
closeness_nt_spatial_raster_lower_res <- aggregate(closeness_nt_spatial_raster, fact = 2, fun = mean)
closeness_t_spatial_raster_lower_res <- aggregate(closeness_t_spatial_raster, fact = 2, fun = mean)

tl_nt_spatial_raster_lower_res <- aggregate(tl_nt_spatial_raster, fact = 2, fun = mean)
tl_t_spatial_raster_lower_res <- aggregate(tl_t_spatial_raster, fact = 2, fun = mean)

#Save aggregated
terra::writeRaster(ivi_nt_spatial_raster_lower_res, "ivi_nt_spatial_raster_lower_res.tif")
terra::writeRaster(ivi_t_spatial_raster_lower_res, "ivi_t_spatial_raster_lower_res.tif")
terra::writeRaster(centrality_nt_spatial_raster_lower_res, "centrality_nt_spatial_raster_lower_res.tif")
terra::writeRaster(centrality_t_spatial_raster_lower_res, "centrality_t_spatial_raster_lower_res.tif")
terra::writeRaster(outdegree_nt_spatial_raster_lower_res, "outdegree_nt_spatial_raster_lower_res.tif")
terra::writeRaster(outdegree_t_spatial_raster_lower_res, "outdegree_t_spatial_raster_lower_res.tif")
terra::writeRaster(indegree_nt_spatial_raster_lower_res, "indegree_nt_spatial_raster_lower_res.tif")
terra::writeRaster(indegree_t_spatial_raster_lower_res, "indegree_t_spatial_raster_lower_res.tif")
terra::writeRaster(closeness_nt_spatial_raster_lower_res, "closeness_nt_spatial_raster_lower_res.tif")
terra::writeRaster(closeness_t_spatial_raster_lower_res, "closeness_t_spatial_raster_lower_res.tif")
terra::writeRaster(tl_nt_spatial_raster_lower_res, "tl_nt_spatial_raster_lower_res.tif")
terra::writeRaster(tl_t_spatial_raster_lower_res, "tl_t_spatial_raster_lower_res.tif")


#Load aggregated
ivi_nt_spatial_raster_lower_res <- terra::rast("ivi_nt_spatial_raster_lower_res.tif")
ivi_t_spatial_raster_lower_res <- terra::rast("ivi_t_spatial_raster_lower_res.tif")
centrality_nt_spatial_raster_lower_res <- terra::rast("centrality_nt_spatial_raster_lower_res.tif")
centrality_t_spatial_raster_lower_res <- terra::rast("centrality_t_spatial_raster_lower_res.tif")
outdegree_nt_spatial_raster_lower_res <- terra::rast("outdegree_nt_spatial_raster_lower_res.tif")
outdegree_t_spatial_raster_lower_res <- terra::rast("outdegree_t_spatial_raster_lower_res.tif")
indegree_nt_spatial_raster_lower_res <- terra::rast("indegree_nt_spatial_raster_lower_res.tif")
indegree_t_spatial_raster_lower_res <- terra::rast("indegree_t_spatial_raster_lower_res.tif")
closeness_nt_spatial_raster_lower_res <- terra::rast("closeness_nt_spatial_raster_lower_res.tif")
closeness_t_spatial_raster_lower_res <- terra::rast("closeness_t_spatial_raster_lower_res.tif")
tl_nt_spatial_raster_lower_res <- terra::rast("tl_nt_spatial_raster_lower_res.tif")
tl_t_spatial_raster_lower_res <- terra::rast("tl_t_spatial_raster_lower_res.tif")


# Function to compute RMSE
rmse <- function(map1, map2) {
  sqrt(mean((map1 - map2)^2, na.rm = TRUE))
}

# Calculate RMSE
ivi_rmse <- rmse(ivi_t_spatial_raster_lower_res, ivi_nt_spatial_raster_lower_res)
#print(ivi_rmse)
plot(ivi_rmse)

################################################################################
#                Fuzzy Jaccard Index (Continuous Jaccard)
################################################################################

library(terra)

'
# Get the minimum and maximum values cell by cell
min_vals_ivi <- pmin(values(ivi_nt_spatial_raster_lower_res), values(ivi_t_spatial_raster_lower_res), na.rm = TRUE)
min_vals_centrality <- pmin(values(centrality_nt_spatial_raster_lower_res), values(centrality_t_spatial_raster_lower_res), na.rm = TRUE)
min_vals_outdegree <- pmin(values(outdegree_nt_spatial_raster_lower_res), values(outdegree_t_spatial_raster_lower_res), na.rm = TRUE)
min_vals_indegree <- pmin(values(indegree_nt_spatial_raster_lower_res), values(indegree_t_spatial_raster_lower_res), na.rm = TRUE)
min_vals_closeness <- pmin(values(closeness_nt_spatial_raster_lower_res), values(closeness_t_spatial_raster_lower_res), na.rm = TRUE)
min_vals_tl <- pmin(values(tl_nt_spatial_raster_lower_res), values(tl_t_spatial_raster_lower_res), na.rm = TRUE)

#Save
saveRDS(min_vals_ivi, "min_vals_ivi.rds")
saveRDS(min_vals_centrality, "min_vals_centrality.rds")
saveRDS(min_vals_outdegree, "min_vals_outdegree.rds")
saveRDS(min_vals_indegree, "min_vals_indegree.rds")
saveRDS(min_vals_closeness, "min_vals_closeness.rds")
saveRDS(min_vals_tl, "min_vals_tl.rds")

#Load
min_vals_ivi <- readRDS("min_vals_ivi.rds")
min_vals_centrality <- readRDS("min_vals_centrality.rds")
min_vals_outdegree <- readRDS("min_vals_outdegree.rds")
min_vals_indegree <- readRDS("min_vals_indegree.rds")
min_vals_closeness <- readRDS("min_vals_closeness.rds")
min_vals_tl <-readRDS("min_vals_tl.rds")

# Calculate the Fuzzy Jaccard Index
fuzzy_jaccard_ivi <- sum(min_vals_ivi, na.rm = TRUE) / sum(max_vals_ivi, na.rm = TRUE)
#fuzzy_jaccard_ivi
'

################################################################################
#          Overlap Coefficient (Sørensen-Dice for Continuous Data)
################################################################################

#ranges from 0 to 1, where 0 indicates that the sets have no common elements 
#and 1 indicates that the sets are identical

#Load packages
library(terra)

min_vals_ivi <- pmin(values(ivi_nt_spatial_raster_lower_res), values(ivi_t_spatial_raster_lower_res), na.rm = TRUE)
min_vals_centrality <- pmin(values(centrality_nt_spatial_raster_lower_res), values(centrality_t_spatial_raster_lower_res), na.rm = TRUE)
min_vals_outdegree <- pmin(values(outdegree_nt_spatial_raster_lower_res), values(outdegree_t_spatial_raster_lower_res), na.rm = TRUE)
min_vals_indegree <- pmin(values(indegree_nt_spatial_raster_lower_res), values(indegree_t_spatial_raster_lower_res), na.rm = TRUE)
min_vals_closeness <- pmin(values(closeness_nt_spatial_raster_lower_res), values(closeness_t_spatial_raster_lower_res), na.rm = TRUE)
min_vals_tl <- pmin(values(tl_nt_spatial_raster_lower_res), values(tl_t_spatial_raster_lower_res), na.rm = TRUE)

#Save
saveRDS(min_vals_ivi, "min_vals_ivi.rds")
saveRDS(min_vals_centrality, "min_vals_centrality.rds")
saveRDS(min_vals_outdegree, "min_vals_outdegree.rds")
saveRDS(min_vals_indegree, "min_vals_indegree.rds")
saveRDS(min_vals_closeness, "min_vals_closeness.rds")
saveRDS(min_vals_tl, "min_vals_tl.rds")

#Load
min_vals_ivi <- readRD("min_vals_ivi.rds")
min_vals_centrality <- readRDS("min_vals_centrality.rds")
min_vals_outdegree <- readRDS("min_vals_outdegree.rds")
min_vals_indegree <- readRDS("min_vals_indegree.rds")
min_vals_closeness <- readRDS("min_vals_closeness.rds")
min_vals_tl <-readRDS("min_vals_tl.rds")

# Calculate the Overlap Coefficient
ivi_overlap_coefficient <- 2 * sum(min_vals_ivi, na.rm = TRUE) / (sum(values(ivi_nt_spatial_raster_lower_res), na.rm = TRUE) + sum(values(ivi_t_spatial_raster_lower_res), na.rm = TRUE))
#ivi_overlap_coefficient
centrality_overlap_coefficient <- 2 * sum(min_vals_centrality, na.rm = TRUE) / (sum(values(centrality_nt_spatial_raster_lower_res), na.rm = TRUE) + sum(values(centrality_t_spatial_raster_lower_res), na.rm = TRUE))
#centrality_overlap_coefficient
outdegree_overlap_coefficient <- 2 * sum(min_vals_outdegree, na.rm = TRUE) / (sum(values(outdegree_nt_spatial_raster_lower_res), na.rm = TRUE) + sum(values(outdegree_t_spatial_raster_lower_res), na.rm = TRUE))
#outdegree_overlap_coefficient
indegree_overlap_coefficient <- 2 * sum(min_vals_indegree, na.rm = TRUE) / (sum(values(indegree_nt_spatial_raster_lower_res), na.rm = TRUE) + sum(values(indegree_t_spatial_raster_lower_res), na.rm = TRUE))
#outdegree_overlap_coefficient
closeness_overlap_coefficient <- 2 * sum(min_vals_closeness, na.rm = TRUE) / (sum(values(closeness_nt_spatial_raster_lower_res), na.rm = TRUE) + sum(values(closeness_t_spatial_raster_lower_res), na.rm = TRUE))
#outdegree_overlap_coefficient
tl_overlap_coefficient <- 2 * sum(min_vals_tl, na.rm = TRUE) / (sum(values(tl_nt_spatial_raster_lower_res), na.rm = TRUE) + sum(values(tl_t_spatial_raster_lower_res), na.rm = TRUE))
#outdegree_overlap_coefficient

#Save
saveRDS(ivi_overlap_coefficient, "ivi_overlap_coefficient.rds")
saveRDS(centrality_overlap_coefficient, "centrality_overlap_coefficient.rds")
saveRDS(outdegree_overlap_coefficient, "outdegree_overlap_coefficient.rds")
saveRDS(indegree_overlap_coefficient, "indegree_overlap_coefficient.rds")
saveRDS(closeness_overlap_coefficient, "closeness_overlap_coefficient.rds")
saveRDS(tl_overlap_coefficient, "tl_overlap_coefficient.rds")

#Load
ivi_overlap_coefficient <- readRDS("~/github/red_listed_networks/ivi_overlap_coefficient.rds")
tl_overlap_coefficient <- readRDS("~/github/red_listed_networks/tl_overlap_coefficient.rds")
closeness_overlap_coefficient <- readRDS("~/github/red_listed_networks/closeness_overlap_coefficient.rds")
indegree_overlap_coefficient <- readRDS("~/github/red_listed_networks/indegree_overlap_coefficient.rds")
outdegree_overlap_coefficient <- readRDS("~/github/red_listed_networks/outdegree_overlap_coefficient.rds")
centrality_overlap_coefficient <- readRDS("~/github/red_listed_networks/centrality_overlap_coefficient.rds")


#################################################################################

#reduce resolution to avoid memory issues
ivi_nt_spatial_raster_lower_res <- aggregate(ivi_nt_spatial_raster, fact = 2, fun = mean)
ivi_t_spatial_raster_lower_res <- aggregate(ivi_t_spatial_raster, fact = 2, fun = mean)
centrality_nt_spatial_raster_lower_res <- aggregate(centrality_nt_spatial_raster, fact = 2, fun = mean)
centrality_t_spatial_raster_lower_res <- aggregate(centrality_t_spatial_raster, fact = 2, fun = mean)
outdegree_nt_spatial_raster_lower_res <- aggregate(outdegree_nt_spatial_raster, fact = 2, fun = mean)
outdegree_t_spatial_raster_lower_res <- aggregate(outdegree_t_spatial_raster, fact = 2, fun = mean)
indegree_nt_spatial_raster_lower_res <- aggregate(indegree_nt_spatial_raster, fact = 2, fun = mean)
indegree_t_spatial_raster_lower_res <- aggregate(indegree_t_spatial_raster, fact = 2, fun = mean)
closeness_nt_spatial_raster_lower_res <- aggregate(closeness_nt_spatial_raster, fact = 2, fun = mean)
closeness_t_spatial_raster_lower_res <- aggregate(closeness_t_spatial_raster, fact = 2, fun = mean)
tl_nt_spatial_raster_lower_res <- aggregate(tl_nt_spatial_raster, fact = 2, fun = mean)
tl_t_spatial_raster_lower_res <- aggregate(tl_t_spatial_raster, fact = 2, fun = mean)

#Save aggregated
terra::writeRaster(ivi_nt_spatial_raster_lower_res, "ivi_nt_spatial_raster_lower_res.tif")
terra::writeRaster(ivi_t_spatial_raster_lower_res, "ivi_t_spatial_raster_lower_res.tif")
terra::writeRaster(centrality_nt_spatial_raster_lower_res, "centrality_nt_spatial_raster_lower_res.tif")
terra::writeRaster(centrality_t_spatial_raster_lower_res, "centrality_t_spatial_raster_lower_res.tif")
terra::writeRaster(outdegree_nt_spatial_raster_lower_res, "outdegree_nt_spatial_raster_lower_res.tif")
terra::writeRaster(outdegree_t_spatial_raster_lower_res, "outdegree_t_spatial_raster_lower_res.tif")
terra::writeRaster(indegree_nt_spatial_raster_lower_res, "indegree_nt_spatial_raster_lower_res.tif")
terra::writeRaster(indegree_t_spatial_raster_lower_res, "indegree_t_spatial_raster_lower_res.tif")
terra::writeRaster(closeness_nt_spatial_raster_lower_res, "closeness_nt_spatial_raster_lower_res.tif")
terra::writeRaster(closeness_t_spatial_raster_lower_res, "closeness_t_spatial_raster_lower_res.tif")
terra::writeRaster(tl_nt_spatial_raster_lower_res, "tl_nt_spatial_raster_lower_res.tif")
terra::writeRaster(tl_t_spatial_raster_lower_res, "tl_t_spatial_raster_lower_res.tif")

################################################################################
#                               Using SpatialPack
################################################################################

##### ivi #####

ivi_nt_spatial_raster_lower_res <- terra::rast("ivi_nt_spatial_raster_lower_res.tif")
ivi_t_spatial_raster_lower_res <- terra::rast("ivi_t_spatial_raster_lower_res.tif")
#
n_cells <- ncell(ivi_nt_spatial_raster_lower_res)

# Get the coordinates of the pixel centroids
ivi_nt_spatial_raster_lower_res_COORDS <- terra::xyFromCell(ivi_nt_spatial_raster_lower_res, 1:n_cells)
ivi_nt_spatial_raster_lower_res_VALUES <- terra::values(ivi_nt_spatial_raster_lower_res)
ivi_t_spatial_raster_lower_res_VALUES <- terra::values(ivi_t_spatial_raster_lower_res)

#Save
saveRDS(ivi_nt_spatial_raster_lower_res_COORDS, "ivi_nt_spatial_raster_lower_res_COORDS.rds")
saveRDS(ivi_nt_spatial_raster_lower_res_VALUES, "ivi_nt_spatial_raster_lower_res_VALUES.rds")
saveRDS(ivi_t_spatial_raster_lower_res_VALUES, "ivi_t_spatial_raster_lower_res_VALUES.rds")

#Load
#ivi_nt_spatial_raster_lower_res_COORDS <- readRDS("ivi_nt_spatial_raster_lower_res_COORDS.rds")
#ivi_nt_spatial_raster_lower_res_VALUES <- readRDS("ivi_nt_spatial_raster_lower_res_VALUES.rds")
#ivi_t_spatial_raster_lower_res_VALUES <- readRDS("ivi_t_spatial_raster_lower_res_VALUES.rds")

#Modified t-test
ivi_ttest <- modified.ttest(ivi_t_spatial_raster_lower_res_VALUES, 
                            ivi_nt_spatial_raster_lower_res_VALUES, 
                            ivi_nt_spatial_raster_lower_res_COORDS
)
#save
saveRDS(ivi_ttest, "ivi_ttest.rds")

##### tl #####

tl_nt_spatial_raster_lower_res <- terra::rast("tl_nt_spatial_raster_lower_res.tif")
tl_t_spatial_raster_lower_res <- terra::rast("tl_t_spatial_raster_lower_res.tif")
#
n_cells <- ncell(tl_nt_spatial_raster_lower_res)

# Get the coordinates of the pixel centroids
tl_nt_spatial_raster_lower_res_COORDS <- terra::xyFromCell(tl_nt_spatial_raster_lower_res, 1:n_cells)
tl_nt_spatial_raster_lower_res_VALUES <- terra::values(tl_nt_spatial_raster_lower_res)
tl_t_spatial_raster_lower_res_VALUES <- terra::values(tl_t_spatial_raster_lower_res)

#Save
saveRDS(tl_nt_spatial_raster_lower_res_COORDS, "tl_nt_spatial_raster_lower_res_COORDS.rds")
saveRDS(tl_nt_spatial_raster_lower_res_VALUES, "tl_nt_spatial_raster_lower_res_VALUES.rds")
saveRDS(tl_t_spatial_raster_lower_res_VALUES, "tl_t_spatial_raster_lower_res_VALUES.rds")

#Load
#tl_nt_spatial_raster_lower_res_COORDS <- readRDS("tl_nt_spatial_raster_lower_res_COORDS.rds")
#tl_nt_spatial_raster_lower_res_VALUES <- readRDS("tl_nt_spatial_raster_lower_res_VALUES.rds")
#tl_t_spatial_raster_lower_res_VALUES <- readRDS("tl_t_spatial_raster_lower_res_VALUES.rds")

#Modified t-test
tl_ttest <- modified.ttest(tl_t_spatial_raster_lower_res_VALUES, 
                           tl_nt_spatial_raster_lower_res_VALUES, 
                           tl_nt_spatial_raster_lower_res_COORDS
)
#save
saveRDS(tl_ttest, "tl_ttest.rds")

##### centrality #####

centrality_nt_spatial_raster_lower_res <- terra::rast("centrality_nt_spatial_raster_lower_res.tif")
centrality_t_spatial_raster_lower_res <- terra::rast("centrality_t_spatial_raster_lower_res.tif")
#
n_cells <- ncell(centrality_nt_spatial_raster_lower_res)

# Get the coordinates of the pixel centroids
centrality_nt_spatial_raster_lower_res_COORDS <- terra::xyFromCell(centrality_nt_spatial_raster_lower_res, 1:n_cells)
centrality_nt_spatial_raster_lower_res_VALUES <- terra::values(centrality_nt_spatial_raster_lower_res)
centrality_t_spatial_raster_lower_res_VALUES <- terra::values(centrality_t_spatial_raster_lower_res)

#Save
saveRDS(centrality_nt_spatial_raster_lower_res_COORDS, "centrality_nt_spatial_raster_lower_res_COORDS.rds")
saveRDS(centrality_nt_spatial_raster_lower_res_VALUES, "centrality_nt_spatial_raster_lower_res_VALUES.rds")
saveRDS(centrality_t_spatial_raster_lower_res_VALUES, "centrality_t_spatial_raster_lower_res_VALUES.rds")

#Load
#centrality_nt_spatial_raster_lower_res_COORDS <- readRDS("centrality_nt_spatial_raster_lower_res_COORDS.rds")
#centrality_nt_spatial_raster_lower_res_VALUES <- readRDS("centrality_nt_spatial_raster_lower_res_VALUES.rds")
#centrality_t_spatial_raster_lower_res_VALUES <- readRDS("centrality_t_spatial_raster_lower_res_VALUES.rds")

#Modified t-test
centrality_ttest <- modified.ttest(centrality_t_spatial_raster_lower_res_VALUES, 
                                   centrality_nt_spatial_raster_lower_res_VALUES, 
                                   centrality_nt_spatial_raster_lower_res_COORDS
)
#save
saveRDS(centrality_ttest, "centrality_ttest.rds")

##### closeness #####

closeness_nt_spatial_raster_lower_res <- terra::rast("closeness_nt_spatial_raster_lower_res.tif")
closeness_t_spatial_raster_lower_res <- terra::rast("closeness_t_spatial_raster_lower_res.tif")
#
n_cells <- ncell(closeness_nt_spatial_raster_lower_res)

# Get the coordinates of the pixel centroids
closeness_nt_spatial_raster_lower_res_COORDS <- terra::xyFromCell(closeness_nt_spatial_raster_lower_res, 1:n_cells)
closeness_nt_spatial_raster_lower_res_VALUES <- terra::values(closeness_nt_spatial_raster_lower_res)
closeness_t_spatial_raster_lower_res_VALUES <- terra::values(closeness_t_spatial_raster_lower_res)

#Save
saveRDS(closeness_nt_spatial_raster_lower_res_COORDS, "closeness_nt_spatial_raster_lower_res_COORDS.rds")
saveRDS(closeness_nt_spatial_raster_lower_res_VALUES, "closeness_nt_spatial_raster_lower_res_VALUES.rds")
saveRDS(closeness_t_spatial_raster_lower_res_VALUES, "closeness_t_spatial_raster_lower_res_VALUES.rds")

#Load
#closeness_nt_spatial_raster_lower_res_COORDS <- readRDS("closeness_nt_spatial_raster_lower_res_COORDS.rds")
#closeness_nt_spatial_raster_lower_res_VALUES <- readRDS("closeness_nt_spatial_raster_lower_res_VALUES.rds")
#closeness_t_spatial_raster_lower_res_VALUES <- readRDS("closeness_t_spatial_raster_lower_res_VALUES.rds")

#Modified t-test
closeness_ttest <- modified.ttest(closeness_t_spatial_raster_lower_res_VALUES, 
                                  closeness_nt_spatial_raster_lower_res_VALUES, 
                                  closeness_nt_spatial_raster_lower_res_COORDS
)
#save
saveRDS(closeness_ttest, "closeness_ttest.rds")

##### indegree #####

indegree_nt_spatial_raster_lower_res <- terra::rast("indegree_nt_spatial_raster_lower_res.tif")
indegree_t_spatial_raster_lower_res <- terra::rast("indegree_t_spatial_raster_lower_res.tif")
#
n_cells <- ncell(indegree_nt_spatial_raster_lower_res)

# Get the coordinates of the pixel centroids
indegree_nt_spatial_raster_lower_res_COORDS <- terra::xyFromCell(indegree_nt_spatial_raster_lower_res, 1:n_cells)
indegree_nt_spatial_raster_lower_res_VALUES <- terra::values(indegree_nt_spatial_raster_lower_res)
indegree_t_spatial_raster_lower_res_VALUES <- terra::values(indegree_t_spatial_raster_lower_res)

#Save
saveRDS(indegree_nt_spatial_raster_lower_res_COORDS, "indegree_nt_spatial_raster_lower_res_COORDS.rds")
saveRDS(indegree_nt_spatial_raster_lower_res_VALUES, "indegree_nt_spatial_raster_lower_res_VALUES.rds")
saveRDS(indegree_t_spatial_raster_lower_res_VALUES, "indegree_t_spatial_raster_lower_res_VALUES.rds")

#Load
#indegree_nt_spatial_raster_lower_res_COORDS <- readRDS("indegree_nt_spatial_raster_lower_res_COORDS.rds")
#indegree_nt_spatial_raster_lower_res_VALUES <- readRDS("indegree_nt_spatial_raster_lower_res_VALUES.rds")
#indegree_t_spatial_raster_lower_res_VALUES <- readRDS("indegree_t_spatial_raster_lower_res_VALUES.rds")

#Modified t-test
indegree_ttest <- modified.ttest(indegree_t_spatial_raster_lower_res_VALUES, 
                                 indegree_nt_spatial_raster_lower_res_VALUES, 
                                 indegree_nt_spatial_raster_lower_res_COORDS
)
#save
saveRDS(indegree_ttest, "indegree_ttest.rds")

##### outdegree #####

outdegree_nt_spatial_raster_lower_res <- terra::rast("outdegree_nt_spatial_raster_lower_res.tif")
outdegree_t_spatial_raster_lower_res <- terra::rast("outdegree_t_spatial_raster_lower_res.tif")
#
n_cells <- ncell(outdegree_nt_spatial_raster_lower_res)

# Get the coordinates of the pixel centroids
outdegree_nt_spatial_raster_lower_res_COORDS <- terra::xyFromCell(outdegree_nt_spatial_raster_lower_res, 1:n_cells)
outdegree_nt_spatial_raster_lower_res_VALUES <- terra::values(outdegree_nt_spatial_raster_lower_res)
outdegree_t_spatial_raster_lower_res_VALUES <- terra::values(outdegree_t_spatial_raster_lower_res)

#Save
saveRDS(outdegree_nt_spatial_raster_lower_res_COORDS, "outdegree_nt_spatial_raster_lower_res_COORDS.rds")
saveRDS(outdegree_nt_spatial_raster_lower_res_VALUES, "outdegree_nt_spatial_raster_lower_res_VALUES.rds")
saveRDS(outdegree_t_spatial_raster_lower_res_VALUES, "outdegree_t_spatial_raster_lower_res_VALUES.rds")

#Load
#outdegree_nt_spatial_raster_lower_res_COORDS <- readRDS("outdegree_nt_spatial_raster_lower_res_COORDS.rds")
#outdegree_nt_spatial_raster_lower_res_VALUES <- readRDS("outdegree_nt_spatial_raster_lower_res_VALUES.rds")
#outdegree_t_spatial_raster_lower_res_VALUES <- readRDS("outdegree_t_spatial_raster_lower_res_VALUES.rds")

#Modified t-test
outdegree_ttest <- modified.ttest(outdegree_t_spatial_raster_lower_res_VALUES, 
                                  outdegree_nt_spatial_raster_lower_res_VALUES, 
                                  outdegree_nt_spatial_raster_lower_res_COORDS
)
#save
saveRDS(outdegree_ttest, "outdegree_ttest.rds")


################################################################################
#                   Are the maps significantly different?    
################################################################################

#From: 
#https://sesync-ci.github.io/blog/raster-change-analysis.html

#indeg_diff_mean <- as.numeric(terra::global(indegree_diff, "mean", na.rm=TRUE))# mean
#indeg_diff_sd <- as.numeric(terra::global(indegree_diff, "sd", na.rm=TRUE))# sd
#indegree_diff_std <- (indegree_diff - indeg_diff_mean)/indeg_diff_sd # standardized image
#terra::writeRaster(indegree_diff_std, "rasters_15JUL\\indegree_diff_std.tif")
#terra::plot(indegree_diff_std)
#hist(indegree_diff_std, main="Standardized difference", xlab="Difference")

#(...)



################################################################################
#    Extract examples for the results & Discussion
################################################################################

network_list <- list.files("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_3")
pat_dataset <- "C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_3"

#Svaldbard

#C422
which(network_list == "C422.csv")
read.csv(paste0(pat_dataset, "\\C422.csv"))

#IVI
terra::global(ivi_diff, "mean", na.rm =TRUE)

#Closeness_diff
terra::global(closeness_diff, "mean", na.rm =TRUE)

#(...)



################################################################################
#    Combine
################################################################################

#indeg_diff
#outdeg_diff
#trophic_level_diff
#closeness_diff
#centrality_diff
#ivi_diff

# Define minimum and maximum values (adjust based on your data)
min_val_indeg_diff <- as.numeric(terra::global(indeg_diff, fun = "min", na.rm = TRUE))
max_val_indeg_diff <- as.numeric(terra::global(indeg_diff, fun = "max", na.rm = TRUE))
# Normalize values (0 to 1)
indeg_diff_normalized_rast <- (indeg_diff - min_val_indeg_diff) / (max_val_indeg_diff - min_val_indeg_diff)
#plot(indeg_diff_normalized_rast)

#Define minimum and maximum values (adjust based on your data)
min_val_outdeg_diff <- as.numeric(terra::global(outdeg_diff, fun = "min", na.rm = TRUE))
max_val_outdeg_diff <- as.numeric(terra::global(outdeg_diff, fun = "max", na.rm = TRUE))
# Normalize values (0 to 1)
outdeg_diff_normalized_rast <- (outdeg_diff - min_val_outdeg_diff) / (max_val_outdeg_diff - min_val_outdeg_diff)
#plot(outdeg_diff_normalized_rast)

#Define minimum and maximum values (adjust based on your data)
min_val_trophic_level_diff <- as.numeric(terra::global(trophic_level_diff, fun = "min", na.rm = TRUE))
max_val_trophic_level_diff <- as.numeric(terra::global(trophic_level_diff, fun = "max", na.rm = TRUE))
# Normalize values (0 to 1)
trophic_level_diff_normalized_rast <- (trophic_level_diff - min_val_trophic_level_diff) / (max_val_trophic_level_diff - min_val_trophic_level_diff)
#plot(trophic_level_diff_normalized_rast)

#Define minimum and maximum values (adjust based on your data)
min_val_closeness_diff <- as.numeric(terra::global(closeness_diff, fun = "min", na.rm = TRUE))
max_val_closeness_diff <- as.numeric(terra::global(closeness_diff, fun = "max", na.rm = TRUE))
# Normalize values (0 to 1)
closeness_diff_normalized_rast <- (closeness_diff - min_val_closeness_diff) / (max_val_closeness_diff - min_val_closeness_diff)
#plot(closeness_diff_normalized_rast)

#Define minimum and maximum values (adjust based on your data)
min_val_centrality_diff <- as.numeric(terra::global(centrality_diff, fun = "min", na.rm = TRUE))
max_val_centrality_diff <- as.numeric(terra::global(centrality_diff, fun = "max", na.rm = TRUE))
# Normalize values (0 to 1)
centrality_diff_normalized_rast <- (centrality_diff - min_val_centrality_diff) / (max_val_centrality_diff - min_val_centrality_diff)
#plot(centrality_diff_normalized_rast)

#Define minimum and maximum values (adjust based on your data)
min_val_ivi_diff <- as.numeric(terra::global(ivi_diff, fun = "min", na.rm = TRUE))
max_val_ivi_diff <- as.numeric(terra::global(ivi_diff, fun = "max", na.rm = TRUE))
# Normalize values (0 to 1)
ivi_diff_normalized_rast <- (ivi_diff - min_val_ivi_diff) / (max_val_ivi_diff - min_val_ivi_diff)
#plot(ivi_diff_normalized_rast)

all_diffs <- c(indeg_diff_normalized_rast,
               outdeg_diff_normalized_rast,
               trophic_level_diff_normalized_rast,
               closeness_diff_normalized_rast,
               centrality_diff_normalized_rast,
               ivi_diff_normalized_rast
)

mean_all_diffs <- terra::app(all_diffs, mean)
plot(mean_all_diffs)
rasterVis::levelplot(mean_all_diffs, par.settings=BuRdTheme())

#Save results
#terra::writeRaster(mean_all_diffs, "mean_all_diffs.tif")
