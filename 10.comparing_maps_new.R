################################################################################
#                               Comparing rasters
################################################################################

#FMestre
#26-07-2024

#Load package
library(terra)
library(SSIMmap)
library(textGrob)
library(gridExtra)
library(grid)
library(RColorBrewer)
#library(spatialEco)
library(SSIMmap)

#################### LOAD RASTERS #################### 

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

#library(terra)

# Get the minimum and maximum values cell by cell
#min_vals_ivi <- pmin(values(ivi_nt_spatial_raster_lower_res), values(ivi_t_spatial_raster_lower_res), na.rm = TRUE)
#max_vals_ivi <- pmax(values(ivi_nt_spatial_raster_lower_res), values(ivi_t_spatial_raster_lower_res), na.rm = TRUE)

# Calculate the Fuzzy Jaccard Index
#fuzzy_jaccard_ivi <- sum(min_vals_ivi, na.rm = TRUE) / sum(max_vals_ivi, na.rm = TRUE)
#fuzzy_jaccard_ivi

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
centralty_overlap_coefficient <- readRDS("~/github/red_listed_networks/centralty_overlap_coefficient.rds")
outdegree_overlap_coefficient <- readRDS("~/github/red_listed_networks/outdegree_overlap_coefficient.rds")
indegree_overlap_coefficient <- readRDS("~/github/red_listed_networks/indegree_overlap_coefficient.rds")
closeness_overlap_coefficient <- readRDS("~/github/red_listed_networks/closeness_overlap_coefficient.rds")
tl_overlap_coefficient <- readRDS("~/github/red_listed_networks/tl_overlap_coefficient.rds")


