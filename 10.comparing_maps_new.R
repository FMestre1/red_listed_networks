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
tl_nt_spatial_raster <- terra::rast("tl_nt_spatial_raster_z_raster.tif")
tl_t_spatial_raster <- terra::rast("tl_t_spatial_raster_z_raster.tif")

ivi_nt_spatial_raster <- terra::rast("ivi_nt_spatial_raster_z_raster.tif")
ivi_t_spatial_raster <- terra::rast("ivi_t_spatial_raster_z_raster.tif")

centrality_nt_spatial_raster <- terra::rast("centrality_nt_spatial_raster_z_raster.tif")
centrality_t_spatial_raster <- terra::rast("centrality_t_spatial_raster_z_raster.tif")

closeness_nt_spatial_raster <- terra::rast("closeness_nt_spatial_raster_z_raster.tif")
closeness_t_spatial_raster <- terra::rast("closeness_t_spatial_raster_z_raster.tif")

indegree_nt_spatial_raster <- terra::rast("indegree_nt_spatial_raster_z_raster.tif")
indegree_t_spatial_raster <- terra::rast("indegree_t_spatial_raster_z_raster.tif")

outdegree_nt_spatial_raster <- terra::rast("outdegree_nt_spatial_raster_z_raster.tif")
outdegree_t_spatial_raster <- terra::rast("outdegree_t_spatial_raster_z_raster.tif")


#tl_ssi_compare <- SSIMmap::ssim_raster(tl_nt_spatial_raster, tl_t_spatial_raster, global = FALSE)
ivi_ssi_compare <- SSIMmap::ssim_raster(ivi_nt_spatial_raster, ivi_t_spatial_raster, global = FALSE)
centrality_ssi_compare <- SSIMmap::ssim_raster(centrality_nt_spatial_raster, centrality_t_spatial_raster, global = FALSE)
closeness_ssi_compare <- SSIMmap::ssim_raster(closeness_nt_spatial_raster, closeness_t_spatial_raster, global = FALSE)
indegree_ssi_compare <- SSIMmap::ssim_raster(indegree_nt_spatial_raster, indegree_t_spatial_raster, global = FALSE)
outdegree_ssi_compare <- SSIMmap::ssim_raster(outdegree_nt_spatial_raster, outdegree_t_spatial_raster, global = FALSE)

#Save...
terra::writeRaster(ivi_ssi_compare[[1]], filename = "ivi_ssi_compare_SSIM.tif")
terra::writeRaster(ivi_ssi_compare[[2]], filename = "ivi_ssi_compare_SIM.tif")
terra::writeRaster(ivi_ssi_compare[[3]], filename = "ivi_ssi_compare_SIV.tif")
terra::writeRaster(ivi_ssi_compare[[4]], filename = "ivi_ssi_compare_SIP.tif")
#
#terra::writeRaster(tl_ssi_compare[[1]], filename = "tl_ssi_compare_SSIM.tif")
#terra::writeRaster(tl_ssi_compare[[2]], filename = "tl_ssi_compare_SIM.tif")
#terra::writeRaster(tl_ssi_compare[[3]], filename = "tl_ssi_compare_SIV.tif")
#terra::writeRaster(tl_ssi_compare[[4]], filename = "tl_ssi_compare_SIP.tif")
#
terra::writeRaster(centrality_ssi_compare[[1]], filename = "centrality_ssi_compare_SSIM.tif")
terra::writeRaster(centrality_ssi_compare[[2]], filename = "centrality_ssi_compare_SIM.tif")
terra::writeRaster(centrality_ssi_compare[[3]], filename = "centrality_ssi_compare_SIV.tif")
terra::writeRaster(centrality_ssi_compare[[4]], filename = "centrality_ssi_compare_SIP.tif")
#
terra::writeRaster(outdegree_ssi_compare[[1]], filename = "outdegree_ssi_compare_SSIM.tif")
terra::writeRaster(outdegree_ssi_compare[[2]], filename = "outdegree_ssi_compare_SIM.tif")
terra::writeRaster(outdegree_ssi_compare[[3]], filename = "outdegree_ssi_compare_SIV.tif")
terra::writeRaster(outdegree_ssi_compare[[4]], filename = "outdegree_ssi_compare_SIP.tif")
#
terra::writeRaster(indegree_ssi_compare[[1]], filename = "indegree_ssi_compare_SSIM.tif")
terra::writeRaster(indegree_ssi_compare[[2]], filename = "indegree_ssi_compare_SIM.tif")
terra::writeRaster(indegree_ssi_compare[[3]], filename = "indegree_ssi_compare_SIV.tif")
terra::writeRaster(indegree_ssi_compare[[4]], filename = "indegree_ssi_compare_SIP.tif")
#
terra::writeRaster(closeness_ssi_compare[[1]], filename = "closeness_ssi_compare_SSIM.tif")
terra::writeRaster(closeness_ssi_compare[[2]], filename = "closeness_ssi_compare_SIM.tif")
terra::writeRaster(closeness_ssi_compare[[3]], filename = "closeness_ssi_compare_SIV.tif")
terra::writeRaster(closeness_ssi_compare[[4]], filename = "closeness_ssi_compare_SIP.tif") 
################################################################################
#                                  Plotting
################################################################################

#https://cran.r-project.org/web/packages/SSIMmap/vignettes/Introduction_to_SSIMmap.html

ivi_ssi_compare_SSIM <- terra::rast("D:\\output_red_listed\\ivi_ssi_compare_SSIM.tif")
ivi_ssi_compare_SIM <- terra::rast("D:\\output_red_listed\\ivi_ssi_compare_SIM.tif")
ivi_ssi_compare_SIV <- terra::rast("D:\\output_red_listed\\ivi_ssi_compare_SIV.tif")
ivi_ssi_compare_SIP <- terra::rast("D:\\output_red_listed\\ivi_ssi_compare_SIP.tif")
#
centrality_ssi_compare_SSIM <- terra::rast("D:\\output_red_listed\\centrality_ssi_compare_SSIM.tif")
centrality_ssi_compare_SIM <- terra::rast("D:\\output_red_listed\\centrality_ssi_compare_SIM.tif")
centrality_ssi_compare_SIV <- terra::rast("D:\\output_red_listed\\centrality_ssi_compare_SIV.tif")
centrality_ssi_compare_SIP <- terra::rast("D:\\output_red_listed\\centrality_ssi_compare_SIP.tif")
#
outdegree_ssi_compare_SSIM <- terra::rast("D:\\output_red_listed\\outdegree_ssi_compare_SSIM.tif")
outdegree_ssi_compare_SIM <- terra::rast("D:\\output_red_listed\\outdegree_ssi_compare_SIM.tif")
outdegree_ssi_compare_SIV <- terra::rast("D:\\output_red_listed\\outdegree_ssi_compare_SIV.tif")
outdegree_ssi_compare_SIP <- terra::rast("D:\\output_red_listed\\outdegree_ssi_compare_SIP.tif")
#
indegree_ssi_compare_SSIM <- terra::rast("D:\\output_red_listed\\indegree_ssi_compare_SSIM.tif")
indegree_ssi_compare_SIM <- terra::rast("D:\\output_red_listed\\indegree_ssi_compare_SIM.tif")
indegree_ssi_compare_SIV <- terra::rast("D:\\output_red_listed\\indegree_ssi_compare_SIV.tif")
indegree_ssi_compare_SIP <- terra::rast("D:\\output_red_listed\\indegree_ssi_compare_SIP.tif")
#
closeness_ssi_compare_SSIM <- terra::rast("D:\\output_red_listed\\closeness_ssi_compare_SSIM.tif")
closeness_ssi_compare_SIM <- terra::rast("D:\\output_red_listed\\closeness_ssi_compare_SIM.tif")
closeness_ssi_compare_SIV <- terra::rast("D:\\output_red_listed\\closeness_ssi_compare_SIV.tif")
closeness_ssi_compare_SIP <- terra::rast("D:\\output_red_listed\\closeness_ssi_compare_SIP.tif")
#
tl_ssi_compare_SSIM <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SSIM.tif")
tl_ssi_compare_SIM <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIM.tif")
tl_ssi_compare_SIV <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIV.tif")
tl_ssi_compare_SIP <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIP.tif")

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

#########################

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

