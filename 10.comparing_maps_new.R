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

#################################################

ivi_ssi_compare <- SSIMmap::ssim_raster(ivi_nt_spatial_raster, ivi_t_spatial_raster, global = FALSE)
centrality_ssi_compare <- SSIMmap::ssim_raster(centrality_nt_spatial_raster, centrality_t_spatial_raster, global = FALSE)
closeness_ssi_compare <- SSIMmap::ssim_raster(closeness_nt_spatial_raster, closeness_t_spatial_raster, global = FALSE)
indegree_ssi_compare <- SSIMmap::ssim_raster(indegree_nt_spatial_raster, indegree_t_spatial_raster, global = FALSE)
outdegree_ssi_compare <- SSIMmap::ssim_raster(outdegree_nt_spatial_raster, outdegree_t_spatial_raster, global = FALSE)
tl_ssi_compare <- SSIMmap::ssim_raster(tl_nt_spatial_raster, tl_t_spatial_raster, global = FALSE)

#Save...

terra::writeRaster(ivi_ssi_compare[[1]], filename = "ivi_ssi_compare_SSIM.tif")
terra::writeRaster(ivi_ssi_compare[[2]], filename = "ivi_ssi_compare_SIM.tif")
terra::writeRaster(ivi_ssi_compare[[3]], filename = "ivi_ssi_compare_SIV.tif")
terra::writeRaster(ivi_ssi_compare[[4]], filename = "ivi_ssi_compare_SIP.tif")
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
#
terra::writeRaster(tl_ssi_compare[[1]], filename = "tl_ssi_compare_SSIM.tif")
terra::writeRaster(tl_ssi_compare[[2]], filename = "tl_ssi_compare_SIM.tif")
terra::writeRaster(tl_ssi_compare[[3]], filename = "tl_ssi_compare_SIV.tif")
terra::writeRaster(tl_ssi_compare[[4]], filename = "tl_ssi_compare_SIP.tif")

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

