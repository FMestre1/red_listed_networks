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
tl_t_spatial_raster <- terra::rast("rasters_15JUL\\t_tl_15JUL.tif")
tl_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_tl_15JUL.tif")
tl_t_spatial_raster <- terra::rast("rasters_15JUL\\t_tl_15JUL.tif")
tl_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_tl_15JUL.tif")
tl_t_spatial_raster <- terra::rast("rasters_15JUL\\t_tl_15JUL.tif")
tl_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_tl_15JUL.tif")
tl_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_tl_15JUL.tif")
tl_t_spatial_raster <- terra::rast("rasters_15JUL\\t_tl_15JUL.tif")

#################################################

ivi_ssi_compare <- SSIMmap::ssim_raster(ivi_nt_spatial_raster, ivi_t_spatial_raster, global = FALSE)
centrality_ssi_compare <- SSIMmap::ssim_raster(centrality_nt_spatial_raster, centrality_t_spatial_raster, global = FALSE)
tl_ssi_compare <- SSIMmap::ssim_raster(tl_nt_spatial_raster, tl_t_spatial_raster, global = FALSE)
tl_ssi_compare <- SSIMmap::ssim_raster(tl_nt_spatial_raster, tl_t_spatial_raster, global = FALSE)
tl_ssi_compare <- SSIMmap::ssim_raster(tl_nt_spatial_raster, tl_t_spatial_raster, global = FALSE)
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
terra::writeRaster(tl_ssi_compare[[1]], filename = "tl_ssi_compare_SSIM.tif")
terra::writeRaster(tl_ssi_compare[[2]], filename = "tl_ssi_compare_SIM.tif")
terra::writeRaster(tl_ssi_compare[[3]], filename = "tl_ssi_compare_SIV.tif")
terra::writeRaster(tl_ssi_compare[[4]], filename = "tl_ssi_compare_SIP.tif")
#
terra::writeRaster(tl_ssi_compare[[1]], filename = "tl_ssi_compare_SSIM.tif")
terra::writeRaster(tl_ssi_compare[[2]], filename = "tl_ssi_compare_SIM.tif")
terra::writeRaster(tl_ssi_compare[[3]], filename = "tl_ssi_compare_SIV.tif")
terra::writeRaster(tl_ssi_compare[[4]], filename = "tl_ssi_compare_SIP.tif")
#
terra::writeRaster(tl_ssi_compare[[1]], filename = "tl_ssi_compare_SSIM.tif")
terra::writeRaster(tl_ssi_compare[[2]], filename = "tl_ssi_compare_SIM.tif")
terra::writeRaster(tl_ssi_compare[[3]], filename = "tl_ssi_compare_SIV.tif")
terra::writeRaster(tl_ssi_compare[[4]], filename = "tl_ssi_compare_SIP.tif")
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
tl_ssi_compare_SSIM <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SSIM.tif")
tl_ssi_compare_SIM <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIM.tif")
tl_ssi_compare_SIV <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIV.tif")
tl_ssi_compare_SIP <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIP.tif")
#
tl_ssi_compare_SSIM <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SSIM.tif")
tl_ssi_compare_SIM <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIM.tif")
tl_ssi_compare_SIV <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIV.tif")
tl_ssi_compare_SIP <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIP.tif")
#
tl_ssi_compare_SSIM <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SSIM.tif")
tl_ssi_compare_SIM <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIM.tif")
tl_ssi_compare_SIV <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIV.tif")
tl_ssi_compare_SIP <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIP.tif")
#
tl_ssi_compare_SSIM <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SSIM.tif")
tl_ssi_compare_SIM <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIM.tif")
tl_ssi_compare_SIV <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIV.tif")
tl_ssi_compare_SIP <- terra::rast("D:\\output_red_listed\\tl_ssi_compare_SIP.tif")




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
                        top=grid::textGrob("tl"))


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
                        top=grid::textGrob("tl"))


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
                        top=grid::textGrob("tl"))

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
                        top=grid::textGrob("TL"))
