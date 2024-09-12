################################################################################
#                               Comparing rasters
################################################################################

#FMestre
#26-07-2024

#Load package
library(terra)
library(SSIMmap)

#################### LOAD RASTERS #################### 

ivi_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_ivi_15JUL.tif")
ivi_t_spatial_raster <- terra::rast("rasters_15JUL\\t_ivi_15JUL.tif")
centrality_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_centrality_15JUL.tif")
centrality_t_spatial_raster <- terra::rast("rasters_15JUL\\t_centrality_15JUL.tif")
indegree_t_spatial_raster <- terra::rast("rasters_15JUL\\t_indegree_15JUL.tif")
indegree_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_indegree_15JUL.tif")
outdegree_t_spatial_raster <- terra::rast("rasters_15JUL\\t_outdegree_15JUL.tif")
outdegree_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_outdegree_15JUL.tif")
closeness_t_spatial_raster <- terra::rast("rasters_15JUL\\t_closeness_15JUL.tif")
closeness_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_closeness_15JUL.tif")
tl_nt_spatial_raster <- terra::rast("rasters_15JUL\\nt_tl_15JUL.tif")
tl_t_spatial_raster <- terra::rast("rasters_15JUL\\t_tl_15JUL.tif")

#################################################

ivi_ssi_compare <- SSIMmap::ssim_raster(ivi_nt_spatial_raster, ivi_t_spatial_raster, global = FALSE)
centrality_ssi_compare <- SSIMmap::ssim_raster(centrality_nt_spatial_raster, centrality_t_spatial_raster, global = FALSE)
indegree_ssi_compare <- SSIMmap::ssim_raster(indegree_nt_spatial_raster, indegree_t_spatial_raster, global = FALSE)
outdegree_ssi_compare <- SSIMmap::ssim_raster(outdegree_nt_spatial_raster, outdegree_t_spatial_raster, global = FALSE)
closeness_ssi_compare <- SSIMmap::ssim_raster(closeness_nt_spatial_raster, closeness_t_spatial_raster, global = FALSE)
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
terra::writeRaster(indegree_ssi_compare[[1]], filename = "indegree_ssi_compare_SSIM.tif")
terra::writeRaster(indegree_ssi_compare[[2]], filename = "indegree_ssi_compare_SIM.tif")
terra::writeRaster(indegree_ssi_compare[[3]], filename = "indegree_ssi_compare_SIV.tif")
terra::writeRaster(indegree_ssi_compare[[4]], filename = "indegree_ssi_compare_SIP.tif")
#
terra::writeRaster(outdegree_ssi_compare[[1]], filename = "outdegree_ssi_compare_SSIM.tif")
terra::writeRaster(outdegree_ssi_compare[[2]], filename = "outdegree_ssi_compare_SIM.tif")
terra::writeRaster(outdegree_ssi_compare[[3]], filename = "outdegree_ssi_compare_SIV.tif")
terra::writeRaster(outdegree_ssi_compare[[4]], filename = "outdegree_ssi_compare_SIP.tif")
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

#Load packages
library(ggplot2)
library(RColorBrewer)

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
indegree_ssi_compare_SSIM <- terra::rast("D:\\output_red_listed\\indegree_ssi_compare_SSIM.tif")
indegree_ssi_compare_SIM <- terra::rast("D:\\output_red_listed\\indegree_ssi_compare_SIM.tif")
indegree_ssi_compare_SIV <- terra::rast("D:\\output_red_listed\\indegree_ssi_compare_SIV.tif")
indegree_ssi_compare_SIP <- terra::rast("D:\\output_red_listed\\indegree_ssi_compare_SIP.tif")
#
outdegree_ssi_compare_SSIM <- terra::rast("D:\\output_red_listed\\outdegree_ssi_compare_SSIM.tif")
outdegree_ssi_compare_SIM <- terra::rast("D:\\output_red_listed\\outdegree_ssi_compare_SIM.tif")
outdegree_ssi_compare_SIV <- terra::rast("D:\\output_red_listed\\outdegree_ssi_compare_SIV.tif")
outdegree_ssi_compare_SIP <- terra::rast("D:\\output_red_listed\\outdegree_ssi_compare_SIP.tif")
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


