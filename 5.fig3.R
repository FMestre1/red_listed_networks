################################################################################
#                          Compute and  plot differences           
################################################################################

#FMestre
#23-10-2023

library(terra)

indeg_diff <- t_indegree - nt_indegree
terra::writeRaster(indeg_diff, filename = "rasters_21OUT\\indeg_diff.tif")
#plot(indeg_diff)
#
outdeg_diff <- t_outdegree - nt_outdegree
terra::writeRaster(outdeg_diff, filename = "rasters_21OUT\\outdeg_diff.tif")
#plot(outdeg_diff)
#
trophic_level_diff <- t_tl - nt_tl
terra::writeRaster(trophic_level_diff, filename = "rasters_21OUT\\trophic_level_diff.tif")
#plot(trophic_level_diff)
#
closeness_diff <- t_closeness - nt_closeness
terra::writeRaster(closeness_diff, filename = "rasters_21OUT\\closeness_diff.tif")
#plot(closeness_diff)
#
centrality_diff <- t_centrality - nt_centrality
terra::writeRaster(centrality_diff, filename = "rasters_21OUT\\centrality_diff.tif")
#plot(centrality_diff)
#
ivi_diff <- t_ivi - nt_ivi
terra::writeRaster(ivi_diff, filename = "rasters_21OUT\\ivi_diff.tif")
#plot(ivi_diff)

