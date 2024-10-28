################################################################################
#       Is there an effect of species richness or proportion of threatened         
################################################################################

#FMestre
#28-10-2024

#1. ##### Upload rasters #####

#Get tthe difference maps
indegree_diff_z_score <- terra::rast("indegree_diff_z_score.tif")
outdegree_diff_z_score <- terra::rast("outdegree_diff_z_score.tif")
t_level_diff_z_score <- terra::rast("t_level_diff_z_score.tif")
closeness_diff_z_score <- terra::rast("closeness_diff_z_score.tif")
centrality_diff_z_score <- terra::rast("centrality_diff_z_score.tif")
ivi_diff_z_score <- terra::rast("ivi_diff_z_score.tif")

#Get the proportion of threatened species
proportion_r <- terra::rast("rasters_15JUL\\proportion_r__15JUL.tif")

#Get the number of species
species_rich_r <- terra::rasterize(x = sp_richness, 
                                   y = indegree_diff_z_score,
                                   field = "spe_igraph")

terra::writeRaster(species_rich_r, "species_rich_r_28OCT24.tif")
species_rich_r <- terra::rast("species_rich_r_28OCT24.tif")


#2. ##### Plot #####
plot(species_rich_r, indegree_diff_z_score)
plot(species_rich_r, outdegree_diff_z_score)
plot(species_rich_r, t_level_diff_z_score)
plot(species_rich_r, closeness_diff_z_score)
plot(species_rich_r, centrality_diff_z_score)
plot(species_rich_r, ivi_diff_z_score)
#
plot(proportion_r, indegree_diff_z_score)
plot(proportion_r, outdegree_diff_z_score)
plot(proportion_r, t_level_diff_z_score)
plot(proportion_r, closeness_diff_z_score)
plot(proportion_r, centrality_diff_z_score)
plot(proportion_r, ivi_diff_z_score)

