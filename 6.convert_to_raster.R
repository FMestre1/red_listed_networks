
################################################################################
#                             Convert to rasters
################################################################################

template_raster <- terra::rast("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\old_results\\outdeg_diff.tif")

#IVI
nt_ivi <- terra::rasterize(x = ivi_nt_spatial, 
                           y = template_raster, 
                           field = "ivi")

#terra::writeRaster(nt_ivi, "rasters_15JUL\\nt_ivi_15JUL.tif")

t_ivi <- terra::rasterize(x = ivi_t_spatial, 
                          y = template_raster, 
                          field = "ivi")

#terra::writeRaster(t_ivi, "rasters_15JUL\\t_ivi_15JUL.tif")


#Centrality
t_centrality <- terra::rasterize(x = centrality_t_spatial, 
                                 y = template_raster, 
                                 field = "centrality")

#terra::writeRaster(t_centrality, "rasters_15JUL\\t_centrality_15JUL.tif")

nt_centrality <- terra::rasterize(x = centrality_nt_spatial, 
                                  y = template_raster, 
                                  field = "centrality")

#terra::writeRaster(nt_centrality, "rasters_15JUL\\nt_centrality_15JUL.tif")

#Indegree
t_indegree <- terra::rasterize(x = indegree_t_spatial, 
                               y = template_raster, 
                               field = "indegree")

#terra::writeRaster(t_indegree, "rasters_15JUL\\t_indegree_15JUL.tif")

nt_indegree <- terra::rasterize(x = indegree_nt_spatial, 
                                y = template_raster, 
                                field = "indegree")

#terra::writeRaster(nt_indegree, "rasters_15JUL\\nt_indegree_15JUL.tif")

#Outdegree
nt_outdegree <- terra::rasterize(x = outdegree_nt_spatial, 
                                 y = template_raster, 
                                 field = "outdegree")

#terra::writeRaster(nt_outdegree, "rasters_15JUL\\nt_outdegree_15JUL.tif")

t_outdegree <- terra::rasterize(x = outdegree_t_spatial, 
                                y = template_raster, 
                                field = "outdegree")

#terra::writeRaster(t_outdegree, "rasters_15JUL\\t_outdegree_15JUL.tif")

#Closeness
t_closeness <- terra::rasterize(x = closeness_t_spatial, 
                                y = template_raster, 
                                field = "closeness")

#terra::writeRaster(t_closeness, "rasters_15JUL\\t_closeness_15JUL.tif")

nt_closeness <- terra::rasterize(x = closeness_nt_spatial, 
                                 y = template_raster, 
                                 field = "closeness")

#terra::writeRaster(nt_closeness, "rasters_15JUL\\nt_closeness_15JUL.tif")

#TL
t_tl <- terra::rasterize(x = tl_t_spatial, 
                         y = template_raster, 
                         field = "trophic_level")

#terra::writeRaster(t_tl, "rasters_15JUL\\t_tl_15JUL.tif")

nt_tl <- terra::rasterize(x = tl_nt_spatial, 
                          y = template_raster, 
                          field = "trophic_level")

#terra::writeRaster(nt_tl, "rasters_15JUL\\nt_tl_15JUL.tif")

#Proportion
proportion_r <- terra::rasterize(x = proportion_spatial, 
                                 y = template_raster, 
                                 field = "proportion")

#terra::writeRaster(proportion_r, "rasters_15JUL\\proportion_r__15JUL.tif")

#Species richness
species_rich_r <- terra::rasterize(x = sp_richness, 
                                 y = template_raster, 
                                 field = "spe_igraph")

#terra::writeRaster(species_rich_r, "rasters_15JUL\\species_rich_r_15JUL.tif")
