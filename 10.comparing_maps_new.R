################################################################################
#                               Comparing rasters
################################################################################

#FMestre
#26-07-2024

#Load package
library(SpatialPack)
library(terra)

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

pearson_corr <- c(cor_tl$correlation[1,2],
cor_indegree$correlation[1,2],
cor_outdegree$correlation[1,2],
cor_centrality$correlation[1,2],
cor_closeness$correlation[1,2],
cor_ivi$correlation[1,2]
)

indexes_compared <- c(
  "trophic level",
  "in degree",
  "out degree",
  "centrality",
  "closeness",
  "IVI"
  )


pearson_df <- data.frame(indexes_compared, pearson_corr)
names(pearson_df) <- c("indexes", "pearson correlation")

base::saveRDS(pearson_df, "pearson_df.RDS")
