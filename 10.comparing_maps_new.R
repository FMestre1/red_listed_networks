################################################################################
#                               Comparing rasters
################################################################################

#FMestre
#26-07-2024

#Load package
library(spatialEco)
library(terra)

#terraOptions(memfrac = 0.5)
?terraOptions
terraOptions()

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

# process in parallel
library(doParallel) 
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

