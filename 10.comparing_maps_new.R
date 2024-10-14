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

################################################################################
#                               Using SpatialPack
################################################################################

#Load package
library(SpatialPack)

##### ivi #####

ivi_nt_spatial_raster_lower_res <- terra::rast("ivi_nt_spatial_raster_lower_res.tif")
ivi_t_spatial_raster_lower_res <- terra::rast("ivi_t_spatial_raster_lower_res.tif")
#
n_cells <- ncell(ivi_nt_spatial_raster_lower_res)

# Get the coordinates of the pixel centroids
ivi_nt_spatial_raster_lower_res_COORDS <- terra::xyFromCell(ivi_nt_spatial_raster_lower_res, 1:n_cells)
ivi_nt_spatial_raster_lower_res_VALUES <- terra::values(ivi_nt_spatial_raster_lower_res)
ivi_t_spatial_raster_lower_res_VALUES <- terra::values(ivi_t_spatial_raster_lower_res)

#Save
saveRDS(ivi_nt_spatial_raster_lower_res_COORDS, "ivi_nt_spatial_raster_lower_res_COORDS.rds")
saveRDS(ivi_nt_spatial_raster_lower_res_VALUES, "ivi_nt_spatial_raster_lower_res_VALUES.rds")
saveRDS(ivi_t_spatial_raster_lower_res_VALUES, "ivi_t_spatial_raster_lower_res_VALUES.rds")

#Load
#ivi_nt_spatial_raster_lower_res_COORDS <- readRDS("ivi_nt_spatial_raster_lower_res_COORDS.rds")
#ivi_nt_spatial_raster_lower_res_VALUES <- readRDS("ivi_nt_spatial_raster_lower_res_VALUES.rds")
#ivi_t_spatial_raster_lower_res_VALUES <- readRDS("ivi_t_spatial_raster_lower_res_VALUES.rds")

#Modified t-test
ivi_ttest <- modified.ttest(ivi_t_spatial_raster_lower_res_VALUES, 
                            ivi_nt_spatial_raster_lower_res_VALUES, 
                            ivi_nt_spatial_raster_lower_res_COORDS
                            )
#save
saveRDS(ivi_ttest, "ivi_ttest.rds")

##### tl #####

tl_nt_spatial_raster_lower_res <- terra::rast("tl_nt_spatial_raster_lower_res.tif")
tl_t_spatial_raster_lower_res <- terra::rast("tl_t_spatial_raster_lower_res.tif")
#
n_cells <- ncell(tl_nt_spatial_raster_lower_res)

# Get the coordinates of the pixel centroids
tl_nt_spatial_raster_lower_res_COORDS <- terra::xyFromCell(tl_nt_spatial_raster_lower_res, 1:n_cells)
tl_nt_spatial_raster_lower_res_VALUES <- terra::values(tl_nt_spatial_raster_lower_res)
tl_t_spatial_raster_lower_res_VALUES <- terra::values(tl_t_spatial_raster_lower_res)

#Save
saveRDS(tl_nt_spatial_raster_lower_res_COORDS, "tl_nt_spatial_raster_lower_res_COORDS.rds")
saveRDS(tl_nt_spatial_raster_lower_res_VALUES, "tl_nt_spatial_raster_lower_res_VALUES.rds")
saveRDS(tl_t_spatial_raster_lower_res_VALUES, "tl_t_spatial_raster_lower_res_VALUES.rds")

#Load
#tl_nt_spatial_raster_lower_res_COORDS <- readRDS("tl_nt_spatial_raster_lower_res_COORDS.rds")
#tl_nt_spatial_raster_lower_res_VALUES <- readRDS("tl_nt_spatial_raster_lower_res_VALUES.rds")
#tl_t_spatial_raster_lower_res_VALUES <- readRDS("tl_t_spatial_raster_lower_res_VALUES.rds")

#Modified t-test
tl_ttest <- modified.ttest(tl_t_spatial_raster_lower_res_VALUES, 
                           tl_nt_spatial_raster_lower_res_VALUES, 
                           tl_nt_spatial_raster_lower_res_COORDS
)
#save
saveRDS(tl_ttest, "tl_ttest.rds")

##### centrality #####

centrality_nt_spatial_raster_lower_res <- terra::rast("centrality_nt_spatial_raster_lower_res.tif")
centrality_t_spatial_raster_lower_res <- terra::rast("centrality_t_spatial_raster_lower_res.tif")
#
n_cells <- ncell(centrality_nt_spatial_raster_lower_res)

# Get the coordinates of the pixel centroids
centrality_nt_spatial_raster_lower_res_COORDS <- terra::xyFromCell(centrality_nt_spatial_raster_lower_res, 1:n_cells)
centrality_nt_spatial_raster_lower_res_VALUES <- terra::values(centrality_nt_spatial_raster_lower_res)
centrality_t_spatial_raster_lower_res_VALUES <- terra::values(centrality_t_spatial_raster_lower_res)

#Save
saveRDS(centrality_nt_spatial_raster_lower_res_COORDS, "centrality_nt_spatial_raster_lower_res_COORDS.rds")
saveRDS(centrality_nt_spatial_raster_lower_res_VALUES, "centrality_nt_spatial_raster_lower_res_VALUES.rds")
saveRDS(centrality_t_spatial_raster_lower_res_VALUES, "centrality_t_spatial_raster_lower_res_VALUES.rds")

#Load
#centrality_nt_spatial_raster_lower_res_COORDS <- readRDS("centrality_nt_spatial_raster_lower_res_COORDS.rds")
#centrality_nt_spatial_raster_lower_res_VALUES <- readRDS("centrality_nt_spatial_raster_lower_res_VALUES.rds")
#centrality_t_spatial_raster_lower_res_VALUES <- readRDS("centrality_t_spatial_raster_lower_res_VALUES.rds")

#Modified t-test
centrality_ttest <- modified.ttest(centrality_t_spatial_raster_lower_res_VALUES, 
                                   centrality_nt_spatial_raster_lower_res_VALUES, 
                                   centrality_nt_spatial_raster_lower_res_COORDS
)
#save
saveRDS(centrality_ttest, "centrality_ttest.rds")

##### closeness #####

closeness_nt_spatial_raster_lower_res <- terra::rast("closeness_nt_spatial_raster_lower_res.tif")
closeness_t_spatial_raster_lower_res <- terra::rast("closeness_t_spatial_raster_lower_res.tif")
#
n_cells <- ncell(closeness_nt_spatial_raster_lower_res)

# Get the coordinates of the pixel centroids
closeness_nt_spatial_raster_lower_res_COORDS <- terra::xyFromCell(closeness_nt_spatial_raster_lower_res, 1:n_cells)
closeness_nt_spatial_raster_lower_res_VALUES <- terra::values(closeness_nt_spatial_raster_lower_res)
closeness_t_spatial_raster_lower_res_VALUES <- terra::values(closeness_t_spatial_raster_lower_res)

#Save
saveRDS(closeness_nt_spatial_raster_lower_res_COORDS, "closeness_nt_spatial_raster_lower_res_COORDS.rds")
saveRDS(closeness_nt_spatial_raster_lower_res_VALUES, "closeness_nt_spatial_raster_lower_res_VALUES.rds")
saveRDS(closeness_t_spatial_raster_lower_res_VALUES, "closeness_t_spatial_raster_lower_res_VALUES.rds")

#Load
#closeness_nt_spatial_raster_lower_res_COORDS <- readRDS("closeness_nt_spatial_raster_lower_res_COORDS.rds")
#closeness_nt_spatial_raster_lower_res_VALUES <- readRDS("closeness_nt_spatial_raster_lower_res_VALUES.rds")
#closeness_t_spatial_raster_lower_res_VALUES <- readRDS("closeness_t_spatial_raster_lower_res_VALUES.rds")

#Modified t-test
closeness_ttest <- modified.ttest(closeness_t_spatial_raster_lower_res_VALUES, 
                                  closeness_nt_spatial_raster_lower_res_VALUES, 
                                  closeness_nt_spatial_raster_lower_res_COORDS
)
#save
saveRDS(closeness_ttest, "closeness_ttest.rds")

##### indegree #####

indegree_nt_spatial_raster_lower_res <- terra::rast("indegree_nt_spatial_raster_lower_res.tif")
indegree_t_spatial_raster_lower_res <- terra::rast("indegree_t_spatial_raster_lower_res.tif")
#
n_cells <- ncell(indegree_nt_spatial_raster_lower_res)

# Get the coordinates of the pixel centroids
indegree_nt_spatial_raster_lower_res_COORDS <- terra::xyFromCell(indegree_nt_spatial_raster_lower_res, 1:n_cells)
indegree_nt_spatial_raster_lower_res_VALUES <- terra::values(indegree_nt_spatial_raster_lower_res)
indegree_t_spatial_raster_lower_res_VALUES <- terra::values(indegree_t_spatial_raster_lower_res)

#Save
saveRDS(indegree_nt_spatial_raster_lower_res_COORDS, "indegree_nt_spatial_raster_lower_res_COORDS.rds")
saveRDS(indegree_nt_spatial_raster_lower_res_VALUES, "indegree_nt_spatial_raster_lower_res_VALUES.rds")
saveRDS(indegree_t_spatial_raster_lower_res_VALUES, "indegree_t_spatial_raster_lower_res_VALUES.rds")

#Load
#indegree_nt_spatial_raster_lower_res_COORDS <- readRDS("indegree_nt_spatial_raster_lower_res_COORDS.rds")
#indegree_nt_spatial_raster_lower_res_VALUES <- readRDS("indegree_nt_spatial_raster_lower_res_VALUES.rds")
#indegree_t_spatial_raster_lower_res_VALUES <- readRDS("indegree_t_spatial_raster_lower_res_VALUES.rds")

#Modified t-test
indegree_ttest <- modified.ttest(indegree_t_spatial_raster_lower_res_VALUES, 
                                 indegree_nt_spatial_raster_lower_res_VALUES, 
                                 indegree_nt_spatial_raster_lower_res_COORDS
)
#save
saveRDS(indegree_ttest, "indegree_ttest.rds")

##### outdegree #####

outdegree_nt_spatial_raster_lower_res <- terra::rast("outdegree_nt_spatial_raster_lower_res.tif")
outdegree_t_spatial_raster_lower_res <- terra::rast("outdegree_t_spatial_raster_lower_res.tif")
#
n_cells <- ncell(outdegree_nt_spatial_raster_lower_res)

# Get the coordinates of the pixel centroids
outdegree_nt_spatial_raster_lower_res_COORDS <- terra::xyFromCell(outdegree_nt_spatial_raster_lower_res, 1:n_cells)
outdegree_nt_spatial_raster_lower_res_VALUES <- terra::values(outdegree_nt_spatial_raster_lower_res)
outdegree_t_spatial_raster_lower_res_VALUES <- terra::values(outdegree_t_spatial_raster_lower_res)

#Save
saveRDS(outdegree_nt_spatial_raster_lower_res_COORDS, "outdegree_nt_spatial_raster_lower_res_COORDS.rds")
saveRDS(outdegree_nt_spatial_raster_lower_res_VALUES, "outdegree_nt_spatial_raster_lower_res_VALUES.rds")
saveRDS(outdegree_t_spatial_raster_lower_res_VALUES, "outdegree_t_spatial_raster_lower_res_VALUES.rds")

#Load
#outdegree_nt_spatial_raster_lower_res_COORDS <- readRDS("outdegree_nt_spatial_raster_lower_res_COORDS.rds")
#outdegree_nt_spatial_raster_lower_res_VALUES <- readRDS("outdegree_nt_spatial_raster_lower_res_VALUES.rds")
#outdegree_t_spatial_raster_lower_res_VALUES <- readRDS("outdegree_t_spatial_raster_lower_res_VALUES.rds")

#Modified t-test
outdegree_ttest <- modified.ttest(outdegree_t_spatial_raster_lower_res_VALUES, 
                                  outdegree_nt_spatial_raster_lower_res_VALUES, 
                                  outdegree_nt_spatial_raster_lower_res_COORDS
)
#save
saveRDS(outdegree_ttest, "outdegree_ttest.rds")
