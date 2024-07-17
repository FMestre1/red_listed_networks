################################################################################
#                       Fig2 - Maps of average metrics           
################################################################################

#FMestre

#Load Rasters
nt_ivi <- terra::rast("rasters_15JUL\\nt_ivi_15JUL.tif")
t_ivi <- terra::rast("rasters_15JUL\\t_ivi_15JUL.tif")
t_centrality <- terra::rast("rasters_15JUL\\t_centrality_15JUL.tif")
nt_centrality <- terra::rast("rasters_15JUL\\nt_centrality_15JUL.tif")
t_indegree <- terra::rast("rasters_15JUL\\t_indegree_15JUL.tif")
nt_indegree <- terra::rast("rasters_15JUL\\nt_indegree_15JUL.tif")
nt_outdegree <- terra::rast("rasters_15JUL\\nt_outdegree_15JUL.tif")
t_outdegree <- terra::rast("rasters_15JUL\\t_outdegree_15JUL.tif")
t_closeness <- terra::rast("rasters_15JUL\\t_closeness_15JUL.tif")
nt_closeness <- terra::rast("rasters_15JUL\\nt_closeness_15JUL.tif")
t_tl <- terra::rast("rasters_15JUL\\t_tl_15JUL.tif")
nt_tl <- terra::rast("rasters_15JUL\\nt_tl_15JUL.tif")
proportion_r <- terra::rast("rasters_15JUL\\proportion_r__15JUL.tif")
#
#Load packages
library(terra)
library(viridis)
library(rasterVis)
library(gridExtra)
library(grid)

#Colours
my.cols <- viridis_pal(option = "D")(100)

#Trophic Level
tl_min_max <- seq(0, 5, length.out = 100)
tl1 <- rasterVis::levelplot(nt_tl, col.regions=my.cols, at=tl_min_max, main = "Not-threatened", scales = list(draw = FALSE))
tl2 <- rasterVis::levelplot(t_tl, col.regions=my.cols, at=tl_min_max, main = "Threatened", scales = list(draw = FALSE))
grid.arrange(tl1, tl2, ncol=2)

#IVI
ivi_min_max <- seq(0, 100, length.out = 100)
ivi1 <- rasterVis::levelplot(nt_ivi, col.regions=my.cols, at=ivi_min_max, main = "Not-threatened", scales = list(draw = FALSE))
ivi2 <- rasterVis::levelplot(t_ivi, col.regions=my.cols, at=ivi_min_max, main = "Threatened", scales = list(draw = FALSE))
grid.arrange(ivi1, ivi2, ncol=2)

#In-degree
ind_min_max <- seq(0, 65, length.out = 100)
ind1 <- rasterVis::levelplot(nt_indegree, col.regions=my.cols, at=ind_min_max, main = "Not-threatened", scales = list(draw = FALSE))
ind2 <- rasterVis::levelplot(t_indegree, col.regions=my.cols, at=ind_min_max, main = "Threatened", scales = list(draw = FALSE))
grid.arrange(ind1, ind2, ncol=2)

#Out-degree
out_min_max <- seq(0, 35, length.out = 100)
outd1 <- rasterVis::levelplot(nt_outdegree, col.regions=my.cols, at=out_min_max, main = "Not-threatened", scales = list(draw = FALSE))
outd2 <- rasterVis::levelplot(t_outdegree, col.regions=my.cols, at=out_min_max, main = "Threatened", scales = list(draw = FALSE))
grid.arrange(outd1, outd2, ncol=2)

#Closeness Centrality
cl_min_max <- seq(0, 1, length.out = 100)
cl1 <- rasterVis::levelplot(nt_closeness, col.regions=my.cols, at=cl_min_max, main = "Not-threatened", scales = list(draw = FALSE))
cl2 <- rasterVis::levelplot(t_closeness, col.regions=my.cols, at=cl_min_max, main = "Threatened", scales = list(draw = FALSE))
grid.arrange(cl1, cl2, ncol=2)

#Betweenness Centrality
bt_min_max <- seq(0, 500, length.out = 100)
bt1 <- rasterVis::levelplot(nt_centrality, col.regions=my.cols, at=bt_min_max, main = "Not-threatened", scales = list(draw = FALSE))
bt2 <- rasterVis::levelplot(t_centrality, col.regions=my.cols, at=bt_min_max, main = "Threatened", scales = list(draw = FALSE))
grid.arrange(bt1, bt2, ncol=2)
