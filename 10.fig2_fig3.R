################################################################################
#                       Fig2 - Maps of average metrics           
################################################################################

#FMestre

#Load packages
library(terra)
library(viridis)
library(rasterVis)
library(gridExtra)
library(grid)

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

#Reclassify to quartiles - values obtained from QGIS (to derive the quantiles here got "Error: std::bad_alloc")

### IVI#########################################################################
nt_ivi_reclass_matrix <- matrix(c(1, 3.496, 1,    # First 
                                  3.496, 4.987, 2,  # Second 
                                  4.987, 6.547, 3,   # Third 
                                  6.547, 72.1, 4), # Fourth
                                  ncol = 3, 
                                  byrow = TRUE)
nt_ivi_reclassified <- terra::classify(nt_ivi, nt_ivi_reclass_matrix)
plot(nt_ivi_reclassified)
terra::writeRaster(nt_ivi_reclassified, "nt_ivi_reclassified.tif", overwrite = TRUE)
#
t_ivi_reclass_matrix <- matrix(c(1, 1.785, 1,    # First 
                                 1.785, 2.446, 2,  # Second 
                                 2.446, 3.753, 3,   # Third 
                                 3.753, 100, 4), # Fourth
                                ncol = 3, 
                                byrow = TRUE)
t_ivi_reclassified <- terra::classify(t_ivi, t_ivi_reclass_matrix)
plot(t_ivi_reclassified)
terra::writeRaster(t_ivi_reclassified, "t_ivi_reclassified.tif", overwrite = TRUE)

###CENTRALITY###################################################################
nt_centrality_reclass_matrix <- matrix(c(-1, 9.903, 1,    # First 
                                         9.903, 22.337, 2,  # Second 
                                         22.337, 33.264, 3,   # Third 
                                         33.264, 111.71, 4), # Fourth
                                       ncol = 3, 
                                       byrow = TRUE)
nt_centrality_reclassified <- terra::classify(nt_centrality, nt_centrality_reclass_matrix)
plot(nt_centrality_reclassified)
terra::writeRaster(nt_centrality_reclassified, "nt_centrality_reclassified.tif", overwrite = TRUE)
#
t_centrality_reclass_matrix <- matrix(c(-1, 0, 1,    # First 
                                        0, 0.880, 2,  # Second 
                                        0.880, 2.667, 3,   # Third 
                                        2.667, 458.62, 4), # Fourth
                                      ncol = 3, 
                                      byrow = TRUE)
t_centrality_reclassified <- terra::classify(t_centrality, t_centrality_reclass_matrix)
plot(t_centrality_reclassified)
terra::writeRaster(t_centrality_reclassified, "t_centrality_reclassified.tif", overwrite = TRUE)


###CLOSENESS####################################################################
nt_closeness_reclass_matrix <- matrix(c(0, 0.000025, 1,    # First 
                                        0.000025, 0.000032, 2,  # Second 
                                        0.000032, 0.000052, 3,   # Third 
                                        0.000052, 1, 4), # Fourth
                                      ncol = 3, 
                                      byrow = TRUE)
nt_closeness_reclassified <- terra::classify(nt_closeness, nt_closeness_reclass_matrix)
plot(nt_closeness_reclassified)
terra::writeRaster(nt_closeness_reclassified, "nt_closeness_reclassified.tif", overwrite = TRUE)
#
t_closeness_reclass_matrix <- matrix(c(0, 0.000023, 1,    # First 
                                        0.000023, 0.000031, 2,  # Second 
                                        0.000031, 0.000049, 3,   # Third 
                                        0.000049, 0.75, 4), # Fourth
                                     ncol = 3, 
                                     byrow = TRUE)
t_closeness_reclassified <- terra::classify(t_closeness, t_closeness_reclass_matrix)
plot(t_closeness_reclassified)
terra::writeRaster(t_closeness_reclassified, "t_closeness_reclassified.tif")


###INDEGREE#####################################################################
nt_indegree_reclass_matrix <- matrix(c(-1, 12.497, 1,    # First 
                                       12.497, 14.915, 2,  # Second 
                                       14.915, 16.442, 3,   # Third 
                                       16.442, 22.8, 4), # Fourth
                                     ncol = 3, 
                                     byrow = TRUE)
nt_indegree_reclassified <- terra::classify(nt_indegree, nt_indegree_reclass_matrix)
plot(nt_indegree_reclassified)
terra::writeRaster(nt_indegree_reclassified, "nt_indegree_reclassified.tif", overwrite = TRUE)
######
t_indegree_reclass_matrix <- matrix(c(-1, 0, 1,    # First 
                                      0, 11, 2,  # Second 
                                      11, 18.25, 3,   # Third 
                                      18.25, 61, 4), # Fourth
                                    ncol = 3, 
                                    byrow = TRUE)
t_indegree_reclassified <- terra::classify(t_indegree, t_indegree_reclass_matrix)
plot(t_indegree_reclassified)
terra::writeRaster(t_indegree_reclassified, "t_indegree_reclassified.tif", overwrite = TRUE)


###OUTDEGREE####################################################################
nt_outdegree_reclass_matrix <- matrix(c(0, 12.429, 1,    # First 
                                        12.429, 14.963, 2,  # Second 
                                        14.963, 16.587, 3,   # Third 
                                        16.587, 22.23, 4), # Fourth
                                      ncol = 3, 
                                      byrow = TRUE)
nt_outdegree_reclassified <- terra::classify(nt_outdegree, nt_outdegree_reclass_matrix)
plot(nt_outdegree_reclassified)
terra::writeRaster(nt_outdegree_reclassified, "nt_outdegree_reclassified.tif", overwrite = TRUE)
#
t_outdegree_reclass_matrix <- matrix(c(0, 8.5, 1,    # First 
                                       8.5, 11.75, 2,  # Second 
                                       11.75, 14.22, 3,   # Third 
                                       14.22, 33, 4), # Fourth
                                     ncol = 3, 
                                     byrow = TRUE)
t_outdegree_reclassified <- terra::classify(t_outdegree, t_outdegree_reclass_matrix)
plot(t_outdegree_reclassified)
terra::writeRaster(t_outdegree_reclassified, "t_outdegree_reclassified.tif", overwrite = TRUE)


####TL
nt_tl_reclass_matrix <- matrix(c(1, 1.376, 1,    # First 
                                 1.376, 1.396, 2,  # Second 
                                 1.396, 1.417, 3,   # Third 
                                 1.417, 2, 4), # Fourth
                               ncol = 3, 
                               byrow = TRUE)
nt_tl_reclassified <- terra::classify(nt_tl, nt_tl_reclass_matrix)
plot(nt_tl_reclassified)
terra::writeRaster(nt_tl_reclassified, "nt_tl_reclassified.tif")
#
t_tl_reclass_matrix <- matrix(c(1, 1, 1,    # First 
                                 1, 1.333, 2,  # Second 
                                1.333, 1.500, 3,   # Third 
                                1.500, 3, 4), # Fourth
                               ncol = 3, 
                               byrow = TRUE)
t_tl_reclassified <- terra::classify(t_tl, t_tl_reclass_matrix)
plot(t_tl_reclassified)
terra::writeRaster(t_tl_reclassified, "t_tl_reclassified.tif")


#Load
nt_ivi_reclassified <- terra::rast("nt_ivi_reclassified.tif")
t_ivi_reclassified <- terra::rast("t_ivi_reclassified.tif")
nt_centrality_reclassified <- terra::rast("nt_centrality_reclassified.tif")
t_centrality_reclassified <- terra::rast("t_centrality_reclassified.tif")
nt_closeness_reclassified <- terra::rast("nt_closeness_reclassified.tif")
t_closeness_reclassified <- terra::rast("t_closeness_reclassified.tif")
nt_indegree_reclassified <- terra::rast("nt_indegree_reclassified.tif")
t_indegree_reclassified <- terra::rast("t_indegree_reclassified.tif")
nt_outdegree_reclassified <- terra::rast("nt_outdegree_reclassified.tif")
t_outdegree_reclassified <- terra::rast("t_outdegree_reclassified.tif")
nt_tl_reclassified <- terra::rast("nt_tl_reclassified.tif")
t_tl_reclassified <- terra::rast("t_tl_reclassified.tif")


#Colours
my.cols <- colorRampPalette(colors = c("#FFA54F", "#EE4000", "#CD2626", "#67001F"))(100)
min_max <- seq(1, 4, length.out = 100)

#Trophic Level
tl1_quant <- rasterVis::levelplot(nt_tl_reclassified, col.regions=my.cols, at=min_max, main = list("Not-threatened", cex=2), scales = list(draw = FALSE))
tl2_quant <- rasterVis::levelplot(t_tl_reclassified, col.regions=my.cols, at=min_max, main = list("Threatened", cex=2), scales = list(draw = FALSE))
grid.arrange(tl1_quant, tl2_quant, ncol=2)

#IVI
ivi1_quant <- rasterVis::levelplot(nt_ivi_reclassified, col.regions=my.cols, at=min_max, main = list("Not-threatened", cex=2), scales = list(draw = FALSE))
ivi2_quant <- rasterVis::levelplot(t_ivi_reclassified, col.regions=my.cols, at=min_max, main = list("Threatened", cex=2), scales = list(draw = FALSE))
grid.arrange(ivi1_quant, ivi2_quant, ncol=2)

#In-degree
ind1_quant <- rasterVis::levelplot(nt_indegree_reclassified, col.regions=my.cols, at=min_max, main = list("Not-threatened", cex=2), scales = list(draw = FALSE))
ind2_quant <- rasterVis::levelplot(t_indegree_reclassified, col.regions=my.cols, at=min_max, main = list("Threatened", cex=2), scales = list(draw = FALSE))
grid.arrange(ind1_quant, ind2_quant, ncol=2)

#Out-degree
outd1_quant <- rasterVis::levelplot(nt_outdegree_reclassified, col.regions=my.cols, at=min_max, main = list("Not-threatened", cex=2), scales = list(draw = FALSE))
outd2_quant <- rasterVis::levelplot(t_outdegree_reclassified, col.regions=my.cols, at=min_max, main = list("Threatened", cex=2), scales = list(draw = FALSE))
grid.arrange(outd1_quant, outd2_quant, ncol=2)

#Closeness Centrality
cl1_quant <- rasterVis::levelplot(nt_closeness_reclassified, col.regions=my.cols, at=min_max, main = list("Not-threatened", cex=2), scales = list(draw = FALSE))
cl2_quant <- rasterVis::levelplot(t_closeness_reclassified, col.regions=my.cols, at=min_max, main = list("Threatened", cex=2), scales = list(draw = FALSE))
grid.arrange(cl1_quant, cl2_quant, ncol=2)

#Betweenness Centrality
bt1_quant <- rasterVis::levelplot(nt_centrality_reclassified, col.regions=my.cols, at=min_max, main = list("Not-threatened", cex=2), scales = list(draw = FALSE))
bt2_quant <- rasterVis::levelplot(t_centrality_reclassified, col.regions=my.cols, at=min_max, main = list("Threatened", cex=2), scales = list(draw = FALSE))
grid.arrange(bt1_quant, bt2_quant, ncol=2)
