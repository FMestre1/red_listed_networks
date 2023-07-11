#FMestre
#10-07-2023

library(raster)
library(terra)
library(rasterVis)
library(cheddar)

#Load all rasters
nt_indegree <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_indegree.tif")
t_indegree <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_indegree.tif")
#
nt_outdegree <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_outdegree.tif")
t_outdegree <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_outdegree.tif")
#
nt_t_level <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_trophic_level.tif")
t_t_level <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_trophic_level.tif")
#
nt_closeness <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_closeness.tif")
t_closeness <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_closeness.tif")
#
nt_centrality <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_centrality.tif")
t_centrality <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_centrality.tif")
#
nt_ivi <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_ivi.tif")
t_ivi <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_ivi.tif")
#
proportion <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\proportion.tif")
#

richness <- terra::vect("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\shapefiles\\richness_2.shp")
terra::plot(richness, "sp_rich")
#richness_2 <- terra::rasterize(richness, proportion, "Count")
#richness$sp_richnes
#rasterVis::levelplot(richness_2)

#Plot with rastervis
rasterVis::levelplot(nt_indegree)
rasterVis::levelplot(t_indegree)
#
rasterVis::levelplot(nt_outdegree)
rasterVis::levelplot(t_outdegree)
#
rasterVis::levelplot(nt_t_level)
rasterVis::levelplot(t_t_level)
#
rasterVis::levelplot(nt_closeness)
rasterVis::levelplot(t_closeness)
#
rasterVis::levelplot(nt_centrality)
rasterVis::levelplot(t_centrality)
#
rasterVis::levelplot(nt_ivi)
rasterVis::levelplot(t_ivi)
#
rasterVis::levelplot(proportion)

#Plot side by side
par(mfrow=c(1, 2))
terra::plot(t_indegree, range = c(0, 65), main = "Threatened")
terra::plot(nt_indegree, range = c(0, 65), main = "Not threatened")
#
par(mfrow=c(1, 2))
terra::plot(t_outdegree, range = c(0, 25), main = "Threatened")
terra::plot(nt_outdegree, range = c(0, 25), main = "Not threatened")
#
par(mfrow=c(1, 2))
terra::plot(t_t_level, range = c(0, 3), main = "Threatened")
terra::plot(nt_t_level, range = c(0, 3), main = "Not threatened")
#
par(mfrow=c(1, 2))
terra::plot(t_centrality, range = c(0, 460), main = "Threatened")
terra::plot(nt_centrality, range = c(0, 460), main = "Not threatened")
#
par(mfrow=c(1, 2))
terra::plot(t_ivi, range = c(0, 100), main = "Threatened")
terra::plot(nt_ivi, range = c(0, 100), main = "Not threatened")
#
par(mfrow=c(1, 2))
terra::plot(t_closeness, range = c(0, 1), main = "Threatened")
terra::plot(nt_closeness, range = c(0, 1), main = "Not threatened")


#Evaluate statistical differences between threatened and non-threatened species
# Extract pixel values from the two rasters
t_indegree_values <- values(t_indegree,  mat=TRUE)
nt_indegree_values <- values(nt_indegree,  mat=TRUE)

# Perform the t-test
in_degree_result <- t.test(t_indegree_values, nt_indegree_values)

# Print the t-test results
print(in_degree_result)

?terra::values

