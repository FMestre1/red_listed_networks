#FMestre
#10-07-2023

library(raster)
library(terra)
library(rasterVis)
library(cheddar)
library(effsize)

################################################################################
#                                  Plot maps
################################################################################

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
#terra::plot(richness, "sp_rich")
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

################################################################################
# Evaluate statistical differences between threatened and non-threatened species
################################################################################

#t test interpretation
#p < 0.05 significance difference

nt_indegree <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/indegree_nt_spatial.shp")
t_indegree <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/indegree_t_spatial.shp")
nt_indegree <- as.data.frame(nt_indegree)
t_indegree <- as.data.frame(t_indegree)
nt_indegree <- nt_indegree[,c(1,4)]
t_indegree <- t_indegree[,c(1,4)]
names(nt_indegree)[2] <- "NT_indegree"
names(t_indegree)[2] <- "T_indegree"
indegree_compare <- merge(nt_indegree, t_indegree)
indegree_compare <- indegree_compare[complete.cases(indegree_compare),]
#View(indegree_compare)
indegree_ttest <- t.test(indegree_compare[,2], indegree_compare[,3])
print(indegree_ttest)
indegree_cohens_d <- effsize::cohen.d(indegree_compare[,2], indegree_compare[,3])
print(indegree_cohens_d)
plot(indegree_compare[,2], indegree_compare[,3])
cor(indegree_compare[,2], indegree_compare[,3], method = "spearman")
#
nt_outdegree <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/outdegree_nt_spatial.shp")
t_outdegree <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/outdegree_t_spatial.shp")
nt_outdegree <- as.data.frame(nt_outdegree)
t_outdegree <- as.data.frame(t_outdegree)
nt_outdegree <- nt_outdegree[,c(1,4)]
t_outdegree <- t_outdegree[,c(1,4)]
names(nt_outdegree)[2] <- "NT_outdegree"
names(t_outdegree)[2] <- "T_outdegree"
outdegree_compare <- merge(nt_outdegree, t_outdegree)
outdegree_compare <- outdegree_compare[complete.cases(outdegree_compare),]
#View(outdegree_compare)
outdegree_ttest <- t.test(outdegree_compare[,2], outdegree_compare[,3])
print(outdegree_ttest)
outdegree_cohens_d <- effsize::cohen.d(outdegree_compare[,2], outdegree_compare[,3])
print(outdegree_cohens_d)
plot(outdegree_compare[,2], outdegree_compare[,3])
cor(outdegree_compare[,2], outdegree_compare[,3], method = "spearman")
#
nt_t_level <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/tl_nt_spatial.shp")
t_t_level <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/tl_t_spatial.shp")
nt_t_level <- as.data.frame(nt_t_level)
t_t_level <- as.data.frame(t_t_level)
nt_t_level <- nt_t_level[,c(1,4)]
t_t_level <- t_t_level[,c(1,4)]
names(nt_t_level)[2] <- "NT_trph_level"
names(t_t_level)[2] <- "T_trph_level"
trophic_level_compare <- merge(nt_t_level, t_t_level)
trophic_level_compare <- trophic_level_compare[complete.cases(trophic_level_compare),]
#View(trophic_level_compare)
trophic_level_ttest <- t.test(trophic_level_compare[,2], trophic_level_compare[,3])
print(trophic_level_ttest)
trophic_level_cohens_d <- effsize::cohen.d(trophic_level_compare[,2], trophic_level_compare[,3])
print(trophic_level_cohens_d)
plot(trophic_level_compare[,2], trophic_level_compare[,3])
cor(trophic_level_compare[,2], trophic_level_compare[,3], method = "spearman")
#
nt_closeness <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/closeness_nt_spatial.shp")
t_closeness <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/closeness_t_spatial.shp")
nt_closeness <- as.data.frame(nt_closeness)
t_closeness <- as.data.frame(t_closeness)
nt_closeness <- nt_closeness[,c(1,4)]
t_closeness <- t_closeness[,c(1,4)]
names(nt_closeness)[2] <- "NT_closeness"
names(t_closeness)[2] <- "T_closeness"
closeness_compare <- merge(nt_closeness, t_closeness)
closeness_compare <- closeness_compare[complete.cases(closeness_compare),]
#View(closeness_compare)
closeness_ttest <- t.test(closeness_compare[,2], closeness_compare[,3])
print(closeness_ttest)
outdegree_cohens_d <- effsize::cohen.d(closeness_compare[,2], closeness_compare[,3])
print(outdegree_cohens_d)
plot(closeness_compare[,2], closeness_compare[,3])
cor(closeness_compare[,2], closeness_compare[,3], method = "spearman")
#
nt_centrality <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/centrality_nt_spatial.shp")
t_centrality <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/centrality_t_spatial.shp")
nt_centrality <- as.data.frame(nt_centrality)
t_centrality <- as.data.frame(t_centrality)
nt_centrality <- nt_centrality[,c(1,4)]
t_centrality <- t_centrality[,c(1,4)]
names(nt_centrality)[2] <- "NT_centrality"
names(t_centrality)[2] <- "T_centrality"
centrality_compare <- merge(nt_centrality, t_centrality)
centrality_compare <- centrality_compare[complete.cases(centrality_compare),]
#View(centrality_compare)
centrality_ttest <- t.test(centrality_compare[,2], centrality_compare[,3])
print(centrality_ttest)
centrality_cohens_d <- effsize::cohen.d(centrality_compare[,2], centrality_compare[,3])
print(centrality_cohens_d)
plot(centrality_compare[,2], centrality_compare[,3])
cor(centrality_compare[,2], centrality_compare[,3], method = "spearman")
#
nt_ivi <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/ivi_nt_spatial_second_version.shp")
t_ivi <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/ivi_t_spatial_second_version.shp")
nt_ivi <- as.data.frame(nt_ivi)
t_ivi <- as.data.frame(t_ivi)
nt_ivi <- nt_ivi[,c(1,4)]
t_ivi <- t_ivi[,c(1,4)]
names(nt_ivi)[2] <- "NT_ivi"
names(t_ivi)[2] <- "T_ivi"
ivi_compare <- merge(nt_ivi, t_ivi)
ivi_compare <- ivi_compare[complete.cases(ivi_compare),]
#View(centrality_compare)
ivi_ttest <- t.test(ivi_compare[,2], ivi_compare[,3])
print(ivi_ttest)
ivi_cohens_d <- effsize::cohen.d(ivi_compare[,2], ivi_compare[,3])
print(ivi_cohens_d)
plot(ivi_compare[,2], ivi_compare[,3])
cor(ivi_compare[,2], ivi_compare[,3], method = "spearman")
#
proportion <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/proportion_spatial.shp")
richness <- terra::vect("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\shapefiles\\richness_2.shp")












