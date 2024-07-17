################################################################################
#                             Evaluate differences           
################################################################################

#FMestre
#23-04-2024

#load packages
library(terra)
library(rasterVis)
library(diffeR)

#Load Rasters
#nt_ivi <- terra::rast("rasters_15JUL\\nt_ivi_15JUL.tif")
#t_ivi <- terra::rast("rasters_15JUL\\t_ivi_15JUL.tif")
#t_centrality <- terra::rast("rasters_15JUL\\t_centrality_15JUL.tif")
#nt_centrality <- terra::rast("rasters_15JUL\\nt_centrality_15JUL.tif")
#t_indegree <- terra::rast("rasters_15JUL\\t_indegree_15JUL.tif")
#nt_indegree <- terra::rast("rasters_15JUL\\nt_indegree_15JUL.tif")
#nt_outdegree <- terra::rast("rasters_15JUL\\nt_outdegree_15JUL.tif")
#t_outdegree <- terra::rast("rasters_15JUL\\t_outdegree_15JUL.tif")
#t_closeness <- terra::rast("rasters_15JUL\\t_closeness_15JUL.tif")
#nt_closeness <- terra::rast("rasters_15JUL\\nt_closeness_15JUL.tif")
#t_tl <- terra::rast("rasters_15JUL\\t_tl_15JUL.tif")
#nt_tl <- terra::rast("rasters_15JUL\\nt_tl_15JUL.tif")
#proportion_r <- terra::rast("rasters_15JUL\\proportion_r__15JUL.tif")
#
#indeg_diff <- terra::rast("rasters_15JUL\\indeg_diff.tif")
#outdeg_diff <- terra::rast("rasters_15JUL\\outdeg_diff.tif")
#trophic_level_diff <- terra::rast("rasters_15JUL\\trophic_level_diff.tif")
#closeness_diff <- terra::rast("rasters_15JUL\\closeness_diff.tif")
#centrality_diff <- terra::rast("rasters_15JUL\\centrality_diff.tif")
#ivi_diff <- terra::rast("rasters_15JUL\\ivi_diff.tif")

################################################################################
#                           Differences between maps
################################################################################

#?diffeR::differenceMR

# Calculate the difference between the two maps
indegree_diff <-t_indegree - nt_indegree
writeRaster(indegree_diff, "rasters_15JUL\\indegree_diff.tif", overwrite=TRUE)

outdegree_diff <- t_outdegree - nt_outdegree
writeRaster(outdegree_diff, "rasters_15JUL\\outdegree_diff.tif", overwrite=TRUE)

t_level_diff <- t_tl - nt_tl
writeRaster(t_level_diff, "rasters_15JUL\\t_level_diff.tif", overwrite=TRUE)

closeness_diff <- t_closeness - nt_closeness
writeRaster(closeness_diff, "rasters_15JUL\\closeness_diff.tif", overwrite=TRUE)

centrality_diff <- t_centrality - nt_centrality
writeRaster(centrality_diff, "rasters_15JUL\\centrality_diff.tif", overwrite=TRUE)

ivi_diff <- t_ivi - nt_ivi
writeRaster(ivi_diff, "rasters_15JUL\\ivi_diff.tif", overwrite=TRUE)

################################################################################
#                   Are the maps significantly different?    
################################################################################

#From: 
#https://sesync-ci.github.io/blog/raster-change-analysis.html

#indeg_diff_mean <- as.numeric(terra::global(indegree_diff, "mean", na.rm=TRUE))# mean
#indeg_diff_sd <- as.numeric(terra::global(indegree_diff, "sd", na.rm=TRUE))# sd
#indegree_diff_std <- (indegree_diff - indeg_diff_mean)/indeg_diff_sd # standardized image
#terra::writeRaster(indegree_diff_std, "rasters_15JUL\\indegree_diff_std.tif")
#terra::plot(indegree_diff_std)
#hist(indegree_diff_std, main="Standardized difference", xlab="Difference")

#(...)

################################################################################
#    Extract examples for the results & Discussion
################################################################################

network_list <- list.files("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_3")
pat_dataset <- "C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_3"

#Svaldbard

#C422
which(network_list == "C422.csv")
read.csv(paste0(pat_dataset, "\\C422.csv"))

#IVI
terra::global(ivi_diff, "mean", na.rm =TRUE)

#Closeness_diff
terra::global(closeness_diff, "mean", na.rm =TRUE)

#(...)

################################################################################
#    Combine
################################################################################

#indeg_diff
#outdeg_diff
#trophic_level_diff
#closeness_diff
#centrality_diff
#ivi_diff

# Define minimum and maximum values (adjust based on your data)
min_val_indeg_diff <- as.numeric(terra::global(indeg_diff, fun = "min", na.rm = TRUE))
max_val_indeg_diff <- as.numeric(terra::global(indeg_diff, fun = "max", na.rm = TRUE))
# Normalize values (0 to 1)
indeg_diff_normalized_rast <- (indeg_diff - min_val_indeg_diff) / (max_val_indeg_diff - min_val_indeg_diff)
#plot(indeg_diff_normalized_rast)

#Define minimum and maximum values (adjust based on your data)
min_val_outdeg_diff <- as.numeric(terra::global(outdeg_diff, fun = "min", na.rm = TRUE))
max_val_outdeg_diff <- as.numeric(terra::global(outdeg_diff, fun = "max", na.rm = TRUE))
# Normalize values (0 to 1)
outdeg_diff_normalized_rast <- (outdeg_diff - min_val_outdeg_diff) / (max_val_outdeg_diff - min_val_outdeg_diff)
#plot(outdeg_diff_normalized_rast)

#Define minimum and maximum values (adjust based on your data)
min_val_trophic_level_diff <- as.numeric(terra::global(trophic_level_diff, fun = "min", na.rm = TRUE))
max_val_trophic_level_diff <- as.numeric(terra::global(trophic_level_diff, fun = "max", na.rm = TRUE))
# Normalize values (0 to 1)
trophic_level_diff_normalized_rast <- (trophic_level_diff - min_val_trophic_level_diff) / (max_val_trophic_level_diff - min_val_trophic_level_diff)
#plot(trophic_level_diff_normalized_rast)

#Define minimum and maximum values (adjust based on your data)
min_val_closeness_diff <- as.numeric(terra::global(closeness_diff, fun = "min", na.rm = TRUE))
max_val_closeness_diff <- as.numeric(terra::global(closeness_diff, fun = "max", na.rm = TRUE))
# Normalize values (0 to 1)
closeness_diff_normalized_rast <- (closeness_diff - min_val_closeness_diff) / (max_val_closeness_diff - min_val_closeness_diff)
#plot(closeness_diff_normalized_rast)

#Define minimum and maximum values (adjust based on your data)
min_val_centrality_diff <- as.numeric(terra::global(centrality_diff, fun = "min", na.rm = TRUE))
max_val_centrality_diff <- as.numeric(terra::global(centrality_diff, fun = "max", na.rm = TRUE))
# Normalize values (0 to 1)
centrality_diff_normalized_rast <- (centrality_diff - min_val_centrality_diff) / (max_val_centrality_diff - min_val_centrality_diff)
#plot(centrality_diff_normalized_rast)

#Define minimum and maximum values (adjust based on your data)
min_val_ivi_diff <- as.numeric(terra::global(ivi_diff, fun = "min", na.rm = TRUE))
max_val_ivi_diff <- as.numeric(terra::global(ivi_diff, fun = "max", na.rm = TRUE))
# Normalize values (0 to 1)
ivi_diff_normalized_rast <- (ivi_diff - min_val_ivi_diff) / (max_val_ivi_diff - min_val_ivi_diff)
#plot(ivi_diff_normalized_rast)

all_diffs <- c(indeg_diff_normalized_rast,
              outdeg_diff_normalized_rast,
              trophic_level_diff_normalized_rast,
              closeness_diff_normalized_rast,
              centrality_diff_normalized_rast,
              ivi_diff_normalized_rast
              )

mean_all_diffs <- terra::app(all_diffs, mean)
plot(mean_all_diffs)
rasterVis::levelplot(mean_all_diffs, par.settings=BuRdTheme())

#Save results
#terra::writeRaster(mean_all_diffs, "mean_all_diffs.tif")
