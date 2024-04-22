################################################################################
#                             Evaluate differences           
################################################################################

#FMestre
#23-04-2024

#load packages
library(terra)
library(rasterVis)

#Load Rasters
nt_ivi <- terra::rast("rasters_21OUT\\nt_ivi_20OUT.tif")
t_ivi <- terra::rast("rasters_21OUT\\t_ivi_20OUT.tif")
t_centrality <- terra::rast("rasters_21OUT\\t_centrality_20OUT.tif")
nt_centrality <- terra::rast("rasters_21OUT\\nt_centrality_20OUT.tif")
t_indegree <- terra::rast("rasters_21OUT\\t_indegree_20OUT.tif")
nt_indegree <- terra::rast("rasters_21OUT\\nt_indegree_20OUT.tif")
nt_outdegree <- terra::rast("rasters_21OUT\\nt_outdegree_20OUT.tif")
t_outdegree <- terra::rast("rasters_21OUT\\t_outdegree_20OUT.tif")
t_closeness <- terra::rast("rasters_21OUT\\t_closeness_20OUT.tif")
nt_closeness <- terra::rast("rasters_21OUT\\nt_closeness_20OUT.tif")
t_tl <- terra::rast("rasters_21OUT\\t_tl_20OUT.tif")
nt_tl <- terra::rast("rasters_21OUT\\nt_tl_20OUT.tif")
proportion_r <- terra::rast("rasters_21OUT\\proportion_r_20OUT.tif")
#
indeg_diff <- terra::rast("rasters_21OUT\\indeg_diff.tif")
outdeg_diff <- terra::rast("rasters_21OUT\\outdeg_diff.tif")
trophic_level_diff <- terra::rast("rasters_21OUT\\trophic_level_diff.tif")
closeness_diff <- terra::rast("rasters_21OUT\\closeness_diff.tif")
centrality_diff <- terra::rast("rasters_21OUT\\centrality_diff.tif")
ivi_diff <- terra::rast("rasters_21OUT\\ivi_diff.tif")



################################################################################
#                   Are the maps significantly different?    
################################################################################

#FMestre
#20-04-2024

#From: 
#https://sesync-ci.github.io/blog/raster-change-analysis.html

indeg_diff_mean <- as.numeric(terra::global(indeg_diff, "mean", na.rm=TRUE))# mean
indeg_diff_sd <- as.numeric(terra::global(indeg_diff, "sd", na.rm=TRUE))# sd
indeg_diff_std <- (indeg_diff - indeg_diff_mean)/indeg_diff_sd # standardized image
terra::writeRaster(indeg_diff_std, "indeg_diff_std.tif")

threshold_val <- c(1.96,1.64)
plot(indeg_diff_std)

hist(indeg_diff_std,
     main="Standardized difference",
     xlab="Difference")

abline(v=indeg_diff_mean + indeg_diff_std,col="red",lty=1)
abline(v=indeg_diff_mean - indeg_diff_std,col="red",lty=1)

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

