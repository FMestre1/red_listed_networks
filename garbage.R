
#################################################
# process in parallel
library(doParallel) 
#library(spatialEco)
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


################################################################################

#Use smaller area to test function

portugal <- terra::vect("C:/Users/mestr/Documents/0. Artigos/IUCN_networks/shapefiles/portugal.shp")
#terra::crs(portugal)
portugal_p <- terra::project(portugal, ivi_nt_spatial_raster)
#terra::crs(portugal_p)
#plot(portugal)
#terra::crs(ivi_nt_spatial_raster)

#plot(ivi_nt_spatial_raster)
#plot(portugal_p, add=TRUE)

ivi_nt_spatial_raster_PT <- terra::crop(ivi_nt_spatial_raster, portugal_p)
#plot(ivi_nt_spatial_raster_PT)
ivi_t_spatial_raster_PT <- terra::crop(ivi_t_spatial_raster, portugal_p)
#plot(ivi_t_spatial_raster_PT)

#ivi_compare <- raster.modified.ttest(
#  ivi_nt_spatial_raster_PT, 
#  ivi_t_spatial_raster_PT, 
#  sample = "random",
#  p = 0.01
#)

####
#Compare maps
#plot(terra::compare(ivi_nt_spatial_raster_PT, ivi_t_spatial_raster_PT, oper = "<"))
#terra::compareGeom(ivi_nt_spatial_raster_PT, ivi_t_spatial_raster_PT)

#Using diffeR
#?terra::compare

#test <- diffeR::differenceMR(ivi_nt_spatial_raster_PT, ivi_t_spatial_raster_PT, eval = "original")
#test <- diffeR::differenceMR(ivi_nt_spatial_raster, ivi_t_spatial_raster)
#test

# process in parallel
#library(doParallel) 
#detectCores()
#cl <- makeCluster(4, type='PSOCK') # number of cores adjusted to the total number (maybe not use all!)
#registerDoParallel(cl)

#test_ivi <- spatialEco::raster.change(ivi_nt_spatial_raster_PT, ivi_t_spatial_raster_PT, stat = "t.test")


#################

#cor.test(as.vector(ivi_nt_spatial_raster_PT), as.vector(ivi_t_spatial_raster_PT))

#################



library(ade4)

tl_nt_spatial_raster
tl_t_spatial_raster

# Reduce resolution by a factor of 2 (or more if needed)
tl_nt_spatial_raster_resampled <- aggregate(tl_nt_spatial_raster, fact = 2)
tl_t_spatial_raster_resampled <- aggregate(tl_t_spatial_raster, fact = 2)

tl_ssi_compare_resample <- SSIMmap::ssim_raster(tl_nt_spatial_raster_resampled, 
                                                tl_t_spatial_raster_resampled, 
                                                global = FALSE)

plot_tl_SSIM <- rasterVis::levelplot(tl_ssi_compare_resample[[1]], contour=FALSE, par.settings = myTheme, main= "SSIM", margin = FALSE)
plot_tl_SIM <- rasterVis::levelplot(tl_ssi_compare_resample[[2]], contour=FALSE, par.settings = myTheme, main= "SIM", margin = FALSE)
plot_tl_SIV <- rasterVis::levelplot(tl_ssi_compare_resample[[3]], contour=FALSE, par.settings = myTheme, main= "SIV", margin = FALSE)
plot_tl_SIP <- rasterVis::levelplot(tl_ssi_compare_resample[[4]], contour=FALSE, par.settings = myTheme, main= "SIP", margin = FALSE)

gridExtra::grid.arrange(plot_tl_SSIM,
                        plot_tl_SIM,
                        plot_tl_SIV,
                        plot_tl_SIP,
                        ncol=2,
                        top=grid::textGrob("Trophic Level"))

# Assuming map1 and map2 are matrices or rasters
#dist_map1 <- dist(as.vector(tl_nt_spatial_raster_resampled))
#dist_map2 <- dist(as.vector(tl_t_spatial_raster_resampled))

# Perform Mantel test
mantel_result <- mantel.rtest(dist_map1, dist_map2, nrepet = 999)

##############################################

library(spdep)

# Create a spatial neighborhood structure (based on map coordinates)
coords <- coordinates(map1)  # Assuming map1 is a raster
nb <- dnearneigh(coords, 0, 1.5)  # Define neighbors within a certain distance
lw <- nb2listw(nb, style = "W")

# Compute the actual correlation or similarity
actual_similarity <- cor(map1_vals, map2_vals)

# Perform a spatially constrained permutation test
n_permutations <- 999
similarity_scores <- numeric(n_permutations)

for (i in 1:n_permutations) {
  # Perform spatially constrained permutation
  permuted_map2_vals <- spsample(map2_vals, lw)
  
  # Compute the similarity for the permuted map
  similarity_scores[i] <- cor(map1_vals, permuted_map2_vals)
}

# Calculate the p-value
p_value <- sum(similarity_scores >= actual_similarity) / n_permutations

print(paste("Spatially constrained permutation test p-value:", p_value))


