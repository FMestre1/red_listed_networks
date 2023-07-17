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

#######################################################################################
#    Evaluate statistical differences between threatened and non-threatened species
#######################################################################################

proportion <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/proportion_spatial.shp")
richness <- terra::vect("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\shapefiles\\richness_2.shp")

#T-test and Cohen's d

#t test interpretation
#p < 0.05 significance difference

#Coehen's d interpretation
#Small effect size: Cohen's d around 0.2 indicates a small effect size. 
#It suggests that there is a small difference between the groups, which may 
#have limited practical significance.

#Medium effect size: Cohen's d around 0.5 indicates a medium effect size. 
#It suggests a moderate difference between the groups, which is typically 
#considered to have a meaningful practical impact.

#Large effect size: Cohen's d around 0.8 or above indicates a large effect size. 
#It suggests a substantial difference between the groups, which is likely to 
#have significant practical implications.

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
indegree_ttest <- t.test(indegree_compare[,2], indegree_compare[,3], paired = TRUE)
print(indegree_ttest)
indegree_cohens_d <- effsize::cohen.d(indegree_compare[,2], indegree_compare[,3])
print(indegree_cohens_d)
plot(indegree_compare[,2], indegree_compare[,3])
cor(indegree_compare[,2], indegree_compare[,3], method = "pearson")
#indegree_chisquare <- chisq.test(indegree_compare[,2:3], correct = FALSE)
#indegree_fisher_test <- fisher.test(indegree_compare[,2:3])

#########################################################################################

#Using sub-sample of the data frame rows in order to avoid the issue of large samples

# Number of sub-samples (e.g., 100 sub-samples)
num_subsamples <- 100

# Size of each sub-sample (e.g., 50% of the original sample size)
subsample_size <- round(length(indegree_compare[,2]) * 0.5)

# Create vectors to store sub-sample t-test results
t_statistic <- rep(NA, num_subsamples)
p_value <- rep(NA, num_subsamples)

# Perform sub-sampling and t-test for each sub-sample
for (i in 1:num_subsamples) {
  # Randomly select subsample_size observations from each group
  sub_sample_table <- indegree_compare[sample(nrow(indegree_compare), subsample_size), ]
  
  # Perform t-test on the sub-sample
  indegree_result <- t.test(sub_sample_table[,2], sub_sample_table[,3], paired = TRUE)
  
  # Store t-test results
  t_statistic[i] <- indegree_result$statistic
  p_value[i] <- indegree_result$p.value
}

#########################################################################################

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
outdegree_ttest <- t.test(outdegree_compare[,2], outdegree_compare[,3], paired = TRUE)
print(outdegree_ttest)
outdegree_cohens_d <- effsize::cohen.d(outdegree_compare[,2], outdegree_compare[,3])
print(outdegree_cohens_d)
plot(outdegree_compare[,2], outdegree_compare[,3])
cor(outdegree_compare[,2], outdegree_compare[,3], method = "pearson")

#########################################################################################

#Using sub-sample of the data frame rows in order to avoid the issue of large samples

# Number of sub-samples (e.g., 100 sub-samples)
num_subsamples <- 100

# Size of each sub-sample (e.g., 50% of the original sample size)
subsample_size <- round(length(outdegree_compare[,2]) * 0.5)

# Create vectors to store sub-sample t-test results
t_statistic <- rep(NA, num_subsamples)
p_value <- rep(NA, num_subsamples)

# Perform sub-sampling and t-test for each sub-sample
for (i in 1:num_subsamples) {
  # Randomly select subsample_size observations from each group
  sub_sample_table <- outdegree_compare[sample(nrow(outdegree_compare), subsample_size), ]
  
  # Perform t-test on the sub-sample
  outdegree_result <- t.test(sub_sample_table[,2], sub_sample_table[,3], paired = TRUE)
  
  # Store t-test results
  t_statistic[i] <- outdegree_result$statistic
  p_value[i] <- outdegree_result$p.value
}

#########################################################################################

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
trophic_level_ttest <- t.test(trophic_level_compare[,2], trophic_level_compare[,3], paired = TRUE)
print(trophic_level_ttest)
trophic_level_cohens_d <- effsize::cohen.d(trophic_level_compare[,2], trophic_level_compare[,3])
print(trophic_level_cohens_d)
plot(trophic_level_compare[,2], trophic_level_compare[,3])
cor(trophic_level_compare[,2], trophic_level_compare[,3], method = "pearson")
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
closeness_ttest <- t.test(closeness_compare[,2], closeness_compare[,4], paired = TRUE)
print(closeness_ttest)
outdegree_cohens_d <- effsize::cohen.d(closeness_compare[,2], closeness_compare[,3])
print(outdegree_cohens_d)
plot(closeness_compare[,2], closeness_compare[,3])
cor(closeness_compare[,2], closeness_compare[,3], method = "pearson")
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
centrality_ttest <- t.test(centrality_compare[,2], centrality_compare[,3], paired = TRUE)
print(centrality_ttest)
centrality_cohens_d <- effsize::cohen.d(centrality_compare[,2], centrality_compare[,3])
print(centrality_cohens_d)
plot(centrality_compare[,2], centrality_compare[,3])
cor(centrality_compare[,2], centrality_compare[,3], method = "pearson")
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
ivi_ttest <- t.test(ivi_compare[,2], ivi_compare[,3], paired = TRUE)
print(ivi_ttest)
ivi_cohens_d <- effsize::cohen.d(ivi_compare[,2], ivi_compare[,3])
print(ivi_cohens_d)
plot(ivi_compare[,2], ivi_compare[,3])
cor(ivi_compare[,2], ivi_compare[,3], method = "pearson")
#

#######################################################################################
#    Where are areas significantly higher or lower in one metrics in T and NT?
#######################################################################################

names(indegree_compare)
names(outdegree_compare)
names(trophic_level_difference)
names(closeness_difference)
names(centrality_difference)
names(ivi_difference)
#
indegree_compare_2 <- data.frame(nt_indegree, t_indegree)
outdegree_compare_2 <- data.frame(nt_outdegree, t_outdegree)
trophic_level_compare_2 <- data.frame(nt_t_level, t_t_level)
closeness_compare_2 <- data.frame(nt_closeness, t_closeness)
centrality_compare_2 <- data.frame(nt_centrality, t_centrality)
ivi_compare_2 <- data.frame(nt_ivi, t_ivi)
#
indegree_difference <- indegree_compare_2[,4] - indegree_compare_2[,2]
outdegree_difference <- outdegree_compare_2[,4] - outdegree_compare_2[,2]
trophic_level_difference <- trophic_level_compare_2[,4] - trophic_level_compare_2[,2]
closeness_difference <- closeness_compare_2[,4] - closeness_compare_2[,2]
centrality_difference <- centrality_compare_2[,4] - centrality_compare_2[,2]
ivi_difference <- ivi_compare_2[,4] - ivi_compare_2[,2]

# Perform a t-test to assess the statistical significance
t_test_indeg <- t.test(indegree_compare_2[,2], indegree_compare_2[,4])
t_test_outdeg <- t.test(outdegree_compare_2[,2], outdegree_compare_2[,4])
t_test_tl <- t.test(trophic_level_compare_2[,2], trophic_level_compare_2[,4])
t_test_closeness <- t.test(closeness_compare_2[,2], closeness_compare_2[,4])
t_test_centrality <- t.test(centrality_compare_2[,2], centrality_compare_2[,4])
t_test_ivi <- t.test(ivi_compare_2[,2], ivi_compare_2[,4])

# Set the significance level
significance_level <- 0.05

# Identify areas where values in shapefile1 are significantly higher than shapefile2
indegree_shape <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/indegree_nt_spatial.shp")
indegree_shape <- indegree_shape[, -4]
indegree_shape$indegree_difference <- indegree_difference
#terra::plot(indegree_shape, "indegree_difference", col=rainbow(10))
indegree_shape$significantly_higher <- ifelse(t_test_indeg$p.value < significance_level & indegree_shape$indegree_difference > 0, "Yes", "No")
terra::plot(indegree_shape, "significantly_higher")
#######################
outdegree_shape <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/outdegree_nt_spatial.shp")
outdegree_shape <- outdegree_shape[, -4]
outdegree_shape$outdegree_difference <- outdegree_difference
#terra::plot(outdegree_shape, "outdegree_difference", col=rainbow(10))
outdegree_shape$significantly_higher <- ifelse(t_test_outdeg$p.value < significance_level & outdegree_shape$outdegree_difference > 0, "Yes", "No")
terra::plot(outdegree_shape, "significantly_higher")
#######################
tl_shape <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/tl_nt_spatial.shp")
tl_shape <- tl_shape[, -4]
tl_shape$trophic_level_difference <- trophic_level_difference
#terra::plot(tl_shape, "trophic_level_difference", col=rainbow(10))
tl_shape$significantly_higher <- ifelse(t_test_tl$p.value < significance_level & tl_shape$trophic_level_difference > 0, "Yes", "No")
terra::plot(tl_shape, "significantly_higher")
#######################
closeness_shape <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/closeness_nt_spatial.shp")
closeness_shape <- closeness_shape[, -4]
closeness_shape$closeness_difference <- closeness_difference
#terra::plot(closeness_shape, "closeness_difference", col=rainbow(10))
closeness_shape$significantly_higher <- ifelse(t_test_closeness$p.value < significance_level & closeness_shape$closeness_difference > 0, "Yes", "No")
terra::plot(closeness_shape, "significantly_higher")
#######################
centrality_shape <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/centrality_nt_spatial.shp")
centrality_shape <- closeness_shape[, -4]
centrality_shape$centrality_difference <- centrality_difference
#terra::plot(centrality_shape, "centrality_difference", col=rainbow(10))
centrality_shape$significantly_higher <- ifelse(t_test_centrality$p.value < significance_level & centrality_shape$centrality_difference > 0, "Yes", "No")
terra::plot(centrality_shape, "significantly_higher")
#######################
ivi_shape <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/ivi_nt_spatial_second_version.shp")
ivi_shape <- ivi_shape[, -4]
ivi_shape$ivi_difference <- ivi_difference
#terra::plot(ivi_shape, "ivi_difference", col=rainbow(10))
ivi_shape$significantly_higher <- ifelse(t_test_ivi$p.value < significance_level & ivi_shape$ivi_difference > 0, "Yes", "No")
terra::plot(ivi_shape, "significantly_higher")

#Save shapefiles
writeVector(indegree_shape, filename = "indegree_shape.shp")
writeVector(outdegree_shape, filename = "outdegree_shape.shp")
writeVector(tl_shape, filename = "tl_shape.shp")
writeVector(closeness_shape, filename = "closeness_shape.shp")
writeVector(centrality_shape, filename = "centrality_shape.shp")
writeVector(ivi_shape, filename = "ivi_shape.shp")

