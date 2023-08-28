#FMestre
#10-07-2023

################################################################################
#                                  Plot maps
################################################################################

################################################################################
#   Where are areas significantly higher or lower in one metrics in T and NT?
################################################################################

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
t_test_indeg <- t.test(indegree_compare_2[,4], indegree_compare_2[,2], paired = TRUE, alternative = "two.sided")
t_test_outdeg <- t.test(outdegree_compare_2[,4], outdegree_compare_2[,2], paired = TRUE, alternative = "two.sided")
t_test_tl <- t.test(trophic_level_compare_2[,4], trophic_level_compare_2[,2], paired = TRUE, alternative = "two.sided")
t_test_closeness <- t.test(closeness_compare_2[,4], closeness_compare_2[,2], paired = TRUE, alternative = "two.sided")
t_test_centrality <- t.test(centrality_compare_2[,4], centrality_compare_2[,2], paired = TRUE, alternative = "two.sided")
t_test_ivi <- t.test(ivi_compare_2[,4], ivi_compare_2[,2], paired = TRUE, alternative = "two.sided")

#Use wilcoxon test - to try and avoid the effect of the large sample size leading to significant differences
  #relate with richness and proportion

richness_proportion_indegree <- merge(richness_proportion, indegree_compare_2)
richness_proportion_indegree <- richness_proportion_indegree[,-5]
head(richness_proportion_indegree)
wilcox.test(richness_proportion_indegree[,2], richness_proportion_indegree[,4], paired = TRUE)
wilcox.test(richness_proportion_indegree[,2], richness_proportion_indegree[,5], paired = TRUE)
wilcox.test(richness_proportion_indegree[,3], richness_proportion_indegree[,4], paired = TRUE)
wilcox.test(richness_proportion_indegree[,3], richness_proportion_indegree[,5], paired = TRUE)
#
richness_proportion_outdegree <- merge(richness_proportion, outdegree_compare_2)
richness_proportion_outdegree <- richness_proportion_outdegree[,-5]
head(richness_proportion_outdegree)
wilcox.test(richness_proportion_outdegree[,2], richness_proportion_outdegree[,4], paired = TRUE)
wilcox.test(richness_proportion_outdegree[,2], richness_proportion_outdegree[,5], paired = TRUE)
wilcox.test(richness_proportion_outdegree[,3], richness_proportion_outdegree[,4], paired = TRUE)
wilcox.test(richness_proportion_outdegree[,3], richness_proportion_outdegree[,5], paired = TRUE)
#
richness_proportion_tl <- merge(richness_proportion, trophic_level_compare_2)
richness_proportion_tl <- richness_proportion_tl[,-5]
head(richness_proportion_tl)
wilcox.test(richness_proportion_tl[,2], richness_proportion_tl[,4], paired = TRUE)
wilcox.test(richness_proportion_tl[,2], richness_proportion_tl[,5], paired = TRUE)
wilcox.test(richness_proportion_tl[,3], richness_proportion_tl[,4], paired = TRUE)
wilcox.test(richness_proportion_tl[,3], richness_proportion_tl[,5], paired = TRUE)
#
richness_proportion_closeness <- merge(richness_proportion, closeness_compare_2)
richness_proportion_closeness <- richness_proportion_closeness[,-5]
head(richness_proportion_closeness)
wilcox.test(richness_proportion_closeness[,2], richness_proportion_closeness[,4], paired = TRUE)
wilcox.test(richness_proportion_closeness[,2], richness_proportion_closeness[,5], paired = TRUE)
wilcox.test(richness_proportion_closeness[,3], richness_proportion_closeness[,4], paired = TRUE)
wilcox.test(richness_proportion_closeness[,3], richness_proportion_closeness[,5], paired = TRUE)
#
richness_proportion_centrality <- merge(richness_proportion, centrality_compare_2)
richness_proportion_centrality <- richness_proportion_centrality[,-5]
head(richness_proportion_centrality)
wilcox.test(richness_proportion_centrality[,2], richness_proportion_centrality[,4], paired = TRUE)
wilcox.test(richness_proportion_centrality[,2], richness_proportion_centrality[,5], paired = TRUE)
wilcox.test(richness_proportion_centrality[,3], richness_proportion_centrality[,4], paired = TRUE)
wilcox.test(richness_proportion_centrality[,3], richness_proportion_centrality[,5], paired = TRUE)
#
richness_proportion_ivi <- merge(richness_proportion, ivi_compare_2)
richness_proportion_ivi <- richness_proportion_ivi[,-5]
head(richness_proportion_ivi)
wilcox.test(richness_proportion_ivi[,2], richness_proportion_ivi[,4], paired = TRUE)
wilcox.test(richness_proportion_ivi[,2], richness_proportion_ivi[,5], paired = TRUE)
wilcox.test(richness_proportion_ivi[,3], richness_proportion_ivi[,4], paired = TRUE)
wilcox.test(richness_proportion_ivi[,3], richness_proportion_ivi[,5], paired = TRUE)

#######

wilcox_indeg <- wilcox.test(indegree_compare_2[,4], indegree_compare_2[,2], paired = TRUE)
wilcox_outdeg <- wilcox.test(outdegree_compare_2[,4], outdegree_compare_2[,2], paired = TRUE)
wilcox_tl <- wilcox.test(trophic_level_compare_2[,4], trophic_level_compare_2[,2], paired = TRUE)
wilcox_closeness <- wilcox.test(closeness_compare_2[,4], closeness_compare_2[,2], paired = TRUE)
wilcox_centrality <- wilcox.test(centrality_compare_2[,4], centrality_compare_2[,2], paired = TRUE)
wilcox_ivi <- wilcox.test(ivi_compare_2[,4], ivi_compare_2[,2], paired = TRUE)
#
wilcox_indeg$p.value
wilcox_outdeg$p.value
wilcox_tl$p.value
wilcox_closeness$p.value
wilcox_centrality$p.value
wilcox_ivi$p.value

###

#CORRELATION

round(cor(indegree_compare_2[,4], indegree_compare_2[,2], method = "spearman", use = "complete"), 3)
round(cor(outdegree_compare_2[,4], outdegree_compare_2[,2], method = "spearman", use = "complete"), 3)
round(cor(trophic_level_compare_2[,4], trophic_level_compare_2[,2], method = "spearman", use = "complete"), 3)
round(cor(closeness_compare_2[,4], closeness_compare_2[,2], method = "spearman", use = "complete"), 3)
round(cor(centrality_compare_2[,4], centrality_compare_2[,2], method = "spearman", use = "complete"), 3)
round(cor(ivi_compare_2[,4], ivi_compare_2[,2], method = "spearman", use = "complete"), 3)

###

#BOXPLOTS

indegree_compare_3_NT <- data.frame(indegree_compare_2[,2], rep("NT", length = length(indegree_compare_2[,2])))
indegree_compare_3_T <- data.frame(indegree_compare_2[,4], rep("T", length = length(indegree_compare_2[,4])))
names(indegree_compare_3_NT) <- c("indeg", "group")
names(indegree_compare_3_T) <- c("indeg", "group")
indegree_compare_3 <- rbind(indegree_compare_3_T, indegree_compare_3_NT)
ggplot(indegree_compare_3, aes(x=group, y=indeg)) + geom_boxplot(notch = TRUE, outlier.colour="grey")
ggplot(indegree_compare_3, aes(x=group, y=indeg)) +
  geom_violin(aes(fill = group)) +
  guides(fill = "none") +
  labs(title = "Average In-degree") +
  scale_fill_manual(values = c("darkgreen", "darkred"))
#####
outdegree_compare_3_NT <- data.frame(outdegree_compare_2[,2], rep("NT", length = length(outdegree_compare_2[,2])))
outdegree_compare_3_T <- data.frame(outdegree_compare_2[,4], rep("T", length = length(outdegree_compare_2[,4])))
names(outdegree_compare_3_NT) <- c("outdeg", "group")
names(outdegree_compare_3_T) <- c("outdeg", "group")
outdegree_compare_3 <- rbind(outdegree_compare_3_T, outdegree_compare_3_NT)
#ggplot(outdegree_compare_3, aes(x=group, y=outdeg)) + geom_boxplot(notch = TRUE, outlier.colour="grey")
ggplot(outdegree_compare_3, aes(x=group, y=outdeg)) +
  geom_violin(aes(fill = group)) +
  guides(fill = "none") +
  labs(title = "Average Out-degree") +
  scale_fill_manual(values = c("darkgreen", "darkred"))
#####
trophic_level_compare_3_NT <- data.frame(trophic_level_compare_2[,2], rep("NT", length = length(trophic_level_compare_2[,2])))
trophic_level_compare_3_T <- data.frame(trophic_level_compare_2[,4], rep("T", length = length(trophic_level_compare_2[,4])))
names(trophic_level_compare_3_NT) <- c("tl", "group")
names(trophic_level_compare_3_T) <- c("tl", "group")
trophic_level_compare_3 <- rbind(trophic_level_compare_3_T, trophic_level_compare_3_NT)
#ggplot(trophic_level_compare_3, aes(x=group, y=tl)) + geom_boxplot(notch = TRUE, outlier.colour="grey")
ggplot(trophic_level_compare_3, aes(x=group, y=tl)) +
  geom_violin(aes(fill = group)) +
  guides(fill = "none") +
  labs(title = "Average Trophic Level") +
  scale_fill_manual(values = c("darkgreen", "darkred"))
#####
closeness_compare_3_NT <- data.frame(closeness_compare_2[,2], rep("NT", length = length(closeness_compare_2[,2])))
closeness_compare_3_T <- data.frame(closeness_compare_2[,4], rep("T", length = length(closeness_compare_2[,4])))
names(closeness_compare_3_NT) <- c("closeness", "group")
names(closeness_compare_3_T) <- c("closeness", "group")
closeness_compare_3 <- rbind(closeness_compare_3_T, closeness_compare_3_NT)
#ggplot(closeness_compare_3, aes(x=group, y=closeness)) + geom_boxplot(notch = TRUE, outlier.colour="grey")
ggplot(closeness_compare_3, aes(x=group, y=closeness)) +
  geom_violin(aes(fill = group)) +
  guides(fill = "none") +
  labs(title = "Average Closeness Centrality") +
  scale_fill_manual(values = c("darkgreen", "darkred"))
#####
centrality_compare_3_NT <- data.frame(centrality_compare_2[,2], rep("NT", length = length(centrality_compare_2[,2])))
centrality_compare_3_T <- data.frame(centrality_compare_2[,4], rep("T", length = length(centrality_compare_2[,4])))
names(centrality_compare_3_NT) <- c("centrality", "group")
names(centrality_compare_3_T) <- c("centrality", "group")
centrality_compare_3 <- rbind(centrality_compare_3_T, centrality_compare_3_NT)
#ggplot(centrality_compare_3, aes(x=group, y=centrality)) + geom_boxplot(notch = TRUE, outlier.colour="grey")
ggplot(centrality_compare_3, aes(x=group, y=centrality)) +
  geom_violin(aes(fill = group)) +
  guides(fill = "none") +
  labs(title = "Average Betweeness Centrality") +
  scale_fill_manual(values = c("darkgreen", "darkred"))
#####
ivi_compare_3_NT <- data.frame(ivi_compare_2[,2], rep("NT", length = length(ivi_compare_2[,2])))
ivi_compare_3_T <- data.frame(ivi_compare_2[,4], rep("T", length = length(ivi_compare_2[,4])))
names(ivi_compare_3_NT) <- c("ivi", "group")
names(ivi_compare_3_T) <- c("ivi", "group")
ivi_compare_3 <- rbind(ivi_compare_3_T, ivi_compare_3_NT)
#ggplot(ivi_compare_3, aes(x=group, y=ivi)) + geom_boxplot(notch = TRUE, outlier.colour="grey")
ggplot(ivi_compare_3, aes(x=group, y=ivi)) +
  geom_violin(aes(fill = group)) +
  guides(fill = "none") +
  labs(title = "Average IVI") +
  scale_fill_manual(values = c("darkgreen", "darkred"))
###

#Cohen's d
library(effsize)
indegree_cohens_d_2 <- effsize::cohen.d(indegree_compare_2[complete.cases(indegree_compare_2),][,4],
                                        indegree_compare_2[complete.cases(indegree_compare_2),][,2], paired = TRUE)
outdegree_cohens_d_2 <- effsize::cohen.d(outdegree_compare_2[complete.cases(outdegree_compare_2),][,4],
                                         outdegree_compare_2[complete.cases(outdegree_compare_2),][,2], paired = TRUE)
tl_cohens_d_2 <- effsize::cohen.d(trophic_level_compare_2[complete.cases(trophic_level_compare_2),][,4],
                                  trophic_level_compare_2[complete.cases(trophic_level_compare_2),][,2], paired = TRUE)
closeness_cohens_d_2 <- effsize::cohen.d(closeness_compare_2[complete.cases(closeness_compare_2),][,4],
                                         closeness_compare_2[complete.cases(closeness_compare_2),][,2], paired = TRUE)
centrality_cohens_d_2 <- effsize::cohen.d(centrality_compare_2[complete.cases(centrality_compare_2),][,4],
                                          centrality_compare_2[complete.cases(centrality_compare_2),][,2], paired = TRUE)
ivi_cohens_d_2 <- effsize::cohen.d(ivi_compare_2[complete.cases(ivi_compare_2),][,4],
                                   ivi_compare_2[complete.cases(ivi_compare_2),][,2], paired = TRUE)

round(indegree_cohens_d_2$estimate, 3)
round(outdegree_cohens_d_2$estimate, 3)
round(tl_cohens_d_2$estimate, 3)
round(closeness_cohens_d_2$estimate, 3)
round(centrality_cohens_d_2$estimate, 3)
round(ivi_cohens_d_2$estimate, 3)

# Set the significance level
significance_level <- 0.05

##########################################################################################
# Identify areas where values in shapefile1 are significantly higher than shapefile2
##########################################################################################

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



################################################################################
#
################################################################################

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

###
#rm(t_indegree, nt_indegree,
#   t_outdegree, nt_outdegree,
#   t_t_level, nt_t_level,
#   t_closeness, nt_closeness,
#   t_centrality, nt_centrality,
#   t_ivi, nt_ivi)
###

#indeg_ttest <- t.test(t_indegree, nt_indegree)
#outdeg_ttest <- t.test(t_outdegree, nt_outdegree)
#tl_ttest <- t.test(t_t_level, nt_t_level)
#closeness_ttest <- t.test(t_closeness, nt_closeness)
#centrality_ttest <- t.test(t_centrality, nt_centrality)
#ivi_ttest <- t.test(t_ivi, nt_ivi)

mask_vect <- terra::vect("C:\\Users\\FMest\\Documents\\github\\red_listed_networks\\europeRaster_poly.shp")

# Add the centroid coordinates
centroids <- terra::centroids(mask_vect, TRUE)

# Add the centroid coordinates to the vector
mask_vect <- terra::cbind2(mask_vect, as.data.frame(terra::crds(centroids)))
longitudes <- seq(min(mask_vect$x), max(mask_vect$x), by = 800000)

#plot(terra::subset(mask_vect, mask_vect$x > longitudes[i] & mask_vect$x  < longitudes[i+1]))

i = 3

mask1 <- terra::subset(mask_vect, mask_vect$x > longitudes[i] & mask_vect$x  < longitudes[i+1])
#plot(mask1)
ext1 <- terra::ext(mask1)

#plot(terra::mask(t_indegree, mask1))
mp1 <- terra::mask(t_indegree, mask1)
mp2 <- terra::mask(nt_indegree, mask1)
#
terra::ext(mp1) <- ext1
terra::ext(mp2) <- ext1
#
plot(mp1)
plot(mp2)
#
tt1 <- t.test(mp1, mp2)
#
rm(mp1, mp2)
