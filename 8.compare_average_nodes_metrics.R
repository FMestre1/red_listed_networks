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

#Wilcoxon's Effects Size ####################################################### START
library("rstatix")
citation("rstatix")

#Indegree
df1_indegree_compare_3 <- data.frame(indegree_compare_2[,4],rep("T_indegree" , length(indegree_compare_2[,4])))
df2_indegree_compare_3 <- data.frame(indegree_compare_2[,2],rep("NT_indegree" , length(indegree_compare_2[,2])))
colnames(df1_indegree_compare_3) <- c("indegree", "group_T_NT")
colnames(df2_indegree_compare_3) <- c("indegree", "group_T_NT")
indegree_compare_3 <- rbind(df1_indegree_compare_3, df2_indegree_compare_3)
indegree_compare_3$group_T_NT <- as.factor(indegree_compare_3$group_T_NT)
#
wilcox_effsize(indegree ~ group_T_NT, paired = TRUE, data=indegree_compare_3)

#Outdegree
df1_outdegree_compare_3 <- data.frame(outdegree_compare_2[,4],rep("T_indegree" , length(outdegree_compare_2[,4])))
df2_outdegree_compare_3 <- data.frame(outdegree_compare_2[,2],rep("NT_indegree" , length(outdegree_compare_2[,2])))
colnames(df1_outdegree_compare_3) <- c("outdegree", "group_T_NT")
colnames(df2_outdegree_compare_3) <- c("outdegree", "group_T_NT")
outdegree_compare_3 <- rbind(df1_outdegree_compare_3, df2_outdegree_compare_3)
outdegree_compare_3$group_T_NT <- as.factor(outdegree_compare_3$group_T_NT)
#
wilcox_effsize(outdegree ~ group_T_NT, paired = TRUE, data=outdegree_compare_3)

#Trophic level
df1_trophic_level_compare_3 <- data.frame(trophic_level_compare_2[,4],rep("T_tlevel" , length(trophic_level_compare_2[,4])))
df2_trophic_level_compare_3 <- data.frame(trophic_level_compare_2[,2],rep("NT_tlevel" , length(trophic_level_compare_2[,2])))
colnames(df1_trophic_level_compare_3) <- c("tlevel", "group_T_NT")
colnames(df2_trophic_level_compare_3) <- c("tlevel", "group_T_NT")
trophic_level_compare_3 <- rbind(df1_trophic_level_compare_3, df2_trophic_level_compare_3)
trophic_level_compare_3$group_T_NT <- as.factor(trophic_level_compare_3$group_T_NT)
#
wilcox_effsize(tlevel ~ group_T_NT, paired = TRUE, data=trophic_level_compare_3)

#Closeness centrality
df1_closeness_compare_3 <- data.frame(closeness_compare_2[,4],rep("T_tlevel" , length(closeness_compare_2[,4])))
df2_closeness_compare_3 <- data.frame(closeness_compare_2[,2],rep("NT_tlevel" , length(closeness_compare_2[,2])))
colnames(df1_closeness_compare_3) <- c("closeness", "group_T_NT")
colnames(df2_closeness_compare_3) <- c("closeness", "group_T_NT")
closeness_compare_3 <- rbind(df1_closeness_compare_3, df2_closeness_compare_3)
closeness_compare_3$group_T_NT <- as.factor(closeness_compare_3$group_T_NT)
#
wilcox_effsize(closeness ~ group_T_NT, paired = TRUE, data=closeness_compare_3)

#Closeness centrality
df1_centrality_compare_3 <- data.frame(centrality_compare_2[,4],rep("T_tlevel" , length(centrality_compare_2[,4])))
df2_centrality_compare_3 <- data.frame(centrality_compare_2[,2],rep("NT_tlevel" , length(centrality_compare_2[,2])))
colnames(df1_centrality_compare_3) <- c("centrality", "group_T_NT")
colnames(df2_centrality_compare_3) <- c("centrality", "group_T_NT")
centrality_compare_3 <- rbind(df1_centrality_compare_3, df2_centrality_compare_3)
centrality_compare_3$group_T_NT <- as.factor(centrality_compare_3$group_T_NT)
#
wilcox_effsize(centrality ~ group_T_NT, paired = TRUE, data=centrality_compare_3)

#IVI
df1_ivi_compare_3 <- data.frame(ivi_compare_2[,4],rep("T_tlevel" , length(ivi_compare_2[,4])))
df2_ivi_compare_3 <- data.frame(ivi_compare_2[,2],rep("NT_tlevel" , length(ivi_compare_2[,2])))
colnames(df1_ivi_compare_3) <- c("ivi", "group_T_NT")
colnames(df2_ivi_compare_3) <- c("ivi", "group_T_NT")
ivi_compare_3 <- rbind(df1_ivi_compare_3, df2_ivi_compare_3)
ivi_compare_3$group_T_NT <- as.factor(ivi_compare_3$group_T_NT)
#
wilcox_effsize(ivi ~ group_T_NT, paired = TRUE, data=ivi_compare_3)

#Wilcoxon's Effects Size ####################################################### END

#Wilcoxon's Significance ####################################################### START

wilcox_indeg$p.value
wilcox_outdeg$p.value
wilcox_tl$p.value
wilcox_closeness$p.value
wilcox_centrality$p.value
wilcox_ivi$p.value

#Wilcoxon's Significance ####################################################### END

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

##########################################################################################
# Identify areas where values in shapefile1 are significantly higher than shapefile2
##########################################################################################

# Set the significance level
significance_level <- 0.05

indegree_shape <- terra::vect("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\indegree_nt_spatial.shp")
indegree_shape <- indegree_shape[, -4]
indegree_shape$indegree_difference <- indegree_difference
#terra::plot(indegree_shape, "indegree_difference", col=rainbow(10))
indegree_shape$significantly_higher <- ifelse(t_test_indeg$p.value < significance_level & indegree_shape$indegree_difference > 0, "Yes", "No")
terra::plot(indegree_shape, "significantly_higher")
#######################
outdegree_shape <- terra::vect("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\outdegree_nt_spatial.shp")
outdegree_shape <- outdegree_shape[, -4]
outdegree_shape$outdegree_difference <- outdegree_difference
#terra::plot(outdegree_shape, "outdegree_difference", col=rainbow(10))
outdegree_shape$significantly_higher <- ifelse(t_test_outdeg$p.value < significance_level & outdegree_shape$outdegree_difference > 0, "Yes", "No")
terra::plot(outdegree_shape, "significantly_higher")
#######################
tl_shape <- terra::vect("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\tl_nt_spatial.shp")
tl_shape <- tl_shape[, -4]
tl_shape$trophic_level_difference <- trophic_level_difference
#terra::plot(tl_shape, "trophic_level_difference", col=rainbow(10))
tl_shape$significantly_higher <- ifelse(t_test_tl$p.value < significance_level & tl_shape$trophic_level_difference > 0, "Yes", "No")
terra::plot(tl_shape, "significantly_higher")
#######################
closeness_shape <- terra::vect("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\closeness_nt_spatial.shp")
closeness_shape <- closeness_shape[, -4]
closeness_shape$closeness_difference <- closeness_difference
#terra::plot(closeness_shape, "closeness_difference", col=rainbow(10))
closeness_shape$significantly_higher <- ifelse(t_test_closeness$p.value < significance_level & closeness_shape$closeness_difference > 0, "Yes", "No")
terra::plot(closeness_shape, "significantly_higher")
#######################
centrality_shape <- terra::vect("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\centrality_nt_spatial.shp")
centrality_shape <- closeness_shape[, -4]
centrality_shape$centrality_difference <- centrality_difference
#terra::plot(centrality_shape, "centrality_difference", col=rainbow(10))
centrality_shape$significantly_higher <- ifelse(t_test_centrality$p.value < significance_level & centrality_shape$centrality_difference > 0, "Yes", "No")
terra::plot(centrality_shape, "significantly_higher")
#######################
ivi_shape <- terra::vect("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\ivi_nt_spatial_second_version.shp")
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
# USING THE RASTERS - Which is much memory consuming...
################################################################################

nt_indegree <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_indegree.tif")
t_indegree <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_indegree.tif")
#
nt_outdegree <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_outdegree.tif")
t_outdegree <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_outdegree.tif")
#
nt_t_level <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_trophic_level.tif")
t_t_level <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_trophic_level.tif")
#
nt_closeness <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_closeness.tif")
t_closeness <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_closeness.tif")
#
nt_centrality <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_centrality.tif")
t_centrality <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_centrality.tif")
#
nt_ivi <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_ivi.tif")
t_ivi <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_ivi.tif")

#####################################################################################

# COMPARING MAPS SPATIALLY
#FMestre
#28-08-2023

#Not used, as raster-based operations are much heavier in terms of memory 

indeg_diff <- t_indegree - nt_indegree
terra::writeRaster(indeg_diff, filename = "indeg_diff.tif")
plot(indeg_diff)
#
outdeg_diff <- t_outdegree - nt_outdegree
terra::writeRaster(outdeg_diff, filename = "outdeg_diff.tif")
plot(outdeg_diff)
#
trophic_level_diff <- t_t_level - nt_t_level
terra::writeRaster(trophic_level_diff, filename = "trophic_level_diff.tif")
plot(trophic_level_diff)
#
closeness_diff <- t_closeness - nt_closeness
terra::writeRaster(closeness_diff, filename = "closeness_diff.tif")
plot(closeness_diff)
#
centrality_diff <- t_centrality - nt_centrality
terra::writeRaster(centrality_diff, filename = "centrality_diff.tif")
plot(centrality_diff)
#
ivi_diff <- t_ivi - nt_ivi
terra::writeRaster(ivi_diff, filename = "ivi_diff.tif")
plot(ivi_diff)

#####################################################################################
