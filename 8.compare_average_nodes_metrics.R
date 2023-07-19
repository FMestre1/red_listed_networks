#FMestre
#10-07-2023

library(raster)
library(terra)
library(rasterVis)
library(cheddar)
library(effsize)
library(gridExtra)
library(viridis)
library(grid)

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
richness_raster <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\richness.tif")
richness_raster2 <- terra::resample(richness_raster, proportion)
#

richness <- terra::vect("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\shapefiles\\richness_2.shp")
proportion_vector <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/proportion_spatial.shp")
#
richness_proportion <- data.frame(merge(richness, proportion_vector))
#head(richness_proportion)
richness_proportion <- richness_proportion[,-c(2,3,4)]

#terra::plot(richness, "sp_rich")
#richness_2 <- terra::rasterize(richness, proportion, "Count")
#richness$sp_richnes
#rasterVis::levelplot(richness_2)

#Plot with rastervis
##
my.col.regions <- rev(terrain.colors(100))
######################
#max(nt_indegree) #max value: 25.28488 
#max(t_indegree) #max value: 61.000 
ind_min_max <- seq(0, 61, length.out = 100)
ind1 <- rasterVis::levelplot(nt_indegree, col.regions=my.col.regions, at=ind_min_max, main = "Not-threatened")
ind2 <- rasterVis::levelplot(t_indegree, col.regions=my.col.regions, at=ind_min_max, main = "Threatened")
ind_title <- textGrob("Average In-degree", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(ind1, ind2, ncol=2, top=ind_title)
######################
#max(nt_outdegree) #max value: 22.2991447
#max(t_outdegree) #max value: 20.4444447
out_min_max <- seq(0, 23, length.out = 100)
out1 <- rasterVis::levelplot(nt_outdegree, col.regions=my.col.regions, at=out_min_max, main = "Not-threatened")
out2 <- rasterVis::levelplot(t_outdegree, col.regions=my.col.regions, at=out_min_max, main = "Threatened")
out_title <- textGrob("Average Out-degree", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(out1, out2, ncol=2, top=out_title)
######################
#max(nt_t_level) #max value: 2
#max(t_t_level) #max value: 3 
tl_min_max <- seq(0, 3, length.out = 100)
tl1 <- rasterVis::levelplot(nt_t_level, col.regions=my.col.regions, at=tl_min_max, main = "Not-threatened")
tl2 <- rasterVis::levelplot(t_t_level, col.regions=my.col.regions, at=tl_min_max, main = "Threatened")
tl_title <- textGrob("Average Trophic Level", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(tl1, tl2, ncol=2, top=tl_title)
######################
#max(nt_closeness) #max value:
#max(t_closeness) #max value:
closeness_min_max <- seq(0, 1, length.out = 100)
cl1 <- rasterVis::levelplot(nt_closeness, col.regions=my.col.regions, at=closeness_min_max, main = "Not-threatened")
cl2 <- rasterVis::levelplot(t_closeness, col.regions=my.col.regions, at=closeness_min_max, main = "Threatened")
cl_title <- textGrob("Average Closeness Centrality", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(cl1, cl2, ncol=2, top=cl_title)
######################
#max(nt_centrality) #max value: 115.24037170
#max(t_centrality) #max value: 459.1808
centrality_min_max <- seq(0, 460, length.out = 100)
bt1 <- rasterVis::levelplot(nt_centrality, col.regions=my.col.regions, at=centrality_min_max, main = "Not-threatened")
bt2 <- rasterVis::levelplot(t_centrality, col.regions=my.col.regions, at=centrality_min_max, main = "Threatened")
bt_title <- textGrob("Average Betweenness Centrality", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(bt1, bt2, ncol=2, top=bt_title)
######################
#max(nt_ivi) #max value: 82.45447
#max(t_ivi) #max value: 100
ivi_min_max <- seq(0, 100, length.out = 100)
ivi1 <- rasterVis::levelplot(nt_ivi, col.regions=my.col.regions, at=ivi_min_max, main = "Not-threatened")
ivi2 <- rasterVis::levelplot(t_ivi, col.regions=my.col.regions, at=ivi_min_max, main = "Threatened")
ivi_title <- textGrob("Average IVI", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(ivi1, ivi2, ncol=2, top=ivi_title)
######################
rasterVis::levelplot(proportion, par.settings = rasterTheme(viridis_pal()(255)), main = "Proportion of threatened species")
rasterVis::levelplot(richness_raster2, par.settings = rasterTheme(viridis_pal()(255)), main = "Number of nodes in each network")

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

#frequency plots
names(indegree_compare)
par(mfrow=c(2, 2))
h_in_nt <- hist(indegree_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "In-degree", breaks = seq(0,65, by=2), ylim = c(0,40000), plot = FALSE)
h_in_t <- hist(indegree_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "In-degree",  breaks = seq(0,65, by=2), ylim = c(0,40000), plot = FALSE)
# Convert the counts to percentages
h_in_nt$counts <- h_in_nt$counts / sum(h_in_nt$counts) * 100
h_in_t$counts <- h_in_t$counts / sum(h_in_t$counts) * 100
#plot(h_in_nt, main = "Non-threatened", col = "darkgreen", xlab = "In-degree")
#plot(h_in_t, main = "Threatened", col = "darkred", xlab = "In-degree")
# Create the line plot
plot(h_in_nt$mids, h_in_nt$counts, type = "n", xlab = "In-degree", ylab = "Frequency (%)", main = "In-degree", ylim = c(0,40))
lines(h_in_nt$mids, h_in_nt$counts, lwd = 3, col = "darkgreen")
lines(h_in_t$mids, h_in_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)
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

#frequency plots
names(outdegree_compare)
par(mfrow=c(1, 2))
h_out_nt <- hist(outdegree_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "Out-degree", breaks = seq(0,65, by=2), ylim = c(0,40000))
h_out_t <- hist(outdegree_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "Out-degree",  breaks = seq(0,65, by=2), ylim = c(0,40000))
# Convert the counts to percentages
h_out_nt$counts <- h_out_nt$counts / sum(h_out_nt$counts) * 100
h_out_t$counts <- h_out_t$counts / sum(h_out_t$counts) * 100
# Create the line plot
plot(h_out_nt$mids, h_out_nt$counts, type = "n", xlab = "Out-degree", ylab = "Frequency (%)", main = "Out-degree", ylim = c(0,40))
lines(h_out_nt$mids, h_out_nt$counts, lwd = 3, col = "darkgreen")
lines(h_out_t$mids, h_out_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)
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

#frequency plots
names(trophic_level_compare)
par(mfrow=c(1, 2))
h_tl_nt <- hist(trophic_level_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "Trophic level", breaks = seq(0,4, by=0.25), ylim = c(0,100000))
h_tl_t <- hist(trophic_level_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "Trophic level",  breaks = seq(0,4, by=0.25), ylim = c(0,100000))
# Convert the counts to percentages
h_tl_nt$counts <- h_tl_nt$counts / sum(h_tl_nt$counts) * 100
h_tl_t$counts <- h_tl_t$counts / sum(h_tl_t$counts) * 100
# Create the line plot
plot(h_tl_nt$mids, h_tl_nt$counts, type = "n", xlab = "Trophic level", ylab = "Frequency (%)", main = "Trophic level", ylim = c(0,100))
lines(h_tl_nt$mids, h_tl_nt$counts, lwd = 3, col = "darkgreen")
lines(h_tl_t$mids, h_tl_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)
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

#frequency plots
names(closeness_compare)
par(mfrow=c(1, 2))
h_clos_nt <- hist(closeness_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "Closeness Centrality", breaks = seq(0,460, by=10), ylim = c(0,150000))
h_clos_t <- hist(closeness_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "Closeness Centrality",  breaks = seq(0,460, by=10), ylim = c(0,150000))
# Convert the counts to percentages
h_clos_nt$counts <- h_clos_nt$counts / sum(h_clos_nt$counts) * 100
h_clos_t$counts <- h_clos_t$counts / sum(h_clos_t$counts) * 100
# Create the line plot
plot(h_clos_nt$mids, h_clos_nt$counts, type = "n", xlab = "Closeness Centrality", ylab = "Frequency (%)", main = "Closeness Centrality", ylim = c(0,100))
lines(h_clos_nt$mids, h_clos_nt$counts, lwd = 3, col = "darkgreen")
lines(h_tl_t$mids, h_tl_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)
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

#frequency plots
names(centrality_compare)
par(mfrow=c(1, 2))
h_bet_nt <- hist(centrality_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "Betwenness Centrality", breaks = seq(0,460, by=10), ylim = c(0,150000))
h_bet_t <- hist(centrality_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "Betwenness Centrality",  breaks = seq(0,460, by=10), ylim = c(0,150000))
# Convert the counts to percentages
h_bet_nt$counts <- h_bet_nt$counts / sum(h_bet_nt$counts) * 100
h_bet_t$counts <- h_bet_t$counts / sum(h_bet_t$counts) * 100
# Create the line plot
plot(h_bet_nt$mids, h_bet_nt$counts, type = "n", xlab = "Betwenness Centrality", ylab = "Frequency (%)", main = "Betwenness Centrality", ylim = c(0,65))
lines(h_bet_nt$mids, h_bet_nt$counts, lwd = 3, col = "darkgreen")
lines(h_tl_t$mids, h_tl_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)
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
#frequency plots
names(ivi_compare)
par(mfrow=c(1, 2))
h_ivi_nt <- hist(ivi_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "IVI", breaks = seq(0,100, by=10), ylim = c(0,150000))
h_ivi_t <- hist(ivi_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "IVI",  breaks = seq(0,100, by=10), ylim = c(0,150000))
# Convert the counts to percentages
h_ivi_nt$counts <- h_ivi_nt$counts / sum(h_ivi_nt$counts) * 100
h_ivi_t$counts <- h_ivi_t$counts / sum(h_ivi_t$counts) * 100
# Create the line plot
plot(h_ivi_nt$mids, h_ivi_nt$counts, type = "n", xlab = "IVI", ylab = "Frequency (%)", main = "IVI", ylim = c(0,100))
lines(h_ivi_nt$mids, h_ivi_nt$counts, lwd = 3, col = "darkgreen")
lines(h_ivi_t$mids, h_ivi_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)
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

################################################################################
################################################################################

