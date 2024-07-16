################################################################################
#                          Plot maps and frequency           
################################################################################

#Load packages
library(rasterVis)
library(grid)
library(gridExtra)
library(scales)
library(grDevices)

#Plot with rastervis
##
my.col.regions <- rev(terrain.colors(100))
######################
#summary(nt_indegree) #max value: 22.71 
#summary(t_indegree) #max value: 52.000 
ind_min_max <- seq(0, 55, length.out = 100)
ind1 <- rasterVis::levelplot(nt_indegree, col.regions=my.col.regions, at=ind_min_max, main = "Not-threatened")
ind2 <- rasterVis::levelplot(t_indegree, col.regions=my.col.regions, at=ind_min_max, main = "Threatened")
ind_title <- textGrob("Average In-degree", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(ind1, ind2, ncol=2, top=ind_title)
######################
#summary(nt_outdegree) #max value: 22.22
#summary(t_outdegree) #max value: 32.00
out_min_max <- seq(0, 35, length.out = 100)
out1 <- rasterVis::levelplot(nt_outdegree, col.regions=my.col.regions, at=out_min_max, main = "Not-threatened")
out2 <- rasterVis::levelplot(t_outdegree, col.regions=my.col.regions, at=out_min_max, main = "Threatened")
out_title <- textGrob("Average Out-degree", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(out1, out2, ncol=2, top=out_title)
######################
#summary(nt_tl) #max value: 2
#summary(t_tl) #max value: 3 
tl_min_max <- seq(0, 3, length.out = 100)
tl1 <- rasterVis::levelplot(nt_tl, col.regions=my.col.regions, at=tl_min_max, main = "Not-threatened")
tl2 <- rasterVis::levelplot(t_tl, col.regions=my.col.regions, at=tl_min_max, main = "Threatened")
tl_title <- textGrob("Average Trophic Level", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(tl1, tl2, ncol=2, top=tl_title)
######################
#summary(nt_closeness) #max value: 1
#summary(t_closeness) #max value: 0.75
closeness_min_max <- seq(0, 1, length.out = 100)
cl1 <- rasterVis::levelplot(nt_closeness, col.regions=my.col.regions, at=closeness_min_max, main = "Not-threatened")
cl2 <- rasterVis::levelplot(t_closeness, col.regions=my.col.regions, at=closeness_min_max, main = "Threatened")
cl_title <- textGrob("Average Closeness Centrality", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(cl1, cl2, ncol=2, top=cl_title)
######################
#summary(nt_centrality) #max value: 111.55
#summary(t_centrality) #max value: 442.03
centrality_min_max <- seq(0, 450, length.out = 100)
bt1 <- rasterVis::levelplot(nt_centrality, col.regions=my.col.regions, at=centrality_min_max, main = "Not-threatened")
bt2 <- rasterVis::levelplot(t_centrality, col.regions=my.col.regions, at=centrality_min_max, main = "Threatened")
bt_title <- textGrob("Average Betweenness Centrality", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(bt1, bt2, ncol=2, top=bt_title)
######################
#summary(nt_ivi) #max value: 68.70
#summary(t_ivi) #max value: 100.00
ivi_min_max <- seq(0, 100, length.out = 100)
ivi1 <- rasterVis::levelplot(nt_ivi, col.regions=my.col.regions, at=ivi_min_max, main = "Not-threatened")
ivi2 <- rasterVis::levelplot(t_ivi, col.regions=my.col.regions, at=ivi_min_max, main = "Threatened")
ivi_title <- textGrob("Average IVI", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(ivi1, ivi2, ncol=2, top=ivi_title)

################################################################################
#                             Supplementary figures
################################################################################

#Figure S1 - Number of species in each local network.
#species_rich_r <- terra::rast("rasters_15JUL\\species_rich_r_15JUL.tif")
#summary(species_rich_r)
png(filename = "Fig_S1.png",width = 500, height = 500)
rasterVis::levelplot(species_rich_r, par.settings = BuRdTheme(), main = "Number of species per network")
dev.off()

#Figure S2 - Proportion of threatened species in each local network.
#proportion_r <- terra::rast("rasters_15JUL\\proportion_r__15JUL.tif")
#summary(proportion_r)
png(filename = "Fig_S2.png",width = 500, height = 500)
rasterVis::levelplot(proportion_r, par.settings = BuRdTheme(), main = "Proportion of threatened species")
dev.off()


