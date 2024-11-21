################################################################################
#                             Evaluate differences           
################################################################################

#FMestre
#23-04-2024

#load packages
library(terra)
library(rasterVis)
library(diffeR)
library(RColorBrewer)

#Load Rasters
nt_ivi <- terra::rast("rasters_15JUL\\nt_ivi_15JUL.tif")
t_ivi <- terra::rast("rasters_15JUL\\t_ivi_15JUL.tif")
t_centrality <- terra::rast("rasters_15JUL\\t_centrality_15JUL.tif")
nt_centrality <- terra::rast("rasters_15JUL\\nt_centrality_15JUL.tif")
t_indegree <- terra::rast("rasters_15JUL\\t_indegree_15JUL.tif")
nt_indegree <- terra::rast("rasters_15JUL\\nt_indegree_15JUL.tif")
nt_outdegree <- terra::rast("rasters_15JUL\\nt_outdegree_15JUL.tif")
t_outdegree <- terra::rast("rasters_15JUL\\t_outdegree_15JUL.tif")
t_closeness <- terra::rast("rasters_15JUL\\t_closeness_15JUL.tif")
nt_closeness <- terra::rast("rasters_15JUL\\nt_closeness_15JUL.tif")
t_tl <- terra::rast("rasters_15JUL\\t_tl_15JUL.tif")
nt_tl <- terra::rast("rasters_15JUL\\nt_tl_15JUL.tif")
proportion_r <- terra::rast("rasters_15JUL\\proportion_r__15JUL.tif")
#
indeg_diff <- terra::rast("rasters_15JUL\\indegree_diff.tif")
outdeg_diff <- terra::rast("rasters_15JUL\\outdegree_diff.tif")
trophic_level_diff <- terra::rast("rasters_15JUL\\t_level_diff.tif")
closeness_diff <- terra::rast("rasters_15JUL\\closeness_diff.tif")
centrality_diff <- terra::rast("rasters_15JUL\\centrality_diff.tif")
ivi_diff <- terra::rast("rasters_15JUL\\ivi_diff.tif")

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
#                           Standartize difference maps
################################################################################

indegree_diff <- terra::rast("rasters_15JUL\\indegree_diff.tif")
outdegree_diff <- terra::rast("rasters_15JUL\\outdegree_diff.tif")
t_level_diff <- terra::rast("rasters_15JUL\\t_level_diff.tif")
closeness_diff <- terra::rast("rasters_15JUL\\closeness_diff.tif")
centrality_diff <- terra::rast("rasters_15JUL\\centrality_diff.tif")
ivi_diff <- terra::rast("rasters_15JUL\\ivi_diff.tif")

# Calculate mean and standard deviation for the raster layer
indegree_diff_mean <- as.numeric(terra::global(indegree_diff, "mean", na.rm = TRUE)[1])
indegree_diff_sd <- as.numeric(terra::global(indegree_diff, "sd", na.rm = TRUE)[1])
# Convert raster to z-scores
indegree_diff_z_score <- (indegree_diff - indegree_diff_mean) /indegree_diff_sd
plot(indegree_diff_z_score)

# Calculate mean and standard deviation for the raster layer
outdegree_diff_mean <- as.numeric(terra::global(outdegree_diff, "mean", na.rm = TRUE)[1])
outdegree_diff_sd <- as.numeric(terra::global(outdegree_diff, "sd", na.rm = TRUE)[1])
# Convert raster to z-scores
outdegree_diff_z_score <- (outdegree_diff - outdegree_diff_mean) /outdegree_diff_sd
plot(outdegree_diff_z_score)

# Calculate mean and standard deviation for the raster layer
t_level_diff_mean <- as.numeric(terra::global(t_level_diff, "mean", na.rm = TRUE)[1])
t_level_diff_sd <- as.numeric(terra::global(t_level_diff, "sd", na.rm = TRUE)[1])
# Convert raster to z-scores
t_level_diff_z_score <- (t_level_diff - t_level_diff_mean) /t_level_diff_sd
plot(t_level_diff_z_score)

# Calculate mean and standard deviation for the raster layer
closeness_diff_mean <- as.numeric(terra::global(closeness_diff, "mean", na.rm = TRUE)[1])
closeness_diff_sd <- as.numeric(terra::global(closeness_diff, "sd", na.rm = TRUE)[1])
# Convert raster to z-scores
closeness_diff_z_score <- (closeness_diff - closeness_diff_mean) /closeness_diff_sd
plot(closeness_diff_z_score)

# Calculate mean and standard deviation for the raster layer
centrality_diff_mean <- as.numeric(terra::global(centrality_diff, "mean", na.rm = TRUE)[1])
centrality_diff_sd <- as.numeric(terra::global(centrality_diff, "sd", na.rm = TRUE)[1])
# Convert raster to z-scores
centrality_diff_z_score <- (centrality_diff - centrality_diff_mean) /centrality_diff_sd
plot(centrality_diff_z_score)

# Calculate mean and standard deviation for the raster layer
ivi_diff_mean <- as.numeric(terra::global(ivi_diff, "mean", na.rm = TRUE)[1])
ivi_diff_sd <- as.numeric(terra::global(ivi_diff, "sd", na.rm = TRUE)[1])
# Convert raster to z-scores
ivi_diff_z_score <- (ivi_diff - ivi_diff_mean) /ivi_diff_sd
plot(ivi_diff_z_score)

#Save
terra::writeRaster(indegree_diff_z_score, "indegree_diff_z_score.tif")
terra::writeRaster(outdegree_diff_z_score, "outdegree_diff_z_score.tif")
terra::writeRaster(t_level_diff_z_score, "t_level_diff_z_score.tif")
terra::writeRaster(closeness_diff_z_score, "closeness_diff_z_score.tif")
terra::writeRaster(centrality_diff_z_score, "centrality_diff_z_score.tif")
terra::writeRaster(ivi_diff_z_score, "ivi_diff_z_score.tif")

################################################################################
# Plot
################################################################################

#Required function obtained from https://gist.github.com/johnbaums/306e4b7e69c87b1826db
diverge0 <- function(p, ramp) {
  # p: a trellis object resulting from rasterVis::levelplot
  # ramp: the name of an RColorBrewer palette (as character), a character 
  #       vector of colour names to interpolate, or a colorRampPalette.
  require(RColorBrewer)
  require(rasterVis)
  if(length(ramp)==1 && is.character(ramp) && ramp %in% 
     row.names(brewer.pal.info)) {
    ramp <- suppressWarnings(colorRampPalette(brewer.pal(11, ramp)))
  } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
    ramp <- colorRampPalette(ramp)
  } else if(!is.function(ramp)) 
    stop('ramp should be either the name of a RColorBrewer palette, ', 
         'a vector of colours to be interpolated, or a colorRampPalette.')
  rng <- range(p$legend[[1]]$args$key$at)
  s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
  i <- findInterval(rng[which.min(abs(rng))], s)
  zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
  p$legend[[1]]$args$key$at <- s[zlim]
  p[[grep('^legend', names(p))]][[1]]$args$key$col <- ramp(1000)[zlim[-length(zlim)]]
  p$panel.args.common$col.regions <- ramp(1000)[zlim[-length(zlim)]]
  p
}

#Load
indegree_diff_z_score <- terra::rast("indegree_diff_z_score.tif")
outdegree_diff_z_score <- terra::rast("outdegree_diff_z_score.tif")
t_level_diff_z_score <- terra::rast("t_level_diff_z_score.tif")
closeness_diff_z_score <- terra::rast("closeness_diff_z_score.tif")
centrality_diff_z_score <- terra::rast("centrality_diff_z_score.tif")
ivi_diff_z_score <- terra::rast("ivi_diff_z_score.tif")

#Study area borders
europe <- terra::vect("C:/Users/mestr/Documents/0. Artigos/IUCN_networks/shapefiles/europe_site.shp")
plot(europe)

#Create plots
indegree_diff_z_score_plot <- rasterVis::levelplot(indegree_diff_z_score, 
                                        contour = FALSE, 
                                        margin = NA,
                                        main = "In-degree",
                                        scales = list(draw = FALSE)
)

outdeg_diff_z_score_plot <- rasterVis::levelplot(outdegree_diff_z_score, 
                                         contour = FALSE, 
                                         margin = NA,
                                         main = "Out-degree",
                                         scales = list(draw = FALSE)
)

trophic_level_diff_z_score_plot <- rasterVis::levelplot(t_level_diff_z_score, 
                                                contour = FALSE, 
                                                margin = NA,
                                                main = "Trophic level",
                                                scales = list(draw = FALSE)
)

closeness_diff_z_score_plot <- rasterVis::levelplot(closeness_diff_z_score, 
                                            contour = FALSE, 
                                            margin = NA,
                                            main = "Closeness centrality",
                                            scales = list(draw = FALSE)
)

centrality_diff_z_score_plot <- rasterVis::levelplot(centrality_diff_z_score, 
                                             contour = FALSE, 
                                             margin = NA,
                                             main = "Betweenness centrality",
                                             scales = list(draw = FALSE)
)

ivi_diff_z_score_plot <- rasterVis::levelplot(ivi_diff_z_score, 
                                      contour = FALSE, 
                                      margin = NA,
                                      main = "IVI",
                                      scales = list(draw = FALSE)
)

#Create custom pallete
palette_fm <- colorRampPalette(colors = c("#67001F", "#CD2626", "#CDAD00", "#FFFFFF", "dodgerblue3", "#051A2E", "black"))(10)

#Create plots
pl1 <- diverge0(indegree_diff_z_score_plot, colorRampPalette(palette_fm))
pl2 <- diverge0(outdeg_diff_z_score_plot, colorRampPalette(palette_fm))
pl3 <- diverge0(trophic_level_diff_z_score_plot, colorRampPalette(palette_fm))
pl4 <- diverge0(closeness_diff_z_score_plot, colorRampPalette(palette_fm))
pl5 <- diverge0(centrality_diff_z_score_plot, colorRampPalette(palette_fm))
pl6 <- diverge0(ivi_diff_z_score_plot, colorRampPalette(palette_fm))

#Compose pairs of maps
gridExtra::grid.arrange(pl3, pl6, ncol=2)
gridExtra::grid.arrange(pl1, pl4, ncol=2)
gridExtra::grid.arrange(pl2, pl5, ncol=2)


#All plots together
#gridExtra::grid.arrange(pl3, pl6, 
#                        pl1, pl4,
#                        pl2, pl5,
#                        ncol=2)



