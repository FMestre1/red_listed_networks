################################################################################
#                          Compute and  plot differences           
################################################################################

#FMestre
#23-10-2023

#load packages
library(terra)
library(rasterVis)
library(gridExtra)

#Load Rasters
#nt_ivi <- terra::rast("rasters_21OUT\\nt_ivi_20OUT.tif")
#t_ivi <- terra::rast("rasters_21OUT\\t_ivi_20OUT.tif")
#t_centrality <- terra::rast("rasters_21OUT\\t_centrality_20OUT.tif")
#nt_centrality <- terra::rast("rasters_21OUT\\nt_centrality_20OUT.tif")
#t_indegree <- terra::rast("rasters_21OUT\\t_indegree_20OUT.tif")
#nt_indegree <- terra::rast("rasters_21OUT\\nt_indegree_20OUT.tif")
#nt_outdegree <- terra::rast("rasters_21OUT\\nt_outdegree_20OUT.tif")
#t_outdegree <- terra::rast("rasters_21OUT\\t_outdegree_20OUT.tif")
#t_closeness <- terra::rast("rasters_21OUT\\t_closeness_20OUT.tif")
#nt_closeness <- terra::rast("rasters_21OUT\\nt_closeness_20OUT.tif")
#t_tl <- terra::rast("rasters_21OUT\\t_tl_20OUT.tif")
#nt_tl <- terra::rast("rasters_21OUT\\nt_tl_20OUT.tif")
#proportion_r <- terra::rast("rasters_21OUT\\proportion_r_20OUT.tif")
#
#indeg_diff <- terra::rast("rasters_21OUT\\indeg_diff.tif")
#outdeg_diff <- terra::rast("rasters_21OUT\\outdeg_diff.tif")
#trophic_level_diff <- terra::rast("rasters_21OUT\\trophic_level_diff.tif")
#closeness_diff <- terra::rast("rasters_21OUT\\closeness_diff.tif")
#centrality_diff <- terra::rast("rasters_21OUT\\centrality_diff.tif")
#ivi_diff <- terra::rast("rasters_21OUT\\ivi_diff.tif")

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

#Create plots
indeg_diff_plot <- rasterVis::levelplot(indegree_diff, 
                     contour = FALSE, 
                     margin = NA,
                     main = "In-degree",
                     scales = list(draw = FALSE)
                     )

outdeg_diff_plot <- rasterVis::levelplot(outdegree_diff, 
                                        contour = FALSE, 
                                        margin = NA,
                                        main = "Out-degree",
                                        scales = list(draw = FALSE)
)

trophic_level_diff_plot <- rasterVis::levelplot(t_level_diff, 
                                         contour = FALSE, 
                                         margin = NA,
                                         main = "Trophic level",
                                         scales = list(draw = FALSE)
)

closeness_diff_plot <- rasterVis::levelplot(closeness_diff, 
                                            contour = FALSE, 
                                            margin = NA,
                                            main = "Closeness centrality",
                                            scales = list(draw = FALSE)
)

centrality_diff_plot <- rasterVis::levelplot(centrality_diff, 
                                            contour = FALSE, 
                                            margin = NA,
                                            main = "Betweenness centrality",
                                            scales = list(draw = FALSE)
)

ivi_diff_plot <- rasterVis::levelplot(ivi_diff, 
                                      contour = FALSE, 
                                      margin = NA,
                                      main = "IVI",
                                      scales = list(draw = FALSE)
                                      )

pl1 <- diverge0(indeg_diff_plot, 'Spectral')
pl2 <- diverge0(outdeg_diff_plot, 'Spectral')
pl3 <- diverge0(trophic_level_diff_plot, 'Spectral')
pl4 <- diverge0(closeness_diff_plot, 'Spectral')
pl5 <- diverge0(centrality_diff_plot, 'Spectral')
pl6 <- diverge0(ivi_diff_plot, 'Spectral')

#Compose pairs of maps
grid.arrange(pl3, pl6, ncol=2)
grid.arrange(pl1, pl4, ncol=2)
grid.arrange(pl2, pl5, ncol=2)
