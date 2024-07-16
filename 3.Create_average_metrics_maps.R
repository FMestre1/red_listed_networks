################################################################################
#                                   Spatialize
################################################################################
#FMestre
#23-10-2023

#Load packages
library(terra)
library(viridis)
library(rasterVis)
library(gridExtra)
library(grid)

#Load objects
load("igraph_list_02SET23.RData")
load("cheddar_list_02SET23.RData")
load("all_species_status_body_mass_amph_13_20OUT.RData")
#
sp_richness <- terra::vect("shapes_20OUT23\\sp_richness_23OUT.shp")
ivi_nt_spatial <- terra::vect("shapes_20OUT23\\ivi_nt_spatial_second_version_20OUT.shp")
ivi_t_spatial <- terra::vect("shapes_20OUT23\\ivi_t_spatial_second_version_20OUT.shp")
centrality_nt_spatial <- terra::vect("shapes_20OUT23\\centrality_nt_spatial_20OUT.shp")
centrality_t_spatial <- terra::vect("shapes_20OUT23\\centrality_t_spatial_20OUT.shp")
indegree_t_spatial <- terra::vect("shapes_20OUT23\\indegree_t_spatial_20OUT.shp")
indegree_nt_spatial <- terra::vect("shapes_20OUT23\\indegree_nt_spatial_20OUT.shp")
outdegree_t_spatial <- terra::vect("shapes_20OUT23\\outdegree_t_spatial_20OUT.shp")
outdegree_nt_spatial <- terra::vect("shapes_20OUT23\\outdegree_nt_spatial_20OUT.shp")
closeness_t_spatial <- terra::vect("shapes_20OUT23\\closeness_t_spatial_20OUT23.shp")
closeness_nt_spatial <- terra::vect("shapes_20OUT23\\closeness_nt_spatial_20OUT.shp")
tl_nt_spatial <- terra::vect("shapes_20OUT23\\tl_nt_spatial_20OUT.shp")
tl_t_spatial <- terra::vect("shapes_20OUT23\\tl_t_spatial_20OUT.shp")
proportion_spatial <- terra::vect("shapes_20OUT23\\proportion_spatial_20OUT.shp")
#
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
#
load("indegree_compare_23OUT.RData")
load("outdegree_compare_23OUT.RData")

#grid <- terra::vect("C:/Users/FMest/Documents/0. Artigos/IUCN_networks/shapefiles/Europa_10km/europe_10km.shp")
#crs(grid)

europe <- terra::vect("C:/Users/asus/Documents/0. Artigos/IUCN_networks/shapefiles/Europe/Europe.shp")
#terra::crs(europe)

europe_coastline_borders <- terra::aggregate(europe, dissolve = TRUE, "merge")
plot(europe_coastline_borders)

#Save shapefile
writeVector(europe_coastline_borders, 
            filename = "europe_coastline_borders.shp", 
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=FALSE, 
            options="ENCODING=UTF-8"
            )

europeRaster <- terra::rast(x="C:\\Users\\asus\\Documents\\github\\red_listed_networks\\mask10k-20230214T144658Z-001\\mask10k\\reference_grid_10km.img")
cells_info <- foreign::read.dbf(file = "C:/Users/asus/Documents/github/red_listed_networks/mask10k-20230214T144658Z-001/mask10k/reference_grid_10km.img.vat.dbf")

#To vector
europeRaster_poly <- terra::as.polygons(europeRaster, values = TRUE, extent=FALSE)
europeRaster_poly <- terra::merge(europeRaster_poly, cells_info)
europeRaster_poly_wgs84 <- terra::project(europeRaster_poly, europe)
#To WGS84
europeRaster_poly_wgs84_coords <- crds(europeRaster_poly_wgs84, df=TRUE)
europeRaster_poly_wgs84_coords <- data.frame(europeRaster_poly_wgs84, europeRaster_poly_wgs84_coords)

#Write vector
#writeVector(europeRaster_poly, filename ="europeRaster_poly.shp", overwrite=TRUE, filetype = "ESRI Shapefile")
#writeVector(europeRaster_poly_wgs84, filename ="europeRaster_poly_wgs84.shp", overwrite=TRUE, filetype = "ESRI Shapefile")

################################################################################

################################################################################
#                               Species Richnesss          
################################################################################

#fw_list_with_status
#cheddar_list
#igraph_list

#length(fw_list_with_status)
#length(cheddar_list)
#length(igraph_list)

#Get number of species, checking if the three sources of information agree

richness <- data.frame(matrix(ncol = 3, nrow = length(igraph_list)))
names(richness) <- c("grid", "spe_igraph", "spe_cheddar")

for(i in 1:length(igraph_list)){
  
  richness$grid[i] <- names(cheddar_list)[i]
  richness$spe_igraph[i] <- length(igraph::V(igraph_list[[i]]))
  richness$spe_cheddar[i] <- cheddar::NumberOfNodes(cheddar_list[[i]])
  
  message(i)

}

sp_richness <- merge(europeRaster_poly, richness, by.x = "PageName", by.y = "grid")
#save(sp_richness, file = "sp_richness.RData")

writeVector(sp_richness, 
            filename = "shape_15JUL24\\sp_richness_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

################################################################################

fw_list_with_status <- fw_list

#Add IUCN status
for(i in 1:length(fw_list_with_status)){
  
  fw3 <- fw_list_with_status[[i]]
  fw3$SP_NAME <- stringr::str_replace(fw3$SP_NAME, "_", " ")
  sp_fw3 <- fw3$SP_NAME
  
  if(any(new_species_iucn$scientificName %in% sp_fw3))
    {
    sp_fw3_redList <- all_species_status_body_mass_amph_15[all_species_status_body_mass_amph_15$species %in% sp_fw3 | all_species_status_body_mass_amph_15$synonym %in% sp_fw3,]
    fw4 <- merge(fw3, sp_fw3_redList, by.x = "SP_NAME", by.y = "species", all.x = TRUE)
    fw_list_with_status[[i]] <- fw4[,c(-2)]
  }


message(i)

}

#fw_list_with_status[[2]]

#Save & Load
#save(fw_list_with_status, file = "fw_list_with_status_15jul2024.RData")
#load("fw_list_with_status_15jul2024.RData")

#IVI - NOT THREATENED ##########################################################

ivi_nt <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  ivi_fw <- fw_list_with_status[[i]]
  fw_nt_ivi <- ivi_fw[ivi_fw$corrected_agreg_ts == "not_threatened",]$ivi
  if(length(fw_nt_ivi)!=0) ivi_nt[i] <- mean(fw_nt_ivi, na.rm=TRUE)
  
  message(i)
  
}

ivi_nt <- data.frame(names(fw_list_with_status), ivi_nt)
names(ivi_nt) <- c("grid", "ivi")
#head(ivi_nt)
#hist(ivi_nt$ivi)

ivi_nt_spatial <- terra::merge(europeRaster_poly, ivi_nt, by.x = "PageName", by.y = "grid")
#save(ivi_nt_spatial, file = "ivi_nt_spatial_15JUL.Rdata")

terra::writeVector(ivi_nt_spatial, 
            filename = "shape_15JUL24\\ivi_nt_spatial_second_version_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#ivi_nt_spatial <- terra::vect("shape_15JUL24\\ivi_nt_spatial_second_version_15JUL.shp")

#IVI - THREATENED ##############################################################

ivi_t <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  ivi_fw2 <- fw_list_with_status[[i]]
  fw_t_ivi2 <- ivi_fw2[ivi_fw2$corrected_agreg_ts == "threatened",]$ivi
  if(length(fw_t_ivi2)!=0) ivi_t[i] <- mean(fw_t_ivi2, na.rm=TRUE)
  
  message(i)
  
}

ivi_t <- data.frame(names(fw_list_with_status), ivi_t)
names(ivi_t) <- c("grid", "ivi")
#head(ivi_t)

ivi_t_spatial <- merge(europeRaster_poly, ivi_t, by.x = "PageName", by.y = "grid")
#save(ivi_t_spatial, file = "ivi_t_spatial_15JUL.Rdata")

writeVector(ivi_t_spatial, 
            filename = "shape_15JUL24\\ivi_t_spatial_second_version_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#ivi_t_spatial <- terra::vect("shape_15JUL24\\ivi_t_spatial_second_version_15JUL.shp")

#CENTRALITY - NON-THREATENED #######################################################

centrality_nt <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  cent_nt_fw2 <- fw_list_with_status[[i]]
  fw_nt_cent2 <- cent_nt_fw2[cent_nt_fw2$corrected_agreg_ts == "not_threatened",]$centrality
  if(length(fw_nt_cent2)!=0) centrality_nt[i] <- mean(fw_nt_cent2, na.rm=TRUE)
  
  message(i)
  
}

centrality_nt <- data.frame(names(fw_list_with_status), centrality_nt)
names(centrality_nt) <- c("grid", "centrality")
#head(centrality_nt)

centrality_nt_spatial <- merge(europeRaster_poly, centrality_nt, by.x = "PageName", by.y = "grid")
#save(centrality_nt_spatial, file = "centrality_nt_spatial_15JUL.Rdata")

writeVector(centrality_nt_spatial, 
            filename = "shape_15JUL24\\centrality_nt_spatial_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#centrality_nt_spatial <- terra::vect("shape_15JUL24\\centrality_nt_spatial_15JUL.shp")

#CENTRALITY - THREATENED #######################################################

centrality_t <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  cent_fw2 <- fw_list_with_status[[i]]
  fw_t_cent2 <- cent_fw2[cent_fw2$corrected_agreg_ts == "threatened",]$centrality
  if(length(fw_t_cent2)!=0) centrality_t[i] <- mean(fw_t_cent2, na.rm=TRUE)
  
  message(i)
  
}

centrality_t <- data.frame(names(fw_list_with_status), centrality_t)
names(centrality_t) <- c("grid", "centrality")
#head(centrality_t)

centrality_t_spatial <- merge(europeRaster_poly, centrality_t, by.x = "PageName", by.y = "grid")
#save(centrality_t_spatial, file = "centrality_t_spatial_15JUL.Rdata")

writeVector(centrality_t_spatial, 
            filename = "shape_15JUL24\\centrality_t_spatial_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#centrality_t_spatial <- terra::vect("shape_15JUL24\\centrality_t_spatial_15JUL.shp")

#IN-DEGREE - THREATENED ########################################################

indegree_t <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  indeg_fw2 <- fw_list_with_status[[i]]
  fw_t_indeg2 <- indeg_fw2[indeg_fw2$corrected_agreg_ts == "threatened",]$indegree
  if(length(fw_t_indeg2)!=0) indegree_t[i] <- mean(fw_t_indeg2, na.rm=TRUE)
  
  message(i)
  
}

indegree_t <- data.frame(names(fw_list_with_status), indegree_t)
names(indegree_t) <- c("grid", "indegree")
#head(indegree_t)

indegree_t_spatial <- merge(europeRaster_poly, indegree_t, by.x = "PageName", by.y = "grid")
#save(indegree_t_spatial, file = "indegree_t_spatial_15JUL.Rdata")

writeVector(indegree_t_spatial, 
            filename = "shape_15JUL24\\indegree_t_spatial_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#indegree_t_spatial <- terra::vect("shape_15JUL24\\indegree_t_spatial_15JUL.shp")

#IN-DEGREE - NON-THREATENED ####################################################

indegree_nt <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  indeg_nt_fw2 <- fw_list_with_status[[i]]
  fw_nt_indeg2 <- indeg_nt_fw2[indeg_nt_fw2$corrected_agreg_ts == "not_threatened",]$indegree
  if(length(fw_nt_indeg2)!=0) indegree_nt[i] <- mean(fw_nt_indeg2, na.rm=TRUE)
  
  message(i)
  
}

indegree_nt <- data.frame(names(fw_list_with_status), indegree_nt)
names(indegree_nt) <- c("grid", "indegree")
#head(indegree_nt)

indegree_nt_spatial <- merge(europeRaster_poly, indegree_nt, by.x = "PageName", by.y = "grid")
#save(indegree_nt_spatial, file = "indegree_nt_spatial_15JUL.Rdata")

writeVector(indegree_nt_spatial, 
            filename = "shape_15JUL24\\indegree_nt_spatial_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#indegree_nt_spatial <- terra::vect("shape_15JUL24\\indegree_nt_spatial_15JUL.shp")

#OUT-DEGREE - THREATENED ########################################################

outegree_t <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  outdeg_fw2 <- fw_list_with_status[[i]]
  fw_t_outdeg2 <- outdeg_fw2[outdeg_fw2$corrected_agreg_ts == "threatened",]$outdegree
  if(length(fw_t_outdeg2)!=0) outegree_t[i] <- mean(fw_t_outdeg2, na.rm=TRUE)
  
  message(i)
  
}

outdegree_t <- data.frame(names(fw_list_with_status), outegree_t)
names(outdegree_t) <- c("grid", "outdegree")
#head(outdegree_t)

outdegree_t_spatial <- merge(europeRaster_poly, outdegree_t, by.x = "PageName", by.y = "grid")
#save(outdegree_t_spatial, file = "outdegree_t_spatial_15JUL.Rdata")

writeVector(outdegree_t_spatial, 
            filename = "shape_15JUL24\\outdegree_t_spatial_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#outdegree_t_spatial <- terra::vect("shape_15JUL24\\outdegree_t_spatial_15JUL.shp")

#OUT-DEGREE - NON-THREATENED ####################################################

outegree_nt <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  outdeg_fw2 <- fw_list_with_status[[i]]
  fw_nt_outdeg2 <- outdeg_fw2[outdeg_fw2$corrected_agreg_ts == "not_threatened",]$outdegree
  if(length(fw_nt_outdeg2)!=0) outegree_nt[i] <- mean(fw_nt_outdeg2, na.rm=TRUE)
  
  message(i)
  
}

outdegree_nt <- data.frame(names(fw_list_with_status), outegree_nt)
names(outdegree_nt) <- c("grid", "outdegree")
#head(outdegree_nt)

outdegree_nt_spatial <- merge(europeRaster_poly, outdegree_nt, by.x = "PageName", by.y = "grid")
#save(outdegree_nt_spatial, file = "outdegree_nt_spatial_15JUL.Rdata")

writeVector(outdegree_nt_spatial, 
            filename = "shape_15JUL24\\outdegree_nt_spatial_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#outdegree_nt_spatial <- terra::vect("shape_15JUL24\\outdegree_nt_spatial_15JUL.shp")

#CLOSENESS - THREATENED ########################################################

closeness_t <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  close_t_fw2 <- fw_list_with_status[[i]]
  fw_t_close2 <- close_t_fw2[close_t_fw2$corrected_agreg_ts == "threatened",]$closeness
  if(length(fw_t_close2)!=0) closeness_t[i] <- mean(fw_t_close2, na.rm=TRUE)
  
  message(i)
  
}

closeness_t <- data.frame(names(fw_list_with_status), closeness_t)
names(closeness_t) <- c("grid", "closeness")
#head(closeness_t)

closeness_t_spatial <- merge(europeRaster_poly, closeness_t, by.x = "PageName", by.y = "grid")
#save(closeness_t_spatial, file = "closeness_t_spatial_15JUL.Rdata")

writeVector(closeness_t_spatial, 
            filename = "shape_15JUL24\\closeness_t_spatial_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#closeness_t_spatial <- terra::vect("shape_15JUL24\\closeness_t_spatial_15JUL.shp")

#CLOSENESS - NON-THREATENED ####################################################

closeness_nt <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  close_nt_fw2 <- fw_list_with_status[[i]]
  fw_nt_close2 <- close_nt_fw2[close_nt_fw2$corrected_agreg_ts == "not_threatened",]$closeness
  if(length(fw_nt_close2)!=0) closeness_nt[i] <- mean(fw_nt_close2, na.rm=TRUE)
  
  message(i)
  
}

closeness_nt <- data.frame(names(fw_list_with_status), closeness_nt)
names(closeness_nt) <- c("grid", "closeness")
#head(closeness_nt)

closeness_nt_spatial <- merge(europeRaster_poly, closeness_nt, by.x = "PageName", by.y = "grid")
#save(closeness_nt_spatial, file = "closeness_nt_spatial_15JUL.Rdata")

writeVector(closeness_nt_spatial, 
            filename = "shape_15JUL24\\closeness_nt_spatial_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#closeness_nt_spatial <- terra::vect("shape_15JUL24\\closeness_nt_spatial_15JUL.shp")

#TROPHIC LEVEL - NON-THREATENED ################################################

TL_nt <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  tl_nt_fw2 <- fw_list_with_status[[i]]
  fw_nt_TL2 <- tl_nt_fw2[tl_nt_fw2$corrected_agreg_ts == "not_threatened",]$TL
  if(length(fw_nt_TL2)!=0) TL_nt[i] <- mean(fw_nt_TL2, na.rm=TRUE)
  
  message(i)
  
}

TL_nt <- data.frame(names(fw_list_with_status), TL_nt)
names(TL_nt) <- c("grid", "trophic_level")
#head(TL_nt)

tl_nt_spatial <- merge(europeRaster_poly, TL_nt, by.x = "PageName", by.y = "grid")
#save(tl_nt_spatial, file = "tl_nt_spatial_15JUL.Rdata")

writeVector(tl_nt_spatial, 
            filename = "shape_15JUL24\\tl_nt_spatial_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#tl_nt_spatial <- terra::vect("shape_15JUL24\\tl_nt_spatial_15JUL.shp")

#TROPHIC LEVEL - THREATENED ####################################################

TL_t <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  tl_t_fw2 <- fw_list_with_status[[i]]
  fw_t_TL2 <- tl_t_fw2[tl_t_fw2$corrected_agreg_ts == "threatened",]$TL
  if(length(fw_t_TL2)!=0) TL_t[i] <- mean(fw_t_TL2, na.rm=TRUE)
  
  message(i)
  
}

TL_t <- data.frame(names(fw_list_with_status), TL_t)
names(TL_t) <- c("grid", "trophic_level")
#head(TL_t)

tl_t_spatial <- merge(europeRaster_poly, TL_t, by.x = "PageName", by.y = "grid")
#save(tl_t_spatial, file = "tl_t_spatial_15JUL.Rdata")

writeVector(tl_t_spatial, 
            filename = "shape_15JUL24\\tl_t_spatial_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#tl_t_spatial <- terra::vect("shape_15JUL24\\tl_t_spatial_15JUL.shp")

#PROPORTION OF THREATENED SPECIES ##############################################

propT <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  propT_fw2 <- fw_list_with_status[[i]]
  prop <- propT_fw2$corrected_agreg_ts
  n_t <- length(prop[prop=="threatened"])
  n_total <- length(prop)
  propT[i] <- n_t/n_total
  
  #if(length(fw_t_TL2)!=0) TL_t[i] <- mean(fw_t_TL2)
  
  message(i)
  
}

propT <- data.frame(names(fw_list_with_status), propT)
names(propT) <- c("grid", "proportion")
#head(propT)

proportion_spatial <- merge(europeRaster_poly, propT, by.x = "PageName", by.y = "grid")
#save(proportion_spatial, file = "proportion_spatial_15JUL.Rdata")

writeVector(proportion_spatial, 
            filename = "shape_15JUL24\\proportion_spatial_15JUL.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#proportion_spatial <- terra::vect("shape_15JUL24\\proportion_spatial_15JUL.shp")
