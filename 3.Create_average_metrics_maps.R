################################################################################
#                                   Spatialize
################################################################################
#FMestre
#09-02-2023

load("C:\\Users\\asus\\Desktop\\igraph_list_02SET23.RData")
load("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\networks_SET23\\cheddar_list_02SET23.RData")
load("all_species_status_body_mass_amph_13_20OUT.RData")
#
terra::vect("shapes_20OUT23\\sp_richness_23OUT.shp")
terra::vect("shapes_20OUT23\\ivi_nt_spatial_second_version_20OUT.shp")
terra::vect("shapes_20OUT23\\ivi_t_spatial_second_version_20OUT.shp")
terra::vect("shapes_20OUT23\\centrality_nt_spatial_20OUT.shp")
terra::vect("shapes_20OUT23\\centrality_t_spatial_20OUT.shp")
terra::vect("shapes_20OUT23\\indegree_t_spatial_20OUT.shp")
terra::vect("shapes_20OUT23\\indegree_nt_spatial_20OUT.shp")
terra::vect("shapes_20OUT23\\outdegree_t_spatial_20OUT.shp")
terra::vect("shapes_20OUT23\\outdegree_nt_spatial_20OUT.shp")
terra::vect("shapes_20OUT23\\closeness_t_spatial_20OUT23.shp")
terra::vect("shapes_20OUT23\\closeness_nt_spatial_20OUT.shp")
terra::vect("shapes_20OUT23\\tl_nt_spatial_20OUT.shp")
terra::vect("shapes_20OUT23\\tl_t_spatial_20OUT.shp")
terra::vect("shapes_20OUT23\\proportion_spatial_20OUT.shp")
#
terra::rast("rasters_21OUT\\nt_ivi_20OUT.tif")
terra::rast("rasters_21OUT\\t_ivi_20OUT.tif")
terra::rast("rasters_21OUT\\t_centrality_20OUT.tif")
terra::rast("rasters_21OUT\\nt_centrality_20OUT.tif")
terra::rast("rasters_21OUT\\t_indegree_20OUT.tif")
terra::rast("rasters_21OUT\\nt_indegree_20OUT.tif")
terra::rast("rasters_21OUT\\nt_outdegree_20OUT.tif")
terra::rast("rasters_21OUT\\t_outdegree_20OUT.tif")
terra::rast("rasters_21OUT\\t_closeness_20OUT.tif")
terra::rast("rasters_21OUT\\nt_closeness_20OUT.tif")
terra::rast("rasters_21OUT\\t_tl_20OUT.tif")
terra::rast("rasters_21OUT\\nt_tl_20OUT.tif")
terra::rast("rasters_21OUT\\proportion_r_20OUT.tif")

#Load packages
library(terra)
library(viridis)
library(rasterVis)
library(gridExtra)

#grid <- terra::vect("C:/Users/FMest/Documents/0. Artigos/IUCN_networks/shapefiles/Europa_10km/europe_10km.shp")
#crs(grid)

europe <- terra::vect("C:/Users/asus/Documents/0. Artigos/IUCN_networks/shapefiles/Europe/Europe.shp")
#terra::crs(europe)

europe_coastline_borders <- aggregate(europe, dissolve = TRUE)
#plot(europe_coastline_borders)

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

length(fw_list_with_status)
length(cheddar_list)
length(igraph_list)

#Get number of species, checking if the three sources of information agree

richness_23OUT <- data.frame(matrix(ncol = 3, nrow = length(igraph_list)))
names(richness_23OUT) <- c("grid", "spe_igraph", "spe_cheddar")

for(i in 1:length(igraph_list)){
  
  richness_23OUT$grid[i] <- names(cheddar_list)[i]
  richness_23OUT$spe_igraph[i] <- length(igraph::V(igraph_list[[i]]))
  richness_23OUT$spe_cheddar[i] <- cheddar::NumberOfNodes(cheddar_list[[i]])
  
  message(i)

}

sp_richness_23OUT <- merge(europeRaster_poly, richness_23OUT, by.x = "PageName", by.y = "grid")

writeVector(sp_richness_23OUT, 
            filename = "shapes_20OUT23\\sp_richness_23OUT.shp",
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
    sp_fw3_redList <- all_species_status_body_mass_amph_13[all_species_status_body_mass_amph_13$species %in% sp_fw3 | all_species_status_body_mass_amph_13$synonym %in% sp_fw3,]
    fw4 <- merge(fw3, sp_fw3_redList, by.x = "SP_NAME", by.y = "species", all.x = TRUE)
    fw_list_with_status[[i]] <- fw4
  }


message(i)

}

#fw_list_with_status[[2]]

#save(fw_list_with_status, file = "fw_list_with_status_20OUT.RData")
#load("fw_list_with_status_20OUT.RData")

#IVI - NOT THREATENED ##########################################################

ivi_nt <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  ivi_fw <- fw_list_with_status[[i]]
  fw_nt_ivi <- ivi_fw[ivi_fw$agreg_ts == "not_threatened",]$ivi
  if(length(fw_nt_ivi)!=0) ivi_nt[i] <- mean(fw_nt_ivi, na.rm=TRUE)
  
  message(i)
  
}

ivi_nt <- data.frame(names(fw_list_with_status), ivi_nt)
names(ivi_nt) <- c("grid", "ivi")
#head(ivi_nt)
#hist(ivi_nt$ivi)

ivi_nt_spatial <- merge(europeRaster_poly, ivi_nt, by.x = "PageName", by.y = "grid")
#save(ivi_nt_spatial, file = "ivi_nt_spatial_20OUT.Rdata")

writeVector(ivi_nt_spatial, 
            filename = "shapes_20OUT23\\ivi_nt_spatial_second_version_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#ivi_nt_spatial <- terra::vect("shapes_20OUT23\\ivi_nt_spatial_second_version_20OUT.shp")

#IVI - THREATENED ##############################################################

ivi_t <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  ivi_fw2 <- fw_list_with_status[[i]]
  fw_t_ivi2 <- ivi_fw2[ivi_fw2$agreg_ts == "threatened",]$ivi
  if(length(fw_t_ivi2)!=0) ivi_t[i] <- mean(fw_t_ivi2, na.rm=TRUE)
  
  message(i)
  
}

ivi_t <- data.frame(names(fw_list_with_status), ivi_t)
names(ivi_t) <- c("grid", "ivi")
#head(ivi_t)

ivi_t_spatial <- merge(europeRaster_poly, ivi_t, by.x = "PageName", by.y = "grid")
#save(ivi_t_spatial, file = "ivi_t_spatial_20OUT.Rdata")

writeVector(ivi_t_spatial, 
            filename = "shapes_20OUT23\\ivi_t_spatial_second_version_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

ivi_t_spatial <- terra::vect("ivi_t_spatial_second_version_20OUT.shp")

#CENTRALITY - NON-THREATENED #######################################################

centrality_nt <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  cent_nt_fw2 <- fw_list_with_status[[i]]
  fw_nt_cent2 <- cent_nt_fw2[cent_nt_fw2$agreg_ts == "not_threatened",]$centrality
  if(length(fw_nt_cent2)!=0) centrality_nt[i] <- mean(fw_nt_cent2, na.rm=TRUE)
  
  message(i)
  
}

centrality_nt <- data.frame(names(fw_list_with_status), centrality_nt)
names(centrality_nt) <- c("grid", "centrality")
#head(centrality_nt)

centrality_nt_spatial <- merge(europeRaster_poly, centrality_nt, by.x = "PageName", by.y = "grid")
#save(centrality_nt_spatial, file = "centrality_nt_spatial_20OUT23.Rdata")

writeVector(centrality_nt_spatial, 
            filename = "shapes_20OUT23\\centrality_nt_spatial_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#centrality_nt_spatial <- terra::vect("centrality_nt_spatial_20OUT.shp")

#CENTRALITY - THREATENED #######################################################

centrality_t <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  cent_fw2 <- fw_list_with_status[[i]]
  fw_t_cent2 <- cent_fw2[cent_fw2$agreg_ts == "threatened",]$centrality
  if(length(fw_t_cent2)!=0) centrality_t[i] <- mean(fw_t_cent2, na.rm=TRUE)
  
  message(i)
  
}

centrality_t <- data.frame(names(fw_list_with_status), centrality_t)
names(centrality_t) <- c("grid", "centrality")
#head(centrality_t)

centrality_t_spatial <- merge(europeRaster_poly, centrality_t, by.x = "PageName", by.y = "grid")
#save(centrality_t_spatial, file = "centrality_t_spatial_20OUT.Rdata")

writeVector(centrality_t_spatial, 
            filename = "shapes_20OUT23\\centrality_t_spatial_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#centrality_t_spatial <- terra::vect("centrality_t_spatial_20OUT.shp")

#IN-DEGREE - THREATENED ########################################################

indegree_t <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  indeg_fw2 <- fw_list_with_status[[i]]
  fw_t_indeg2 <- indeg_fw2[indeg_fw2$agreg_ts == "threatened",]$indegree
  if(length(fw_t_indeg2)!=0) indegree_t[i] <- mean(fw_t_indeg2, na.rm=TRUE)
  
  message(i)
  
}

indegree_t <- data.frame(names(fw_list_with_status), indegree_t)
names(indegree_t) <- c("grid", "indegree")
#head(indegree_t)

indegree_t_spatial <- merge(europeRaster_poly, indegree_t, by.x = "PageName", by.y = "grid")
#save(indegree_t_spatial, file = "indegree_t_spatial_20OUT.Rdata")

writeVector(indegree_t_spatial, 
            filename = "shapes_20OUT23\\indegree_t_spatial_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#indegree_t_spatial <- terra::vect("indegree_t_spatial_20OUT.shp")

#IN-DEGREE - NON-THREATENED ####################################################

indegree_nt <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  indeg_nt_fw2 <- fw_list_with_status[[i]]
  fw_nt_indeg2 <- indeg_nt_fw2[indeg_nt_fw2$agreg_ts == "not_threatened",]$indegree
  if(length(fw_nt_indeg2)!=0) indegree_nt[i] <- mean(fw_nt_indeg2, na.rm=TRUE)
  
  message(i)
  
}

indegree_nt <- data.frame(names(fw_list_with_status), indegree_nt)
names(indegree_nt) <- c("grid", "indegree")
#head(indegree_nt)

indegree_nt_spatial <- merge(europeRaster_poly, indegree_nt, by.x = "PageName", by.y = "grid")
#save(indegree_nt_spatial, file = "indegree_nt_spatial_20OUT.Rdata")

writeVector(indegree_nt_spatial, 
            filename = "shapes_20OUT23\\indegree_nt_spatial_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#indegree_nt_spatial <- terra::vect("indegree_nt_spatial_20OUT.shp")

#OUT-DEGREE - THREATENED ########################################################

outegree_t <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  outdeg_fw2 <- fw_list_with_status[[i]]
  fw_t_outdeg2 <- outdeg_fw2[outdeg_fw2$agreg_ts == "threatened",]$outdegree
  if(length(fw_t_outdeg2)!=0) outegree_t[i] <- mean(fw_t_outdeg2, na.rm=TRUE)
  
  message(i)
  
}

outdegree_t <- data.frame(names(fw_list_with_status), outegree_t)
names(outdegree_t) <- c("grid", "outdegree")
#head(outdegree_t)

outdegree_t_spatial <- merge(europeRaster_poly, outdegree_t, by.x = "PageName", by.y = "grid")
#save(outdegree_t_spatial, file = "outdegree_t_spatial_20OUT.Rdata")

writeVector(outdegree_t_spatial, 
            filename = "shapes_20OUT23\\outdegree_t_spatial_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#outdegree_t_spatial <- terra::vect("outdegree_t_spatial_20OUT.shp")

#OUT-DEGREE - NON-THREATENED ####################################################

outegree_nt <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  outdeg_fw2 <- fw_list_with_status[[i]]
  fw_nt_outdeg2 <- outdeg_fw2[outdeg_fw2$agreg_ts == "not_threatened",]$outdegree
  if(length(fw_nt_outdeg2)!=0) outegree_nt[i] <- mean(fw_nt_outdeg2, na.rm=TRUE)
  
  message(i)
  
}

outdegree_nt <- data.frame(names(fw_list_with_status), outegree_nt)
names(outdegree_nt) <- c("grid", "outdegree")
#head(outdegree_nt)

outdegree_nt_spatial <- merge(europeRaster_poly, outdegree_nt, by.x = "PageName", by.y = "grid")
#save(outdegree_nt_spatial, file = "outdegree_nt_spatial_20OUT.Rdata")

writeVector(outdegree_nt_spatial, 
            filename = "shapes_20OUT23\\outdegree_nt_spatial_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#outdegree_nt_spatial <- terra::vect("outdegree_nt_spatial_20OUT.shp")

#CLOSENESS - THREATENED ########################################################

closeness_t <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  close_t_fw2 <- fw_list_with_status[[i]]
  fw_t_close2 <- close_t_fw2[close_t_fw2$agreg_ts == "threatened",]$closeness
  if(length(fw_t_close2)!=0) closeness_t[i] <- mean(fw_t_close2, na.rm=TRUE)
  
  message(i)
  
}

closeness_t <- data.frame(names(fw_list_with_status), closeness_t)
names(closeness_t) <- c("grid", "closeness")
#head(closeness_t)

closeness_t_spatial <- merge(europeRaster_poly, closeness_t, by.x = "PageName", by.y = "grid")
#save(closeness_t_spatial, file = "closeness_t_spatial_20OUT.Rdata")

writeVector(closeness_t_spatial, 
            filename = "shapes_20OUT23\\closeness_t_spatial_20OUT23.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#closeness_t_spatial <- terra::vect("closeness_t_spatial_20OUT23.shp")

#CLOSENESS - NON-THREATENED ####################################################

closeness_nt <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  close_nt_fw2 <- fw_list_with_status[[i]]
  fw_nt_close2 <- close_nt_fw2[close_nt_fw2$agreg_ts == "not_threatened",]$closeness
  if(length(fw_nt_close2)!=0) closeness_nt[i] <- mean(fw_nt_close2, na.rm=TRUE)
  
  message(i)
  
}

closeness_nt <- data.frame(names(fw_list_with_status), closeness_nt)
names(closeness_nt) <- c("grid", "closeness")
#head(closeness_nt)

closeness_nt_spatial <- merge(europeRaster_poly, closeness_nt, by.x = "PageName", by.y = "grid")
#save(closeness_nt_spatial, file = "closeness_nt_spatial_20OUT.Rdata")

writeVector(closeness_nt_spatial, 
            filename = "shapes_20OUT23\\closeness_nt_spatial_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#closeness_nt_spatial <- terra::vect("closeness_nt_spatial_20OUT.shp")

#TROPHIC LEVEL - NON-THREATENED ################################################

TL_nt <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  tl_nt_fw2 <- fw_list_with_status[[i]]
  fw_nt_TL2 <- tl_nt_fw2[tl_nt_fw2$agreg_ts == "not_threatened",]$TL
  if(length(fw_nt_TL2)!=0) TL_nt[i] <- mean(fw_nt_TL2, na.rm=TRUE)
  
  message(i)
  
}

TL_nt <- data.frame(names(fw_list_with_status), TL_nt)
names(TL_nt) <- c("grid", "trophic_level")
#head(TL_nt)

tl_nt_spatial <- merge(europeRaster_poly, TL_nt, by.x = "PageName", by.y = "grid")
#save(tl_nt_spatial, file = "tl_nt_spatial_20OUT.Rdata")

writeVector(tl_nt_spatial, 
            filename = "shapes_20OUT23\\tl_nt_spatial_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#tl_nt_spatial <- terra::vect("tl_nt_spatial_20OUT.shp")

#TROPHIC LEVEL - THREATENED ####################################################

TL_t <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  tl_t_fw2 <- fw_list_with_status[[i]]
  fw_t_TL2 <- tl_t_fw2[tl_t_fw2$agreg_ts == "threatened",]$TL
  if(length(fw_t_TL2)!=0) TL_t[i] <- mean(fw_t_TL2, na.rm=TRUE)
  
  message(i)
  
}

TL_t <- data.frame(names(fw_list_with_status), TL_t)
names(TL_t) <- c("grid", "trophic_level")
#head(TL_t)

tl_t_spatial <- merge(europeRaster_poly, TL_t, by.x = "PageName", by.y = "grid")
#save(tl_t_spatial, file = "tl_t_spatial_20OUT.Rdata")

writeVector(tl_t_spatial, 
            filename = "shapes_20OUT23\\tl_t_spatial_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#tl_t_spatial <- terra::vect("tl_t_spatial_20OUT.shp")

#PROPORTION OF THREATENED SPECIES ##############################################

propT <- rep(NA, length(fw_list_with_status))

for(i in 1:length(fw_list_with_status)){
  
  propT_fw2 <- fw_list_with_status[[i]]
  prop <- propT_fw2$agreg_ts
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
#save(proportion_spatial, file = "proportion_spatial_20OUT.Rdata")

writeVector(proportion_spatial, 
            filename = "shapes_20OUT23\\proportion_spatial_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#proportion_spatial <- terra::vect("proportion_spatial_20OUT.shp")

##################################################################################################################
#                                                #STANDARDIZE MAPS
##################################################################################################################
#FMestre
#09-03-2023

#Load packages
library(vegan)
library(terra)

###########
#   IVI
###########

#NT

ivi_nt_std_vector <- as.vector(vegan::decostand(ivi_nt_spatial$ivi, method = "standardize", na.rm = TRUE))
ivi_nt_spatial_STD <- ivi_nt_spatial
ivi_nt_spatial_STD$ivi_STD <- ivi_nt_std_vector
#
writeVector(ivi_nt_spatial_STD, 
            filename = "shapes_20OUT23\\ivi_nt_spatial_STD_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#T

ivi_t_std_vector <- as.vector(vegan::decostand(ivi_t_spatial$ivi, method = "standardize", na.rm = TRUE))
ivi_t_spatial_STD <- ivi_t_spatial
ivi_t_spatial_STD$ivi_STD <- ivi_t_std_vector
#
writeVector(ivi_t_spatial_STD, 
            filename = "shapes_20OUT23\\ivi_t_spatial_STD_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

###########
#CENTRALITY
###########

#NT

centrality_nt_std_vector <- as.vector(vegan::decostand(centrality_nt_spatial$centrality, method = "standardize", na.rm = TRUE))
centrality_nt_std_vector_STD <- centrality_nt_spatial
centrality_nt_std_vector_STD$centrality_STD <- centrality_nt_std_vector
#
writeVector(centrality_nt_std_vector_STD, 
            filename = "shapes_20OUT23\\centrality_nt_std_vector_STD_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#T

centrality_t_std_vector <- as.vector(vegan::decostand(centrality_t_spatial$centrality, method = "standardize", na.rm = TRUE))
centrality_t_std_vector_STD <- centrality_t_spatial
centrality_t_std_vector_STD$centrality_STD <- centrality_t_std_vector
#
writeVector(centrality_t_std_vector_STD, 
            filename = "shapes_20OUT23\\centrality_t_std_vector_STD_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


###########
#IN-DEGREE
###########

#NT

indegree_nt_std_vector <- as.vector(vegan::decostand(indegree_nt_spatial$indegree, method = "standardize", na.rm = TRUE))
indegree_nt_std_vector_STD <- indegree_nt_spatial
indegree_nt_std_vector_STD$indegree_STD <- indegree_nt_std_vector
#
writeVector(indegree_nt_std_vector_STD, 
            filename = "shapes_20OUT23\\indegree_nt_std_vector_STD_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#T

indegree_t_std_vector <- as.vector(vegan::decostand(indegree_t_spatial$indegree, method = "standardize", na.rm = TRUE))
indegree_t_std_vector_STD <- indegree_t_spatial
indegree_t_std_vector_STD$indegree_STD <- indegree_t_std_vector
#
writeVector(indegree_t_std_vector_STD, 
            filename = "shapes_20OUT23\\indegree_t_std_vector_STD_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


###########
#OUT-DEGREE
###########

#NT

outdegree_nt_std_vector <- as.vector(vegan::decostand(outdegree_nt_spatial$outdegree, method = "standardize", na.rm = TRUE))
outdegree_nt_std_vector_STD <- outdegree_nt_spatial
outdegree_nt_std_vector_STD$outdegree_STD <- outdegree_nt_std_vector
#
writeVector(outdegree_nt_std_vector_STD, 
            filename = "shapes_20OUT23\\outdegree_nt_std_vector_STD_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#T

outdegree_t_std_vector <- as.vector(vegan::decostand(outdegree_t_spatial$outdegree, method = "standardize", na.rm = TRUE))
outdegree_t_std_vector_STD <- outdegree_t_spatial
outdegree_t_std_vector_STD$outdegree_STD <- outdegree_t_std_vector
#
writeVector(outdegree_t_std_vector_STD, 
            filename = "shapes_20OUT23\\outdegree_t_std_vector_STD_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


###########
#CLOSENESS
###########

#NT

closeness_nt_std_vector <- as.vector(vegan::decostand(closeness_nt_spatial$closeness, method = "standardize", na.rm = TRUE))
closeness_nt_std_vector_STD <- closeness_nt_spatial
closeness_nt_std_vector_STD$closeness_STD <- closeness_nt_std_vector
#
writeVector(closeness_nt_std_vector_STD, 
            filename = "shapes_20OUT23\\closeness_nt_std_vector_STD_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#T

closeness_t_std_vector <- as.vector(vegan::decostand(closeness_t_spatial$closeness, method = "standardize", na.rm = TRUE))
closeness_t_std_vector_STD <- closeness_t_spatial
closeness_t_std_vector_STD$closeness_STD <- closeness_t_std_vector
#
writeVector(closeness_t_std_vector_STD, 
            filename = "shapes_20OUT23\\closeness_t_std_vector_STD_20OUT.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


################################################################################
#                            COMPARE MAPS - WILCOXON
################################################################################

#FMestre
#15-03-2023

#IVI
ivi_t_spatial$ivi
ivi_nt_spatial$ivi

cor.test(ivi_t_spatial$ivi, ivi_nt_spatial$ivi, method = "spearman", use = "complete.obs")
#
ivi_test <- data.frame(ivi_t_spatial$ivi, ivi_nt_spatial$ivi)
ivi_test <- ivi_test[complete.cases(ivi_test),]
wilcx_ivi <- wilcox.test(ivi_test[,1], ivi_test[,2], paired = TRUE)

#CLOSENESS
closeness_t_spatial$closeness
closeness_nt_spatial$closeness

cor.test(closeness_t_spatial$closeness, closeness_nt_spatial$closeness, method = "spearman", use = "complete.obs")
#
close_test <- data.frame(closeness_t_spatial$closeness, closeness_nt_spatial$closeness)
close_test <- close_test[complete.cases(close_test),]
wilcx_close <- wilcox.test(close_test[,1], close_test[,2], paired = TRUE)

#CENTRALITY
centrality_t_spatial$centrality
centrality_nt_spatial$centrality

cor.test(centrality_t_spatial$centrality, centrality_nt_spatial$centrality, method = "spearman", use = "complete.obs")
#
ctrl_test <- data.frame(centrality_t_spatial$centrality, centrality_nt_spatial$centrality)
ctrl_test <- ctrl_test[complete.cases(ctrl_test),]
ctrl_ivi <- wilcox.test(ctrl_test[,1], ctrl_test[,2], paired = TRUE)


################################################################################
#                             Convert to rasters
################################################################################

template_raster <- terra::rast("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\old_results\\outdeg_diff.tif")

#IVI
nt_ivi <- terra::rasterize(x = ivi_nt_spatial, 
                 y = template_raster, 
                 field = "ivi")

#terra::writeRaster(nt_ivi, "rasters_21OUT\\nt_ivi_20OUT.tif")

t_ivi <- terra::rasterize(x = ivi_t_spatial, 
                           y = template_raster, 
                           field = "ivi")

#terra::writeRaster(t_ivi, "rasters_21OUT\\t_ivi_20OUT.tif")


#Centrality
t_centrality <- terra::rasterize(x = centrality_t_spatial, 
                          y = template_raster, 
                          field = "centrality")

#terra::writeRaster(t_centrality, "rasters_21OUT\\t_centrality_20OUT.tif")

nt_centrality <- terra::rasterize(x = centrality_nt_spatial, 
                                 y = template_raster, 
                                 field = "centrality")

#terra::writeRaster(nt_centrality, "rasters_21OUT\\nt_centrality_20OUT.tif")

#Indegree
t_indegree <- terra::rasterize(x = indegree_t_spatial, 
                                 y = template_raster, 
                                 field = "indegree")

#terra::writeRaster(t_indegree, "rasters_21OUT\\t_indegree_20OUT.tif")

nt_indegree <- terra::rasterize(x = indegree_nt_spatial, 
                                  y = template_raster, 
                                  field = "indegree")

#terra::writeRaster(nt_indegree, "rasters_21OUT\\nt_indegree_20OUT.tif")

#Outdegree
nt_outdegree <- terra::rasterize(x = outdegree_nt_spatial, 
                                y = template_raster, 
                                field = "outdegree")

#terra::writeRaster(nt_outdegree, "rasters_21OUT\\nt_outdegree_20OUT.tif")

t_outdegree <- terra::rasterize(x = outdegree_t_spatial, 
                                 y = template_raster, 
                                 field = "outdegree")

#terra::writeRaster(t_outdegree, "rasters_21OUT\\t_outdegree_20OUT.tif")


#Closeness
t_closeness <- terra::rasterize(x = closeness_t_spatial, 
                                y = template_raster, 
                                field = "closeness")

#terra::writeRaster(t_closeness, "rasters_21OUT\\t_closeness_20OUT.tif")

nt_closeness <- terra::rasterize(x = closeness_nt_spatial, 
                                y = template_raster, 
                                field = "closeness")

#terra::writeRaster(nt_closeness, "rasters_21OUT\\nt_closeness_20OUT.tif")

#TL
t_tl <- terra::rasterize(x = tl_t_spatial, 
                                 y = template_raster, 
                                 field = "trophic_level")

#terra::writeRaster(t_tl, "rasters_21OUT\\t_tl_20OUT.tif")

nt_tl <- terra::rasterize(x = tl_nt_spatial, 
                          y = template_raster, 
                          field = "trophic_level")

#terra::writeRaster(nt_tl, "rasters_21OUT\\nt_tl_20OUT.tif")


#Proportion
proportion_r <- terra::rasterize(x = proportion_spatial, 
                          y = template_raster, 
                          field = "proportion")

#terra::writeRaster(proportion_r, "rasters_21OUT\\proportion_r_20OUT.tif")

################################################################################
#                          Plot maps and frequency           
################################################################################

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
rasterVis::levelplot(proportion_r, par.settings = rasterTheme(viridis_pal()(255)), main = "Proportion of threatened species")
rasterVis::levelplot(xxx, par.settings = rasterTheme(viridis_pal()(255)), main = "Number of nodes in each network")

################################################################################
#                             Histogram plots
################################################################################

################################### indegree ###################################

nt_indegree_df <- as.data.frame(indegree_nt_spatial)
t_indegree_df <- as.data.frame(indegree_t_spatial)
#head(nt_indegree_df)
#head(t_indegree_df)
nt_indegree_df <- nt_indegree_df[,c(1,4)]
t_indegree_df <- t_indegree_df[,c(1,4)]
#
names(nt_indegree_df)[2] <- "NT_indegree"
names(t_indegree_df)[2] <- "T_indegree"
#
indegree_compare <- merge(nt_indegree_df, t_indegree_df)
indegree_compare <- indegree_compare[complete.cases(indegree_compare),]
#View(indegree_compare)
indegree_ttest <- t.test(indegree_compare[,2], indegree_compare[,3], paired = TRUE)
#print(indegree_ttest)
indegree_cohens_d <- effsize::cohen.d(indegree_compare[,2], indegree_compare[,3])
#print(indegree_cohens_d)
#plot(indegree_compare[,2], indegree_compare[,3])
cor(indegree_compare[,2], indegree_compare[,3], method = "pearson")

####### Plot

#frequency plots
min(indegree_compare[,2])
max(indegree_compare[,2])
min(indegree_compare[,3])
max(indegree_compare[,3])

#names(indegree_compare)
par(mfrow=c(2, 1))
h_in_nt <- hist(indegree_compare[,2])
h_in_t <- hist(indegree_compare[,3])
# Convert the counts to percentages
h_in_nt$counts <- h_in_nt$counts / sum(h_in_nt$counts) * 100
h_in_t$counts <- h_in_t$counts / sum(h_in_t$counts) * 100
# Create the line plot
plot(h_in_nt$mids, h_in_nt$counts, type = "n", xlab = "In-degree", ylab = "Frequency (%)", main = "In-degree", ylim = c(0,45))
lines(h_in_nt$mids, h_in_nt$counts, lwd = 3, col = "darkgreen")
lines(h_in_t$mids, h_in_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)

################################## outdegree ###################################

nt_outdegree_df <- as.data.frame(outdegree_nt_spatial)
t_outegree_df <- as.data.frame(outdegree_t_spatial)
#head(nt_outdegree_df)
#head(t_outdegree_df)
nt_outdegree_df <- nt_outdegree_df[,c(1,4)]
t_outdegree_df <- t_outdegree_df[,c(1,4)]
#
names(nt_outdegree_df)[2] <- "NT_outdegree"
names(t_outdegree_df)[2] <- "T_outdegree"
#
outdegree_compare <- merge(nt_outdegree_df, t_outdegree_df)
outdegree_compare <- outdegree_compare[complete.cases(outdegree_compare),]
#View(outdegree_compare)
outdegree_ttest <- t.test(outdegree_compare[,2], outdegree_compare[,3], paired = TRUE)
#print(outdegree_ttest)
outdegree_cohens_d <- effsize::cohen.d(outdegree_compare[,2], outdegree_compare[,3])
#print(outdegree_cohens_d)
#plot(outdegree_compare[,2], outdegree_compare[,3])
cor(outdegree_compare[,2], outdegree_compare[,3], method = "pearson")

####### Plot

#frequency plots
min(outdegree_compare[,2])
max(outdegree_compare[,2])
min(outdegree_compare[,3])
max(outdegree_compare[,3])

#names(indegree_compare)
par(mfrow=c(2, 1))
h_out_nt <- hist(outdegree_compare[,2])
h_out_t <- hist(outdegree_compare[,3])
# Convert the counts to percentages
h_out_nt$counts <- h_out_nt$counts / sum(h_out_nt$counts) * 100
h_out_t$counts <- h_out_t$counts / sum(h_out_t$counts) * 100
# Create the line plot
plot(h_out_nt$mids, h_out_nt$counts, type = "n", xlab = "In-degree", ylab = "Frequency (%)", main = "Out-degree", ylim = c(0,25))
lines(h_out_nt$mids, h_out_nt$counts, lwd = 3, col = "darkgreen")
lines(h_out_t$mids, h_out_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)







################################################################################
#                          Compute differences           
################################################################################

indeg_diff <- t_indegree - nt_indegree
terra::writeRaster(indeg_diff, filename = "rasters_21OUT\\indeg_diff.tif")
#plot(indeg_diff)
#
outdeg_diff <- t_outdegree - nt_outdegree
terra::writeRaster(outdeg_diff, filename = "rasters_21OUT\\outdeg_diff.tif")
#plot(outdeg_diff)
#
trophic_level_diff <- t_tl - nt_tl
terra::writeRaster(trophic_level_diff, filename = "rasters_21OUT\\trophic_level_diff.tif")
#plot(trophic_level_diff)
#
closeness_diff <- t_closeness - nt_closeness
terra::writeRaster(closeness_diff, filename = "rasters_21OUT\\closeness_diff.tif")
#plot(closeness_diff)
#
centrality_diff <- t_centrality - nt_centrality
terra::writeRaster(centrality_diff, filename = "rasters_21OUT\\centrality_diff.tif")
#plot(centrality_diff)
#
ivi_diff <- t_ivi - nt_ivi
terra::writeRaster(ivi_diff, filename = "rasters_21OUT\\ivi_diff.tif")
#plot(ivi_diff)

