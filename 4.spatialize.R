#FMestre
#09-02-2023

library(terra)

grid <- terra::vect("C:/Users/FMest/Documents/0. Artigos/IUCN_networks/shapefiles/Europa_10km/europe_10km.shp")
crs(grid)

europe <- terra::vect("C:/Users/FMest/Documents/0. Artigos/IUCN_networks/shapefiles/Europe/Europe.shp")
crs(europe)

#Convert CRS
grid_wgs84 <- terra::project(grid, europe)
#
#plot(europe)
#plot(grid_wgs84, add = TRUE)

europe_coastline_borders <- aggregate(europe, dissolve = TRUE)
#plot(europe_coastline_borders)
#crs(europe_coastline_borders)

#Save shapefile
writeVector(europe_coastline_borders, 
            filename = "europe_coastline_borders.shp", 
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=FALSE, 
            options="ENCODING=UTF-8"
            )

#grid_europe <- intersect(grid_wgs84, europe_coastline_borders) #takes toooo long! I'll do this in ArcGis.

grid_europe <- terra::vect("C:/Users/FMest/Documents/0. Artigos/IUCN_networks/shapefiles/grid_10_EUROPE.shp")
crs(grid_europe)

grid_europe_wgs84 <- terra::project(grid_europe, europe)
crs(grid_europe_wgs84)
#plot(europe_coastline_borders)
#plot(grid_europe_wgs84, add = TRUE)

#Codes from the grid downloaded
#codes_from_grids <- grid_europe_wgs84$CellCode

#Codes from Núria´s dataset
#codes_from_nuria <- names(fw_list)

#The codes do not match. Asked Núria for the shapefile.

europeRaster <- terra::rast(x="C:/Users/FMest/Documents/github/red_listed_networks/mask10k-20230214T144658Z-001/mask10k/reference_grid_10km.img")
cells_info <- foreign::read.dbf(file = "C:/Users/FMest/Documents/github/red_listed_networks/mask10k-20230214T144658Z-001/mask10k/reference_grid_10km.img.vat.dbf")
#head(cells_info)
nrow(cells_info)

#To vector
europeRaster_poly <- terra::as.polygons(europeRaster, values = TRUE, extent=FALSE)
europeRaster_poly <- merge(europeRaster_poly, cells_info)
#plot(europeRaster_poly)
#head(europeRaster_poly)
europeRaster_poly_wgs84 <- terra::project(europeRaster_poly, europe)
#
europeRaster_poly_wgs84_coords <- crds(europeRaster_poly_wgs84, df=TRUE)
europeRaster_poly_wgs84_coords <- data.frame(europeRaster_poly_wgs84, europeRaster_poly_wgs84_coords)
#head(europeRaster_poly_wgs84_coords)

#Write vector
#writeVector(europeRaster_poly, filename ="europeRaster_poly.shp", overwrite=TRUE, filetype = "ESRI Shapefile")

################################################################################

red_listed_3
fw_list

#as.numeric(lapply(fw_list, nrow))

fw_list_with_status <- fw_list

#Add IUCN status
for(i in 1:length(fw_list_with_status)){
  
  fw3 <- fw_list_with_status[[i]]
  fw3$SP_NAME <- stringr::str_replace(fw3$SP_NAME, "_", " ")
  sp_fw3 <- stringr::str_replace(fw3$SP_NAME, "_", " ")
  
  if(any(red_listed_3$full_name %in% sp_fw3))
    {
    sp_fw3_redList <- red_listed_3[red_listed_3$full_name %in% sp_fw3,]
    fw4 <- merge(fw3, sp_fw3_redList, by.x = "SP_NAME", by.y = "full_name")
    fw_list_with_status[[i]] <- fw4
  }


message(i)

}

fw_list_with_status_aggreg <- fw_list_with_status

#Add threatened/non-threatened
for(i in 1:length(fw_list_with_status)){
  
  fw5 <- fw_list_with_status[[i]]
  
  if(nrow(fw5)!=0){
  
  categories_fw5 <- fw5$europeanRegionalRedListCategory
  #
  categories_fw5 <- stringr::str_replace(categories_fw5, "VU", "threatened")
  categories_fw5 <- stringr::str_replace(categories_fw5, "EN", "threatened")
  categories_fw5 <- stringr::str_replace(categories_fw5, "CR", "threatened")
  #
  categories_fw5 <- stringr::str_replace(categories_fw5, "LC", "non-threatened")
  categories_fw5 <- stringr::str_replace(categories_fw5, "NT", "non-threatened")
  #
  categories_fw5 <- stringr::str_replace(categories_fw5, "DD", "others")
  categories_fw5 <- stringr::str_replace(categories_fw5, "NE", "others")
  categories_fw5 <- stringr::str_replace(categories_fw5, "RE", "others")
  #
  
  fw5 <- data.frame(fw5, categories_fw5)
  names(fw5)[12] <- "aggreg_IUCN"
  
  fw_list_with_status_aggreg[[i]] <- fw5
  }
  message(i)
  
}


#IVI - NOT THREATENED ##########################################################

ivi_nt <- rep(NA, length(fw_list_with_status_aggreg))

for(i in 1:length(fw_list_with_status_aggreg)){
  
  ivi_fw <- fw_list_with_status_aggreg[[i]]
  fw_nt_ivi <- ivi_fw[ivi_fw$aggreg_IUCN == "non-threatened",]$ivi
  if(length(fw_nt_ivi)!=0) ivi_nt[i] <- mean(fw_nt_ivi, na.rm=TRUE)
  
  message(i)
  
}

ivi_nt <- data.frame(names(fw_list_with_status_aggreg), ivi_nt)
names(ivi_nt) <- c("grid", "ivi")
#head(ivi_nt)
#hist(ivi_nt$ivi)

ivi_nt_spatial <- merge(europeRaster_poly, ivi_nt, by.x = "PageName", by.y = "grid")
#plot(ivi_nt_spatial)
#plot(europeRaster_poly)

writeVector(ivi_nt_spatial, 
            filename = "ivi_nt_spatial_second_version.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#IVI - THREATENED ##############################################################

ivi_t <- rep(NA, length(fw_list_with_status_aggreg))

for(i in 1:length(fw_list_with_status_aggreg)){
  
  ivi_fw2 <- fw_list_with_status_aggreg[[i]]
  fw_t_ivi2 <- ivi_fw2[ivi_fw2$aggreg_IUCN == "threatened",]$ivi
  if(length(fw_t_ivi2)!=0) ivi_t[i] <- mean(fw_t_ivi2, na.rm=TRUE)
  
  message(i)
  
}

ivi_t <- data.frame(names(fw_list_with_status_aggreg), ivi_t)
names(ivi_t) <- c("grid", "ivi")
#head(ivi_t)
#hist(ivi_t$ivi)

ivi_t_spatial <- merge(europeRaster_poly, ivi_t, by.x = "PageName", by.y = "grid")
#plot(ivi_t_spatial)

writeVector(ivi_t_spatial, 
            filename = "ivi_t_spatial_second_version.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#CENTRALITY - NON-THREATENED #######################################################

centrality_nt <- rep(NA, length(fw_list_with_status_aggreg))

for(i in 1:length(fw_list_with_status_aggreg)){
  
  cent_nt_fw2 <- fw_list_with_status_aggreg[[i]]
  fw_nt_cent2 <- cent_nt_fw2[cent_nt_fw2$aggreg_IUCN == "non-threatened",]$centrality
  if(length(fw_nt_cent2)!=0) centrality_nt[i] <- mean(fw_nt_cent2, na.rm=TRUE)
  
  message(i)
  
}

centrality_nt <- data.frame(names(fw_list_with_status_aggreg), centrality_nt)
names(centrality_nt) <- c("grid", "centrality")
#head(centrality_nt)
#hist(centrality_nt$centrality)

centrality_nt_spatial <- merge(europeRaster_poly, centrality_nt, by.x = "PageName", by.y = "grid")
#plot(centrality_nt_spatial)

writeVector(centrality_nt_spatial, 
            filename = "centrality_nt_spatial.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#CENTRALITY - THREATENED #######################################################

centrality_t <- rep(NA, length(fw_list_with_status_aggreg))

for(i in 1:length(fw_list_with_status_aggreg)){
  
  cent_fw2 <- fw_list_with_status_aggreg[[i]]
  fw_t_cent2 <- cent_fw2[cent_fw2$aggreg_IUCN == "threatened",]$centrality
  if(length(fw_t_cent2)!=0) centrality_t[i] <- mean(fw_t_cent2, na.rm=TRUE)
  
  message(i)
  
}

centrality_t <- data.frame(names(fw_list_with_status_aggreg), centrality_t)
names(centrality_t) <- c("grid", "centrality")
#head(centrality_t)
#hist(centrality_t$centrality)

centrality_t_spatial <- merge(europeRaster_poly, centrality_t, by.x = "PageName", by.y = "grid")
#plot(centrality_t_spatial)

writeVector(centrality_t_spatial, 
            filename = "centrality_t_spatial.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


#IN-DEGREE - THREATENED ########################################################

indegree_t <- rep(NA, length(fw_list_with_status_aggreg))

for(i in 1:length(fw_list_with_status_aggreg)){
  
  indeg_fw2 <- fw_list_with_status_aggreg[[i]]
  fw_t_indeg2 <- indeg_fw2[indeg_fw2$aggreg_IUCN == "threatened",]$indegree
  if(length(fw_t_indeg2)!=0) indegree_t[i] <- mean(fw_t_indeg2, na.rm=TRUE)
  
  message(i)
  
}

indegree_t <- data.frame(names(fw_list_with_status_aggreg), indegree_t)
names(indegree_t) <- c("grid", "indegree")
#head(indegree_t)
#hist(indegree_t$indegree)

indegree_t_spatial <- merge(europeRaster_poly, indegree_t, by.x = "PageName", by.y = "grid")
#plot(indegree_t_spatial)

writeVector(indegree_t_spatial, 
            filename = "indegree_t_spatial.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


#IN-DEGREE - NON-THREATENED ####################################################

indegree_nt <- rep(NA, length(fw_list_with_status_aggreg))

for(i in 1:length(fw_list_with_status_aggreg)){
  
  indeg_nt_fw2 <- fw_list_with_status_aggreg[[i]]
  fw_nt_indeg2 <- indeg_nt_fw2[indeg_nt_fw2$aggreg_IUCN == "non-threatened",]$indegree
  if(length(fw_nt_indeg2)!=0) indegree_nt[i] <- mean(fw_nt_indeg2, na.rm=TRUE)
  
  message(i)
  
}

indegree_nt <- data.frame(names(fw_list_with_status_aggreg), indegree_nt)
names(indegree_nt) <- c("grid", "indegree")
#head(indegree_nt)
#hist(indegree_nt$indegree)

indegree_nt_spatial <- merge(europeRaster_poly, indegree_nt, by.x = "PageName", by.y = "grid")
#plot(indegree_nt_spatial)

writeVector(indegree_nt_spatial, 
            filename = "indegree_nt_spatial.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


#CLOSENESS - THREATENED ########################################################

closeness_t <- rep(NA, length(fw_list_with_status_aggreg))

for(i in 1:length(fw_list_with_status_aggreg)){
  
  close_t_fw2 <- fw_list_with_status_aggreg[[i]]
  fw_t_close2 <- close_t_fw2[close_t_fw2$aggreg_IUCN == "threatened",]$closeness
  if(length(fw_t_close2)!=0) closeness_t[i] <- mean(fw_t_close2, na.rm=TRUE)
  
  message(i)
  
}

closeness_t <- data.frame(names(fw_list_with_status_aggreg), closeness_t)
names(closeness_t) <- c("grid", "closeness")
#head(closeness_t)
#hist(closeness_t$closeness)

closeness_t_spatial <- merge(europeRaster_poly, closeness_t, by.x = "PageName", by.y = "grid")
#plot(closeness_t_spatial)

writeVector(closeness_t_spatial, 
            filename = "closeness_t_spatial.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#CLOSENESS - NON-THREATENED ####################################################

closeness_nt <- rep(NA, length(fw_list_with_status_aggreg))

for(i in 1:length(fw_list_with_status_aggreg)){
  
  close_nt_fw2 <- fw_list_with_status_aggreg[[i]]
  fw_nt_close2 <- close_nt_fw2[close_nt_fw2$aggreg_IUCN == "non-threatened",]$closeness
  if(length(fw_nt_close2)!=0) closeness_nt[i] <- mean(fw_nt_close2, na.rm=TRUE)
  
  message(i)
  
}

closeness_nt <- data.frame(names(fw_list_with_status_aggreg), closeness_nt)
names(closeness_nt) <- c("grid", "closeness")
#head(closeness_nt)
#hist(closeness_nt$closeness)

closeness_nt_spatial <- merge(europeRaster_poly, closeness_nt, by.x = "PageName", by.y = "grid")
#plot(closeness_t_spatial)

writeVector(closeness_nt_spatial, 
            filename = "closeness_nt_spatial.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


#TROPHIC LEVEL - NON-THREATENED ################################################

TL_nt <- rep(NA, length(fw_list_with_status_aggreg))

for(i in 1:length(fw_list_with_status_aggreg)){
  
  tl_nt_fw2 <- fw_list_with_status_aggreg[[i]]
  fw_nt_TL2 <- tl_nt_fw2[tl_nt_fw2$aggreg_IUCN == "non-threatened",]$TL
  if(length(fw_nt_TL2)!=0) TL_nt[i] <- mean(fw_nt_TL2, na.rm=TRUE)
  
  message(i)
  
}

TL_nt <- data.frame(names(fw_list_with_status_aggreg), TL_nt)
names(TL_nt) <- c("grid", "trophic_level")
#head(TL_nt)
#hist(TL_nt$trophic_level)

tl_nt_spatial <- merge(europeRaster_poly, TL_nt, by.x = "PageName", by.y = "grid")
#plot(closeness_t_spatial)

writeVector(tl_nt_spatial, 
            filename = "tl_nt_spatial.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


#TROPHIC LEVEL - THREATENED ####################################################

TL_t <- rep(NA, length(fw_list_with_status_aggreg))

for(i in 1:length(fw_list_with_status_aggreg)){
  
  tl_t_fw2 <- fw_list_with_status_aggreg[[i]]
  fw_t_TL2 <- tl_t_fw2[tl_t_fw2$aggreg_IUCN == "threatened",]$TL
  if(length(fw_t_TL2)!=0) TL_t[i] <- mean(fw_t_TL2, na.rm=TRUE)
  
  message(i)
  
}

TL_t <- data.frame(names(fw_list_with_status_aggreg), TL_t)
names(TL_t) <- c("grid", "trophic_level")
#head(TL_t)
#hist(TL_t$trophic_level)

tl_t_spatial <- merge(europeRaster_poly, TL_t, by.x = "PageName", by.y = "grid")
#plot(closeness_t_spatial)

writeVector(tl_t_spatial, 
            filename = "tl_t_spatial.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#PROPORTION OF THREATENED SPECIES ##############################################

propT <- rep(NA, length(fw_list_with_status_aggreg))

for(i in 1:length(fw_list_with_status_aggreg)){
  
  propT_fw2 <- fw_list_with_status_aggreg[[i]]
  prop <- propT_fw2$aggreg_IUCN
  n_t <- length(prop[prop=="threatened"])
  n_total <- length(prop)
  propT[i] <- n_t/n_total
  
  #if(length(fw_t_TL2)!=0) TL_t[i] <- mean(fw_t_TL2)
  
  message(i)
  
}

propT <- data.frame(names(fw_list_with_status_aggreg), propT)
names(propT) <- c("grid", "proportion")
#head(propT)
#hist(propT$proportion)

proportion_spatial <- merge(europeRaster_poly, propT, by.x = "PageName", by.y = "grid")
#plot(proportion_spatial)

writeVector(proportion_spatial, 
            filename = "proportion_spatial.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


#df123 <- plyr::ldply(fw_list_with_status_aggreg, data.frame)
#table(df123$aggreg_IUCN)
