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
plot(europeRaster_poly)
head(europeRaster_poly)
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

#Add ICUN status
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
  if(length(fw_nt_ivi)!=0) ivi_nt[i] <- mean(fw_nt_ivi)
  
  message(i)
  
}

ivi_nt <- data.frame(names(fw_list_with_status_aggreg), ivi_nt)
names(ivi_nt) <- c("grid", "ivi")

ivi_nt_spatial <- merge(europeRaster_poly, ivi_nt, by.x = "PageName", by.y = "grid")
plot(ivi_nt_spatial)
plot(europeRaster_poly)

writeVector(ivi_nt_spatial, 
            filename = "ivi_nt_spatial.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=FALSE, 
            options="ENCODING=UTF-8"
)
