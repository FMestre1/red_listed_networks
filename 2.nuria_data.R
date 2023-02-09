#FMestre
#07-02-2023

library(terra)

grid <- terra::vect("C:/Users/FMest/Documents/0. Artigos/IUCN_networks/shapefiles/Europa_10km/europe_10km.shp")
crs(grid)

europe <- terra::vect("C:/Users/FMest/Documents/0. Artigos/IUCN_networks/shapefiles/Europe/Europe.shp")
crs(europe)

#Convert CRS
grid_wgs84 <- terra::project(grid, europe)

plot(europe)
plot(grid_wgs84, add=TRUE)


####Load data
#C:\Users\FMest\Documents\0. Artigos\IUCN_networks\data\data_nuria

path1 <- "C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria\\"
files_folder <- list.files(path1)
#files_folder[[i]]

fw_list <- list()
#
for(i in 1:length(files_folder)){
  fw_list[[i]] <- read.table(paste0(path1, files_folder[i]), sep=",", header = 1)
  names(fw_list)[i] <- stringr::str_split(files_folder[i], ".csv")[[1]][1]
  message(i)
}
#
save(fw_list, file = "fw_list.RData")

#Get the species names
species_names <- c()

for(i in 1:length(fw_list)){
  
  df1 <- fw_list[[i]]
  sp1 <- df1$SP_NAME
  species_names <- c(species_names, sp1)
  
  message(i)
}

species_names2 <- unique(species_names)
length(species_names)
length(species_names2)

save(species_names2, file = "species_names2.RData")
