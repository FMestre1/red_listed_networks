#FMestre
#03-05-2023

#Paper: https://doi.org/10.1111/ele.12002

#library(bipartite)
#library(combinat)
#library(factoextra)
#library(terra)

########Incomplete code... bellow!

#df_S_combined <- as.data.frame(matrix(ncol = length(fw_list_with_status_aggreg_BS), 
#                                      nrow = length(fw_list_with_status_aggreg_BS))
#)


#rownames(df_S_combined) <- names(fw_list_with_status_aggreg_BS)
#colnames(df_S_combined) <- names(fw_list_with_status_aggreg_BS)

#df_OS_combined <- df_S_combined 
#df_WN_combined <- df_S_combined 
#df_ST_combined <- df_S_combined 

#grid_names <- names(iberian_fw_new_comm_collection_new_version_28SET2022)

#matrix_comparisons_combined <- data.frame(t(combn(grid_names, 2, simplify = TRUE)))
#nrow(matrix_comparisons_combined)

##

#Save these i
#save_i_combined <- seq(from = 1, to = nrow(matrix_comparisons_combined), by = round(nrow(matrix_comparisons_combined)/1000, 0))
#save_i_combined <- c(save_i_combined, nrow(matrix_comparisons_combined))

#for(i in 1:nrow(matrix_comparisons_combined)){
  
#  row1 <-  matrix_comparisons_combined[i,]
  
#  nt1 <- as.character(row1[1])
#  nt2 <- as.character(row1[2])
  
#  adj1 <- as.matrix(combined_dataset_globi_maiorano_tg_ADJACENCY_MATRIX[[nt1]])
#  adj2 <- as.matrix(combined_dataset_globi_maiorano_tg_ADJACENCY_MATRIX[[nt2]])
  
#  try(out1 <- betalinkr(webs2array(list(adj1=adj1, adj2=adj2)), 
#                        partitioning="commondenom", 
#                        partition.st=FALSE
#  ))
  
#  df_S_combined[nt1,nt2] <- as.numeric(out1["S"])
#  df_OS_combined[nt1,nt2] <- as.numeric(out1["OS"])
#  df_WN_combined[nt1,nt2] <- as.numeric(out1["WN"])
#  df_ST_combined[nt1,nt2] <- as.numeric(out1["ST"])
  
#  message(round((i*100)/nrow(matrix_comparisons_combined), 2))
  
#  if(i !=0) rm(out1)
  
#  if(any(i == save_i_combined)) {
#    save(df_S_combined, file = paste0("matrix_out_combined/df_S_combined_16dez23_", i, ".RData"))
#    save(df_OS_combined, file = paste0("matrix_out_combined/df_OS_combined_16dez23_",i, ".RData"))
#    save(df_WN_combined, file = paste0("matrix_out_combined/df_WN_combined_16dez23_", i, ".RData"))
#    save(df_ST_combined, file = paste0("matrix_out_combined/df_ST_combined_16dez23_", i, ".RData"))
    
#  }
  
#}#Ran this at the cluster 

########Incomplete code... up!

################################################################################
################################################################################
################################################################################
#FMestre
#05-04-2023

library(terra)
library(factoextra)

europe <- terra::vect("C:/Users/FMest/Documents/0. Artigos/IUCN_networks/shapefiles/Europe/Europe.shp")
crs(europe)

#Add modularity to the grid ####################################################

  #where's the grid? ##########
europeRaster <- terra::rast(x="C:/Users/FMest/Documents/github/red_listed_networks/mask10k-20230214T144658Z-001/mask10k/reference_grid_10km.img")
cells_info <- foreign::read.dbf(file = "C:/Users/FMest/Documents/github/red_listed_networks/mask10k-20230214T144658Z-001/mask10k/reference_grid_10km.img.vat.dbf")
#head(cells_info)
#nrow(cells_info)
#To vector
europeRaster_poly_0 <- terra::as.polygons(europeRaster, values = TRUE, extent=FALSE)
europeRaster_poly <- merge(europeRaster_poly_0, cells_info)
plot(europeRaster_poly)
class(europeRaster_poly)
#
plot(europeRaster_poly_0)
class(europeRaster_poly_0)

  #where's the modularity (and other metrics)? ##########

path2 <- "C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_2\\"

files_folder2 <- list.files(path2)
#length(files_folder2)

fw_list2 <- data.frame(rep(NA, length(files_folder2)), 
                       rep(NA, length(files_folder2)),
                       rep(NA, length(files_folder2)),
                       rep(NA, length(files_folder2)), 
                       rep(NA, length(files_folder2)),
                       rep(NA, length(files_folder2))
                       )
colnames(fw_list2) <- c("network", "modularity", "connectance", "basal", "top", "intermediate")
head(fw_list2)

#
for(i in 1:length(files_folder2)){
  fw_list2[i,2] <- read.table(paste0(path2, files_folder2[i]), sep=",", header = 1)$modularity
  fw_list2[i,1] <- stringr::str_split(files_folder2[i], ".csv")[[1]][1]
  fw_list2[i,3] <- read.table(paste0(path2, files_folder2[i]), sep=",", header = 1)$C
  fw_list2[i,4] <- read.table(paste0(path2, files_folder2[i]), sep=",", header = 1)$basal
  fw_list2[i,5] <- read.table(paste0(path2, files_folder2[i]), sep=",", header = 1)$top
  fw_list2[i,6] <- read.table(paste0(path2, files_folder2[i]), sep=",", header = 1)$intermediate
  
  message(i)
}
#

  #Let's bring them together ... ##########
modularity_spatial <- merge(europeRaster_poly_0, fw_list2, by.x = "PageName", by.y = "network")
modularity_spatial_wgs84 <- terra::project(modularity_spatial, europe)
plot(modularity_spatial_wgs84)

#Write vector
#writeVector(modularity_spatial_wgs84, filename ="modularity_spatial_wgs84.shp", overwrite=TRUE, filetype = "ESRI Shapefile")

  #... and plot it... ##########
#plot(modularity_spatial_wgs84)

#Clustering ####################################################################

modularity_spatial_2 <- as.data.frame(modularity_spatial)
modularity_spatial_2$modularity

# Compute the WSS for different numbers of clusters
#wss <- numeric(10)
#for (k in 1:10) {
#  kmeans_model <- kmeans(modularity_spatial_2$modularity, centers = k)
#  wss[k] <- kmeans_model$tot.withinss
#}

# Plot the elbow curve

modularity_spatial_2_scaled <- scale(modularity_spatial_2[sample(nrow(modularity_spatial_2), round(nrow(modularity_spatial_2)/5)), 4:8])
#?fviz_nbclust

#Next code takes too long...
#fviz_nbclust(modularity_spatial_2_scaled, kmeans, method = "silhouette", k.max = 10, verbose =TRUE)

# Initialize the k-means model with the number of clusters you want
kmeans_model <- kmeans(modularity_spatial_2[,4:8], centers = 2)

# Get the cluster labels and centroids
labels1 <- kmeans_model$cluster
centroids1 <- kmeans_model$centers
#
labels1
centroids1

length(labels1)
nrow(europeRaster_poly_0)

#Create df
clusters_df <- data.frame(modularity_spatial_2$PageName, labels1)
#View(clusters_df)
names(clusters_df) <- c("grid", "clusters")
#Save
write.csv(clusters_df, file = "clusters_df.csv")


