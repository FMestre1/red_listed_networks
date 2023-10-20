################################################################################
#                    FIRST NÚRIA DATASET (node metrics)
################################################################################
#FMestre
#07-02-2023

library(igraph)
library(cheddar)

####Load data
path1 <- "C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria\\"
files_folder <- list.files(path1)

fw_list <- list()
#
for(i in 1:length(files_folder)){
  fw_list[[i]] <- read.table(paste0(path1, files_folder[i]), sep=",", header = 1)
  names(fw_list)[i] <- stringr::str_split(files_folder[i], ".csv")[[1]][1]
  message(i)
}
#
#save(fw_list, file = "fw_list.RData")
#load("fw_list.RData")

#Get the species names
species_names <- c()

for(i in 1:length(fw_list)){
  
  df1 <- fw_list[[i]]
  sp1 <- df1$SP_NAME
  species_names <- c(species_names, sp1)
  
  message(i)
}

species_names2 <- unique(species_names)

#Load & Save
#load("species_names2.RData")
#save(species_names2, file = "species_names2.RData")

################################################################################
#                   SECOND NÚRIA DATASET (network metrics)
################################################################################

dataset2_nuria <- list.files("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_2")
#length(dataset2_nuria)

metrics_dataset_2 <- data.frame()
#head(metrics_dataset_2)

for(i in 1:length(dataset2_nuria)){
  
  m123 <- read.csv(paste0("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_2\\", dataset2_nuria[i]))
  grid_code <- strsplit(dataset2_nuria[i], ".csv")[[1]]
  metrics_dataset_2 <- rbind(metrics_dataset_2, data.frame(grid_code,m123))
  #
  message(paste0("Did ", i, "!"))
  
}
#View(metrics_dataset_2)
#nrow(metrics_dataset_2)

metrics_dataset_3 <- unique(metrics_dataset_2)
#View(metrics_dataset_3)
#nrow(metrics_dataset_3)

#Load & Save
#load("metrics_dataset_FINAL.RData")
#save(metrics_dataset_3, file = "metrics_dataset_FINAL_26SET2023.RData")

################################################################################
#                       THIRD NÚRIA DATASET (networks)
################################################################################

#FMestre
#22-05-23

network_list <- list()
path3 <- "C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_3"

#length(network_list)

dataset3_nuria <- list.files("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_3")
#length(dataset3_nuria)

for(i in 1:length(dataset3_nuria)){
  network_list[[i]] <- read.table(paste0(path3, "\\",dataset3_nuria[i]), sep=",", header = 1)
  names(network_list)[i] <- stringr::str_split(dataset3_nuria[i], ".csv")[[1]][1]
  message(i)
}

#Load & Save
#save(network_list, file = "network_list24SET23.RData")

#Convert to igraph and cheddar

## in the cluster ## START
network_list_igraph <- vector(mode='list', length=length(network_list))
names(network_list_igraph) <- names(network_list)
#
network_list_cheddar <- network_list_igraph

for(i in 1:length(network_list)){
  
  fw_a <- network_list[[i]]
  rownames(fw_a) <- fw_a$X
  fw_a <- fw_a[,-1]
  
  nodes1 <- colnames(fw_a)
  nodes1 <- stringr::str_replace_all(nodes1, pattern = "_", replacement = " ")
  nodes1 <- all_species_status_body_mass_amph_13[all_species_status_body_mass_amph_13$species %in% nodes1,]
  #
  edges1 <- as.data.frame(which(fw_a==1, arr.ind = TRUE))
  rownames(edges1) <- 1:nrow(edges1)
  
  for(j in 1:nrow(edges1)){
  row_a <- edges1[j,]
  fw_a
  edges1[j,2] <- colnames(fw_a)[as.numeric(row_a[2])]
  edges1[j,1] <- rownames(fw_a)[as.numeric(row_a[1])]
  }
  names(edges1) <- c("resource", "consumer")
  edges1$resource <- stringr::str_replace_all(edges1$resource, pattern = "_", replacement = " ")
  edges1$consumer <- stringr::str_replace_all(edges1$consumer, pattern = "_", replacement = " ")
  
  g1 <- graph_from_data_frame(edges1, directed = TRUE, vertices = nodes1)
  
  #Insert igraph
  network_list_igraph[[i]] <- g1
  
  #... no the cheddar version
  names(nodes1)[1] <- "node"
  ch1 <- cheddar::Community(nodes = nodes1, properties = list(title= paste0("Interactions in grid ", names(network_list[i]))), trophic.links = edges1)
  network_list_cheddar[[i]] <- ch1
  
  message(i)
  
}

#I had to run this in chunks...
#...combining partial lists after:

#igraph list

#load("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\new_sept23_dataset\\network_list_igraph_until_38636.RData")
list1_igraph <- get("network_list_igraph")
rm(network_list_igraph)
list1_igraph[[1]]
list1_igraph[[38635]]

#load("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\new_sept23_dataset\\network_list_igraph_until_from_38636_until_54275.RData")
list2_igraph <- get("network_list_igraph")
rm(network_list_igraph)
list2_igraph[[38636]]
list2_igraph[[54275]]

#load("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\new_sept23_dataset\\network_list_igraph_until_from_54276_until_68933.RData")
list3_igraph <- get("network_list_igraph")
rm(network_list_igraph)
list3_igraph[[54276]]
list3_igraph[[68933]]

#load("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\new_sept23_dataset\\network_list_igraph_from_68934_until_70129.RData")
list4_igraph <- get("network_list_igraph")
rm(network_list_igraph)
list4_igraph[[68934]]
list4_igraph[[70129]]

#load("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\new_sept23_dataset\\network_list_igraph_from_70130_until_118292.RData")
list5_igraph <- get("network_list_igraph")
rm(network_list_igraph)
list5_igraph[[70130]]
list5_igraph[[118292]]

#Gather all the graph in one list
igraph_list <- c(
  list1_igraph[1:38635],
  list2_igraph[38636:54275],
  list3_igraph[54276:68933],
  list4_igraph[68934:70129],
  list5_igraph[70130:118292]
)

length(igraph_list)
igraph_class <- data.frame(unlist(lapply(igraph_list, class)))
igraph_class <- igraph_class$unlist.lapply.igraph_list..class..
table(igraph_class)

#Save
#load("C:\\Users\\asus\\Desktop\\igraph_list_02SET23.RData")
#save(igraph_list, file = "C:\\Users\\asus\\Desktop\\igraph_list_02SET23.RData")

#Delete objects
rm(igraph_list,
   igraph_class,
   list1_igraph,
   list2_igraph,
   list3_igraph,
   list4_igraph,
   list5_igraph
)

################################################################################
#cheddar list

#load("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\new_sept23_dataset\\network_list_cheddar_until_38636.RData")
list1_cheddar <- get("network_list_cheddar")
rm(network_list_cheddar)
list1_cheddar[[1]]
list1_cheddar[[38635]]
#
#load("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\new_sept23_dataset\\network_list_cheddar_until_from_38636_until_54275.RData")
list2_cheddar <- get("network_list_cheddar")
rm(network_list_cheddar)
list2_cheddar[[38636]]
list2_cheddar[[54275]]
#
#load("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\new_sept23_dataset\\network_list_cheddar_until_from_54276_until_68933.RData")
list3_cheddar <- get("network_list_cheddar")
rm(network_list_cheddar)
list3_cheddar[[54276]]
list3_cheddar[[68933]]
#
#load("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\new_sept23_dataset\\network_list_cheddar_from_68934_until_70129.RData")
list4_cheddar <- get("network_list_cheddar")
rm(network_list_cheddar)
list4_cheddar[[68934]]
list4_cheddar[[70129]]
#
#load("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\new_sept23_dataset\\network_list_cheddar_from_70130_until_118292.RData")
list5_cheddar <- get("network_list_cheddar")
rm(network_list_cheddar)
list5_cheddar[[70130]]
list5_cheddar[[118292]]

#Gather all the graph in one list
cheddar_list <- vector(mode='list', length=118292)
for(i in 1:118292){
  if(i >= 1 && i < 38636) cheddar_list[[i]] <- list1_cheddar[[i]]
  if(i >= 38636 && i < 54276) cheddar_list[[i]] <- list1_cheddar[[i]]
  if(i >= 54276 && i < 68934) cheddar_list[[i]] <- list1_cheddar[[i]]
  if(i >= 68934 && i < 70130) cheddar_list[[i]] <- list1_cheddar[[i]]
  if(i >= 70130 && i < 118293) cheddar_list[[i]] <- list1_cheddar[[i]]
  message(i)
}

#Gather all the graph in one list
cheddar_list <- c(
  list1_cheddar[1:38635],
  list2_cheddar[38636:54275],
  list3_cheddar[54276:68933],
  list4_cheddar[68934:70129],
  list5_cheddar[70130:118292]
)

length(cheddar_list)
cheddar_class <- data.frame(unlist(lapply(cheddar_list, class)))
cheddar_class <- cheddar_class$unlist.lapply.cheddar_list..class..
table(cheddar_class)

#Save
#load("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\networks_SET23\\cheddar_list_02SET23.RData")
#save(cheddar_list, file = "C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\networks_SET23\\cheddar_list_02SET23.RData")

#Delete objects
rm(cheddar_list,
   cheddar_class,
   list1_cheddar,
   list2_cheddar,
   list3_cheddar,
   list4_cheddar,
   list5_cheddar
)

#Loading the newly created lists
#load("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\networks_SET23\\igraph_list_02SET23.RData")
#load("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\networks_SET23\\cheddar_list_02SET23.RData")

cheddar_list[[1]]$nodes
