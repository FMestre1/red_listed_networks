#FMestre
#07-02-2023

library(igraph)
library(cheddar)

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
length(species_names)
length(species_names2)

save(species_names2, file = "species_names2.RData")

################################################################################
#                   SECOND NÚRIA DATASET (network metrics)
################################################################################

dataset2_nuria <- list.files("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_2")
length(dataset2_nuria)

metrics_dataset_2 <- data.frame()
head(metrics_dataset_2)

for(i in 1:length(dataset2_nuria)){
  
  m123 <- read.csv(paste0("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_2\\", dataset2_nuria[i]))
  grid_code <- strsplit(dataset2_nuria[i], ".csv")[[1]]
  metrics_dataset_2 <- rbind(metrics_dataset_2, data.frame(grid_code,m123))
  #
  message(paste0("Did ", i, "!"))
  
}

View(metrics_dataset_2)

#nrow(metrics_dataset_2)
#View(metrics_dataset_2)

metrics_dataset_3 <- unique(metrics_dataset_2)
#View(metrics_dataset_3)
#nrow(metrics_dataset_3)

#Save
save(metrics_dataset_3, file = "metrics_dataset_FINAL.RData")


################################################################################
#                       THIRD NÚRIA DATASET (networks)
################################################################################

#FMestre
#22-05-23

#### RAN IN CLUSTER #### START

network_list <- list()
path3 <- "C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_3"

dataset3_nuria <- list.files("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_3")
#length(dataset3_nuria)

for(i in 1:length(dataset3_nuria)){
  network_list[[i]] <- read.table(paste0(path3, "\\",dataset3_nuria[i]), sep=",", header = 1)
  names(network_list)[i] <- stringr::str_split(dataset3_nuria[i], ".csv")[[1]][1]
  message(i)
}

#Save
#save(network_list, file = "network_list.RData")
#load("network_list.RData")

#### RAN IN CLUSTER #### END

#Coming from the cluster - 31-05-2023
load("network_list_25_may_23.RData")

#Check
#length(network_list)
#names(network_list)

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
  nodes1 <- all_species_status_body_mass_amph_12[all_species_status_body_mass_amph_12$species %in% nodes1,]
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


#Save
save(network_list_igraph, file = "network_list_igraph.RData")
save(network_list_cheddar, file = "network_list_cheddar.RData")

## in the cluster ## END

##### Getting what ran in the cluster #####

library(cheddar)
library(igraph)

load("C:/Users/FMest/Documents/github/red_listed_networks/05_jun_2023/network_list_cheddar_06jun23.RData")
load("C:/Users/FMest/Documents/github/red_listed_networks/05_jun_2023/eur_comm_collection_06jun23.RData")
#
load("C:/Users/FMest/Documents/github/red_listed_networks/05_jun_2023/network_list_igraph_2_ate_72649.RData")
load("C:/Users/FMest/Documents/github/red_listed_networks/05_jun_2023/network_list_igraph_2_de_72650_a_118292.RData")



network_list_igraph_2_ate_72649 <- network_list_igraph_2
rm(network_list_igraph_2)
rm(network_list_igraph_2_ate_72649)
#
network_list_igraph_2_de_72650_a_118292 <- network_list_igraph_2
rm(network_list_igraph_2)
rm(network_list_igraph_2_de_72650_a_118292)
#
list_a <- network_list_igraph_2_ate_72649[1:72649]
#
list_b <- network_list_igraph_2_de_72650_a_118292[72650:118292]

#Finally, concatenate both
network_list_igraph_2 <- c(list_a, list_b)
rm(list_a)
rm(list_b)
#length(network_list_igraph_2)
save(network_list_igraph_2, file = "network_list_igraph_2_all_06JUN2023.RData")
