#Loading data

load("all_species_status_body_mass_amph_12.RData")
load("fw_list_with_status_aggreg_BS.RData")
load("species_names_fw.RData")

load("species_list_centrality.RData")
load("species_list_presence_absence.RData")
load("trophic_level.RData")

load("network_list_25_may_23.RData")

load("network_list_cheddar_06jun23.RData")
load("network_list_igraph_2.RData")

#load("network_list_igraph.RData")

load("copiar/network_list_cheddar_06jun23.RData")
load("copiar/network_list_igraph_2_all_13JUN2023.RData")
load("copiar/eur_comm_collection_06jun23.RData")

load("igraph_node_attrib_df_summarized_4.RData")


################################################################################
#                       RED LISTED SPECIES
################################################################################

getwd()
setwd("C:\\Users\\fmestre\\red_list")

library(igraph)
library(cheddar)
library(ggplot2)
library(NetIndices)
library(hrbrthemes)
library(dplyr)
library(data.table)

#Code
species_list_centrality <- vector(mode = "list", length = length(species_names_fw))
names(species_list_centrality) <- species_names_fw

species_list_presence_absence <- species_list_centrality
trophic_level <- species_list_centrality

for(i in 1:length(fw_list_with_status_aggreg_BS)){
  
  
  fw_grid <- fw_list_with_status_aggreg_BS[[i]]
  species_fw_grid <- fw_grid$SP_NAME
  nr_species_fw_grid <- length(species_fw_grid)
  
  #species_fw_grid %in% names(species_list_centrality)
  if(nr_species_fw_grid !=0){  
    for(j in 1:nr_species_fw_grid){
      species_list_centrality[species_fw_grid][[j]] <- c(species_list_centrality[species_fw_grid][[j]], fw_grid$centrality[j])
      trophic_level[species_fw_grid][[j]] <- c(trophic_level[species_fw_grid][[j]], fw_grid$TL[j])
      species_list_presence_absence[species_fw_grid][[j]] <- c(species_list_presence_absence[species_fw_grid][[j]], 1)
    }}
  message(i)
}

#save
#save(species_list_centrality, file = "species_list_centrality.RData")
#save(species_list_presence_absence, file = "species_list_presence_absence.RData")
#save(trophic_level, file = "trophic_level.RData")

View(species_list_centrality)
View(species_list_presence_absence)
View(trophic_level)

################################################################################
#                       THIRD NÚRIA DATASET (networks)
################################################################################
#FMestre
#23-05-23

path3 <- "C:\\Users\\fmestre\\red_list\\local-networks"
dataset3_nuria <- list.files("C:\\Users\\fmestre\\red_list\\local-networks")
#length(dataset3_nuria)

network_list <- vector(mode='list', length=length(dataset3_nuria))

for(i in 1:length(dataset3_nuria)){
  network_list[[i]] <- read.table(paste0(path3, "\\",dataset3_nuria[i]), sep=",", header = 1)
  names(network_list)[i] <- stringr::str_split(dataset3_nuria[i], ".csv")[[1]][1]
  message(i)
}

#Save
#save(network_list, file = "network_list_25_may_23.RData")

################################################################################

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
  
  
  #... now the cheddar version
  names(nodes1)[1] <- "node"
  ch1 <- cheddar::Community(nodes = nodes1, properties = list(title= paste0("Interactions in grid ", names(network_list[i]))), trophic.links = edges1)
  network_list_cheddar[[i]] <- ch1
  
  message(i)
  
}

#Save
#save(network_list_igraph, file = "network_list_igraph.RData")
#save(network_list_cheddar, file = "network_list_cheddar.RData")

################################################################################

#upload directly as igraph

path3 <- "C:\\Users\\fmestre\\red_list\\local-networks"

dataset3_nuria <- list.files("C:\\Users\\fmestre\\red_list\\local-networks")
#length(dataset3_nuria)

#IGRAPH #############################################

network_list_igraph <- vector(mode='list', length=length(dataset3_nuria))
length(network_list_igraph)

for(i in 1:length(dataset3_nuria)){
  fw_a <- read.table(paste0(path3, "\\",dataset3_nuria[i]), sep=",", header = 1)
  #names(network_list)[i] <- stringr::str_split(dataset3_nuria[i], ".csv")[[1]][1]
  #message(i)
  
  #fw_a <- network_list[[i]]
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
  
  message(i)
  
}

#Save
#save(network_list_igraph, file = "network_list_igraph.RData")

#CHEDDAR ##########################################

grid_names <- rep(NA, length(dataset3_nuria))

for(i in 1:length(dataset3_nuria)){
  grid_names[i] <- stringr::str_split(dataset3_nuria[i], "\\.")[[1]][1]
}

network_list_cheddar <- vector(mode='list', length=length(dataset3_nuria))
#names(network_list_cheddar) <- names(network_list)

for(i in 1:length(dataset3_nuria)){
  fw_a <- read.table(paste0(path3, "\\",dataset3_nuria[i]), sep=",", header = 1)
  #names(network_list)[i] <- stringr::str_split(dataset3_nuria[i], ".csv")[[1]][1]
  #message(i)
  
  #fw_a <- network_list[[i]]
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
  
  names(nodes1)[1] <- "node"
  
  ch1 <- cheddar::Community(nodes = nodes1, properties = list(title= paste0("Interactions in grid ", grid_names[i])), trophic.links = edges1)
  network_list_cheddar[[i]] <- ch1
  
  message(i)
  
}

#Save
#save(network_list_cheddar, file = "network_list_cheddar_06jun23.RData")

#Create community collection
eur_comm_collection <- cheddar::CommunityCollection(network_list_cheddar)
#class(eur_comm_collection)
#save(eur_comm_collection, file = "eur_comm_collection_06jun23.RData")

#Using toigraph ##################################

path3 <- "C:\\Users\\fmestre\\red_list\\local-networks"
dataset3_nuria <- list.files("C:\\Users\\fmestre\\red_list\\local-networks")

network_list_igraph_2 <- vector(mode='list', length=length(dataset3_nuria))

#Required function
ToIgraph <- function(community, weight=NULL)
{
  if(is.null(TLPS(community)))
  {
    stop('The community has no trophic links')
  }
  else
  {
    tlps <- TLPS(community, link.properties=weight)
    if(!is.null(weight))
    {
      tlps$weight <- tlps[,weight]
    }
    return (graph.data.frame(tlps,
                             vertices=NPS(community),
                             directed=TRUE))
  }
}

for(i in 1:length(network_list_cheddar)){
  
  network_list_igraph_2[[i]] <- ToIgraph(network_list_cheddar[[i]])
  
  names(network_list_igraph_2)[i] <- cheddar::CPS(network_list_cheddar[[i]])$title
  
  message(i)
}

###############################

#Merging 

#load("copiar/network_list_igraph_2_ate_72649.RData")
network_list_igraph_2_ate_72649 <- network_list_igraph_2
rm(network_list_igraph_2)
#
#load("copiar/network_list_igraph_2_de_72650_a_118292.RData")
network_list_igraph_2_de_72650_a_118292 <- network_list_igraph_2
rm(network_list_igraph_2)
#
list_a <- network_list_igraph_2_ate_72649[1:72649]
#
list_b <- network_list_igraph_2_de_72650_a_118292[72650:118292]

#Finally, concatenate both
network_list_igraph_2 <- c(list_a, list_b)
rm(list_a)
rm(list_b)
#length(network_list_igraph_2)

#Save
#save(network_list_igraph_2, file = "network_list_igraph_2_all_13JUN2023.RData")

rm(network_list_igraph_2)
rm(network_list_igraph_2_ate_72649)
rm(network_list_igraph_2_de_72650_a_118292)

#####################################
#  Looking into igraph objects
#####################################

#plot(network_list_cheddar[[200]])
#cheddar::PlotWebByLevel(network_list_cheddar[[200]])
#cheddar::TrophicSpecies(network_list_cheddar[[200]])

cheddar_node_attrib <- vector(mode='list', length=length(network_list_cheddar))

#CHEDDAR
for(i in 1:length(network_list_cheddar)){
  
  cheddar_node_attrib[[i]] <- cbind(cheddar::NPS(network_list_cheddar[[i]]), cheddar::TrophicLevels(network_list_cheddar[[i]]))
  
  names(cheddar_node_attrib)[i] <- cheddar::summary.Community(network_list_cheddar[[i]])$title
  
  message(i)
  
  gc()
  
}

#save(cheddar_node_attrib, file = "cheddar_node_attrib_from_24743_to_30281_with_errors.RData")

#IGRAPH

igraph_node_attrib <- vector(mode='list', length=length(network_list_igraph_2))

for(i in 1:length(igraph_node_attrib)){
  
  igraph_0 <- network_list_igraph_2[[i]]
  
  df_0 <- as.data.frame(do.call(cbind, igraph::vertex.attributes(igraph_0)))  
  
  test.graph.adj_igraph_0 <- get.adjacency(igraph_0, sparse = TRUE)
  
  df_1 <- NetIndices::TrophInd(as.matrix(test.graph.adj_igraph_0))
  
  df_1 <- data.frame(rownames(df_1), df_1)
  
  names(df_1)[1] <- "name"
  
  igraph_node_attrib[[i]] <- merge(df_0, df_1)
  
  message(i)
  
  gc()
  
}

#save(igraph_node_attrib, file = "igraph_node_attrib.RData")

#####################################

#igraph_node_attrib[[25]]

igraph_node_attrib_df <- do.call("rbind", igraph_node_attrib)

#head(igraph_node_attrib_df)
#View(igraph_node_attrib_df)
#nrow(igraph_node_attrib_df)

#save(igraph_node_attrib_df, file = "igraph_node_attrib_df.RData")

# Load ggplot2
# The mtcars dataset is natively available
# head(mtcars)

# A really basic boxplot.
ggplot(igraph_node_attrib_df, aes(x=as.factor(status), y=TL)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  xlab("IUCN Status")

#######

# mtcars dataset is natively available in R
# head(mtcars)

# A basic scatterplot with color depending on Species
#ggplot(igraph_node_attrib_df, aes(x=as.factor(status), y=TL, color=body_size)) + 
#  geom_point(size=6) +
#  theme_ipsum()

igraph_node_attrib_df_summarized <- igraph_node_attrib_df %>%
  group_by(name) %>%
  summarize(average_TL = mean(TL))

igraph_node_attrib_df_summarized <- data.frame(igraph_node_attrib_df_summarized)
names(igraph_node_attrib_df_summarized)
head(igraph_node_attrib_df_summarized)

#save(igraph_node_attrib_df_summarized, file = "igraph_node_attrib_df_summarized.RData")

# Convert data frames to data.table objects
dt1 <- data.table(igraph_node_attrib_df_summarized)
dt2 <- data.table(igraph_node_attrib_df[,1:4])

names(dt1)
names(dt2)

# Merge the data.tables
igraph_node_attrib_df_summarized_3 <- merge(dt1, dt2, by = "name")
#head(igraph_node_attrib_df_summarized_3)
igraph_node_attrib_df_summarized_3 <- as.data.frame(igraph_node_attrib_df_summarized_3)
#save(igraph_node_attrib_df_summarized_3, file = "igraph_node_attrib_df_summarized_3.RData")

igraph_node_attrib_df_summarized_4 <- unique(igraph_node_attrib_df_summarized_3)

# A basic scatterplot with color depending on Species
names(igraph_node_attrib_df_summarized_4)

ggplot(igraph_node_attrib_df_summarized_4, aes(x=as.factor(status), y=average_TL, color=body_size)) + 
  geom_point(size=6) +
  theme_ipsum()

#save(igraph_node_attrib_df_summarized_4, file = "igraph_node_attrib_df_summarized_4.RData")

