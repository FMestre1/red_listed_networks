#FMestre
#04-01-2023

library(rredlist) #API - xxx

load("iberian_fw_new_GLOBI_comm_collection_IGRAPH.RData")
load("iberian_fw_MAIORANO_comm_collection_IGRAPH.RData")
load("iberian_fw_new_comm_collection_new_version_28SET2022_IGRAPH.RData")
#
load("iberian_fw_new_comm_collection_new_version_28SET2022_ADJACENCY_MATRICES.RData")
load("iberian_fw_MAIORANO_comm_collection_IGRAPH_ADJACENCY_MATRICES.RData")
load("iberian_fw_new_GLOBI_comm_collection_IGRAPH_ADJACENCY_MATRICES.RData")
#
load("iberian_fw_new_comm_collection_new_version_28SET2022.RData")
load("iberian_fw_new_GLOBI_comm_collection.RData")
load("iberian_fw_MAIORANO_comm_collection.RData")

#Which species?
species <- lapply(iberian_fw_new_comm_collection_new_version_28SET2022_ADJACENCY_MATRICES, colnames)
species <- unique(unlist(species))

rl_version()

rl_sp_country(country = "Portugal",  key = NULL)

################################################################################
#                               Red-listed species 
################################################################################
#Source: https://www.eea.europa.eu/data-and-maps/data/european-red-lists-7

red_listed <- read.csv("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\data\\European_Red_List_2017_December_csv\\European_Red_List_2017_December.csv", 
                       sep = ","
                       )
#head(red_listed)
#View(red_listed)
#unique(red_listed$speciesGroup)

red_listed_mammals_birds <- red_listed[red_listed$speciesGroup == c("Birds", "Mammals"),]
#unique(red_listed_mammals_birds$speciesGroup)
#nrow(red_listed_mammals_birds)

View(red_listed_mammals_birds)
nrow(red_listed_mammals_birds)

#Create vector with genus+species name

full_name <- c()

for(i in 1:nrow(red_listed_mammals_birds)){
  
  #BIRD
  if(red_listed_mammals_birds$speciesGroup[i] == "Birds"){
  genus <- red_listed_mammals_birds$taxonomicRankGenus[i]
  species1 <- red_listed_mammals_birds$taxonomicRankSpecies[i]
  full_name[i] <- paste0(genus, ".", species1)
  }
  
  #MAMMAL
  if(red_listed_mammals_birds$speciesGroup[i] == "Mammals"){
    species2 <- red_listed_mammals_birds$scientificName[i]
    full_name[i] <- stringr::str_replace(species2, " ", ".")
    
  }
  
  
  
  
  
}

#combine with the previous data frame
red_listed_mammals_birds2 <- data.frame(full_name, red_listed_mammals_birds)
rownames(red_listed_mammals_birds2) <- 1:nrow(red_listed_mammals_birds2)
View(red_listed_mammals_birds2)

#Remove unwanted columns
names(red_listed_mammals_birds2)
#

red_listed_mammals_birds3 <- data.frame(
red_listed_mammals_birds2$speciesGroup,
red_listed_mammals_birds2$full_name,
red_listed_mammals_birds2$europeanRegionalRedListCategory,
red_listed_mammals_birds2$endemicToEurope,
red_listed_mammals_birds2$europeanRegionalRedListCategory
)

names(red_listed_mammals_birds3) <- c("group", "full_name", "europeanRegionalRedListCategory",
                                      "endemic_to_europe", "europeanRegionalRedListCategory")

View(red_listed_mammals_birds3)
save(red_listed_mammals_birds3, file = "red_listed_mammals_birds3.RData")

################################################################################