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

