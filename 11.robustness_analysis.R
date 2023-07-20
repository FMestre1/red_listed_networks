################################################################################
#                           Robustness Analysis
################################################################################

#FMestre
#20-07-2023

#Loading the necessary packages
#install.packages("devtools")
library(devtools) 
install_github("FMestre1/fw_package")
library(FWebs)
library(igraph)

#Load the modified functions
source("12.modified_fw_functions.R")

#Loading the list of igraph networks
load("from_cluster/network_list_igraph_2_all_13JUN2023.RData")

network_list_igraph_3 <- vector(mode='list', length=length(network_list_igraph_2))
names(network_list_igraph_3) <- names(network_list_igraph_2)

for(i in 1:length(network_list_igraph_2)) {
  
  network_list_igraph_3[[i]]  <- igraph::upgrade_graph(network_list_igraph_2[[i]])
  message(i)
  
  }

#Save the new, updated, igraphs
#save(network_list_igraph_3, file = "network_list_igraph_3_20JUL.RData")
