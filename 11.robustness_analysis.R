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

#Loading the list of igraph networks
load("network_list_igraph_2.RData")
#network_list_igraph

#Convert to list a list of graph objects
graph_list1 <- convert.to.graph.list(mg1)

#Create a vector with the values for the Intentionality Index (I)
i_index <- seq(from = 0, to = 1, by =0.01)
i_index <- head(i_index,-1)

#Extract one food web as example
fw1 <- graph_list1[[40]]

#Compute the probability to remove each species 
prob_exp <- exponent.removal(fw1, i_index)

#Simulate the extraction of species to evaluate how many primary extinctions are required to have 50% of the total species extinguished
it1 <- iterate(fw_to_attack=fw1, prob_exp, alpha1=50, iter=10, i_index, plot = TRUE)