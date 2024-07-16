################################################################################
#                            DATASET ON NODE METRICS
################################################################################

#FMestre

#Load packages
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

#Load & Save
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
#                                IGRAPH LIST
################################################################################

#load("igraph_list.RData")

################################################################################
#                                 CHEDDAR LIST
################################################################################

#load("cheddar_list_02SET23.RData")
