#FMestre
#07-02-2023

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

#fw_list[[1]]


################################################################################
#                   SECOND NÚRIA DATASET (network metrics)
################################################################################


dataset2_nuria <- list.files("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_2")
length(dataset2_nuria)

metrics_dataset_2 <- data.frame()

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

