#FMestre
#27-04-2023


#Reptiles database

#reptile_id <- read.csv("C:\\Users\\FMest\\Documents\\0. Posdoc\\CONTRATO\\species_databases\\Life-history trait database of European reptile species\\Species.csv", sep = ";")
#reptile_bs <- read.csv("C:\\Users\\FMest\\Documents\\0. Posdoc\\CONTRATO\\species_databases\\Life-history trait database of European reptile species\\mass.csv", sep = ";")
#head(reptile_id)
#head(reptile_bs)
#unique(reptile_bs$Who)
#reptile_bs_2 <- reptile_bs[reptile_bs$Who %in% c("adult females", "adult males"),]
#View(reptile_bs_2)
#reptile_id_unique <- unique(reptile_bs_2$Species_ID)

#reptile_bs_3 <- data.frame()

#reptile_bs_matched <- merge(x = reptile_id, 
#                            all.x = TRUE, 
#                            y = reptile_bs, 
#                            by.x = "Species.ID", 
#                            by.y = "Species_ID"
#                            )


#lacerta <- read.csv("C:\\Users\\FMest\\Documents\\0. Posdoc\\CONTRATO\\species_databases\\Traits of lizards\\Appendix S1 - Lizard data version 1.0.csv")
#rm(lacerta)


##Why are there repeated species names?

#all_species_status_body_mass_amph_12
nrow(all_species_status_body_mass_amph_12)
nrow(unique(all_species_status_body_mass_amph_12))
all_species_status_body_mass_amph_12[all_species_status_body_mass_amph_12$species == "Ursus maritimus",]

#all_species_status_body_mass_amph_11
nrow(all_species_status_body_mass_amph_11)
nrow(unique(all_species_status_body_mass_amph_11))
all_species_status_body_mass_amph_11[all_species_status_body_mass_amph_11$species == "Ursus maritimus",]

#all_species_status_body_mass_amph_10
nrow(all_species_status_body_mass_amph_10)
nrow(unique(all_species_status_body_mass_amph_10))
all_species_status_body_mass_amph_10[all_species_status_body_mass_amph_10$species == "Ursus maritimus",]

#all_species_status_body_mass_amph_9
nrow(all_species_status_body_mass_amph_9)
nrow(unique(all_species_status_body_mass_amph_9))
all_species_status_body_mass_amph_9[all_species_status_body_mass_amph_9$species == "Ursus maritimus",]

#all_species_status_body_mass_amph_8
nrow(all_species_status_body_mass_amph_8)
nrow(unique(all_species_status_body_mass_amph_8))
all_species_status_body_mass_amph_8[all_species_status_body_mass_amph_8$species == "Ursus maritimus",]

#all_species_status_body_mass_amph_7
nrow(all_species_status_body_mass_amph_7)
nrow(unique(all_species_status_body_mass_amph_7))
all_species_status_body_mass_amph_7[all_species_status_body_mass_amph_7$species == "Ursus maritimus",]

#all_species_status_body_mass_amph_6
nrow(all_species_status_body_mass_amph_6)
nrow(unique(all_species_status_body_mass_amph_6))
all_species_status_body_mass_amph_6[all_species_status_body_mass_amph_6$species == "Ursus maritimus",]

#all_species_status_body_mass_amph_5
nrow(all_species_status_body_mass_amph_5)
nrow(unique(all_species_status_body_mass_amph_5))
all_species_status_body_mass_amph_5[all_species_status_body_mass_amph_5$species == "Ursus maritimus",]

#all_species_status_body_mass_amph_4
nrow(all_species_status_body_mass_amph_4)
nrow(unique(all_species_status_body_mass_amph_4))
all_species_status_body_mass_amph_4[all_species_status_body_mass_amph_4$species == "Ursus maritimus",]

#all_species_status_body_mass_amph_3
nrow(all_species_status_body_mass_amph_3)
nrow(unique(all_species_status_body_mass_amph_3))
all_species_status_body_mass_amph_3[all_species_status_body_mass_amph_3$species == "Ursus maritimus",]

#all_species_status_body_mass_amph_2
nrow(all_species_status_body_mass_amph_2)
nrow(unique(all_species_status_body_mass_amph_2))
all_species_status_body_mass_amph_2[all_species_status_body_mass_amph_2$species == "Ursus maritimus",]

#all_species_status_body_mass_amph
nrow(all_species_status_body_mass_amph)
nrow(unique(all_species_status_body_mass_amph))
all_species_status_body_mass_amph[all_species_status_body_mass_amph$all_species_names == "Ursus maritimus",]

#all_species_status_body_mass
nrow(all_species_status_body_mass)
nrow(unique(all_species_status_body_mass))
all_species_status_body_mass[all_species_status_body_mass$all_species_names == "Ursus maritimus",]

#traits
nrow(traits)
nrow(unique(traits))
traits[traits$traits.species == "Ursus maritimus",]
traits_0[traits_0$species == "Ursus maritimus",]

#all_species_status
nrow(all_species_status)
nrow(unique(all_species_status))
all_species_status[all_species_status$all_species_names == "Ursus maritimus",]

#traits_0
names(traits_0)

traits_0[traits_0$species == "Ursus maritimus",]

#########################################################################################################
#########################################################################################################
#########################################################################################################

# Required libraries
library(terra)

# Load the two shapefiles as spatial objects
#nt_centrality3 <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/centrality_nt_spatial.shp")
#t_centrality3 <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/centrality_t_spatial.shp")
nt_outdegree3 <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/outdegree_nt_spatial.shp")
t_outdegree3 <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/outdegree_t_spatial.shp")

# Extract the attribute values from the shapefiles
#values_map1 <- t_centrality3$centrality
#values_map2 <- nt_centrality3$centrality
values_map1 <- t_outdegree3$outdegree
values_map2 <- nt_outdegree3$outdegree

# Calculate the observed difference between the two maps
obs_diff <- values_map1 - values_map2

# Generate a null distribution of differences using a permutation test

#n_permutations <- 1000  # Number of permutations
n_permutations <- length(nt_centrality3$centrality)  # Number of permutations
null_diff <- rep(NA, n_permutations)

for (i in 1:n_permutations) {
  # Randomly permute the values of the second map
  permuted_values_map2 <- sample(values_map2)
  
  # Calculate the difference between the first map and permuted second map
  perm_diff <- abs(values_map1 - permuted_values_map2)
  
  # Store the difference in the null distribution
  null_diff[i] <- mean(perm_diff, na.rm = TRUE)
  message(i)
}

# Calculate the p-value by comparing the observed difference to the null distribution
p_value <- sum(null_diff >= obs_diff, na.rm=TRUE)/n_permutations
print(p_value)


################################################################################
################################################################################
################################################################################

#FMestre
#18-07-2023

library(terra)

# Read the two rasters
nt_centrality <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_centrality.tif")
t_centrality <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_centrality.tif")

# Calculate the difference between the two rasters
difference <- t_centrality - nt_centrality
plot(difference)
rm(t_centrality, nt_centrality)

# Use a statistical test to determine if the differences are significant
tt1 <- t.test(difference)

# Plot the difference raster
plot(difference)

################################################################################
################################################################################

#Plot side by side
colfunc <- colorRampPalette(c("lightyellow", "tan4", "brown", "black"))

par(mfrow=c(1, 2))
terra::plot(t_indegree, range = c(0, 65), main = "Threatened", col=colfunc(20))
terra::plot(nt_indegree, range = c(0, 65), main = "Not threatened", col=colfunc(20))
#
par(mfrow=c(1, 2))
terra::plot(t_outdegree, range = c(0, 25), main = "Threatened", col=colfunc(20))
terra::plot(nt_outdegree, range = c(0, 25), main = "Not threatened", col=colfunc(20))
#
par(mfrow=c(1, 2))
terra::plot(t_t_level, range = c(0, 3), main = "Threatened", col=colfunc(20))
terra::plot(nt_t_level, range = c(0, 3), main = "Not threatened", col=colfunc(20))
#
par(mfrow=c(1, 2))
terra::plot(t_centrality, range = c(0, 460), main = "Threatened", col=colfunc(20))
terra::plot(nt_centrality, range = c(0, 460), main = "Not threatened", col=colfunc(20))
#
par(mfrow=c(1, 2))
terra::plot(t_ivi, range = c(0, 100), main = "Threatened", col=colfunc(20))
terra::plot(nt_ivi, range = c(0, 100), main = "Not threatened", col=colfunc(20))
#
par(mfrow=c(1, 2))
terra::plot(t_closeness, range = c(0, 1), main = "Threatened", col=colfunc(20))
terra::plot(nt_closeness, range = c(0, 1), main = "Not threatened", col=colfunc(20))

#######################################################################################
#######################################################################################

#Create frequency plots

values_nt_indegree <- terra::values(nt_indegree)
frequency_nt_indegree <- table(values_nt_indegree)
df <- data.frame(value = as.numeric(names(frequency_nt_indegree)), frequency = as.numeric(frequency_nt_indegree))
ggplot(df, aes(x = value, y = frequency)) +
  geom_line() +
  xlab("Value") +
  ylab("Frequency") +
  ggtitle("Frequency of Values in Raster")

######################################################

#########################################################################################

#Using sub-sample of the data frame rows in order to avoid the issue of large samples

# Number of sub-samples (e.g., 100 sub-samples)
num_subsamples <- 100

# Size of each sub-sample (e.g., 50% of the original sample size)
subsample_size <- round(length(outdegree_compare[,2]) * 0.5)

# Create vectors to store sub-sample t-test results
t_statistic <- rep(NA, num_subsamples)
p_value <- rep(NA, num_subsamples)

# Perform sub-sampling and t-test for each sub-sample
for (i in 1:num_subsamples) {
  # Randomly select subsample_size observations from each group
  sub_sample_table <- outdegree_compare[sample(nrow(outdegree_compare), subsample_size), ]
  
  # Perform t-test on the sub-sample
  outdegree_result <- t.test(sub_sample_table[,2], sub_sample_table[,3], paired = TRUE)
  
  # Store t-test results
  t_statistic[i] <- outdegree_result$statistic
  p_value[i] <- outdegree_result$p.value
}

#########################################################################################

#########################################################################################

#Using sub-sample of the data frame rows in order to avoid the issue of large samples

# Number of sub-samples (e.g., 100 sub-samples)
num_subsamples <- 100

# Size of each sub-sample (e.g., 50% of the original sample size)
subsample_size <- round(length(indegree_compare[,2]) * 0.5)

# Create vectors to store sub-sample t-test results
t_statistic <- rep(NA, num_subsamples)
p_value <- rep(NA, num_subsamples)

# Perform sub-sampling and t-test for each sub-sample
for (i in 1:num_subsamples) {
  # Randomly select subsample_size observations from each group
  sub_sample_table <- indegree_compare[sample(nrow(indegree_compare), subsample_size), ]
  
  # Perform t-test on the sub-sample
  indegree_result <- t.test(sub_sample_table[,2], sub_sample_table[,3], paired = TRUE)
  
  # Store t-test results
  t_statistic[i] <- indegree_result$statistic
  p_value[i] <- indegree_result$p.value
}

#########################################################################################


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

################################################################################
################################################################################

#install.packages("diffeR")
library(diffeR)

# Calculate the difference between the two maps
indegree_diff <- differenceMR(t_indegree, nt_indegree)
outdegree_diff <- differenceMR(t_outdegree, nt_outdegree)
t_level_diff <- differenceMR(t_t_level, nt_t_level)
closeness_diff <- differenceMR(t_closeness, nt_closeness)
centrality_diff <- differenceMR(t_centrality, nt_centrality)
ivi_diff <- differenceMR(t_ivi, nt_ivi)

# Plot the difference map
plot(indegree_diff)
plot(outdegree_diff)
plot(t_level_diff)
plot(closeness_diff)
plot(centrality_diff)
plot(ivi_diff)

#Write
writeRaster(indegree_diff, "indegree_diff.tif", overwrite=TRUE)
writeRaster(outdegree_diff, "outdegree_diff.tif", overwrite=TRUE)
writeRaster(t_level_diff, "t_level_diff.tif", overwrite=TRUE)
writeRaster(closeness_diff, "closeness_diff.tif", overwrite=TRUE)
writeRaster(centrality_diff, "centrality_diff.tif", overwrite=TRUE)
writeRaster(ivi_diff, "ivi_diff.tif", overwrite=TRUE)


################################################################################
################################################################################



