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



###
#rm(t_indegree, nt_indegree,
#   t_outdegree, nt_outdegree,
#   t_t_level, nt_t_level,
#   t_closeness, nt_closeness,
#   t_centrality, nt_centrality,
#   t_ivi, nt_ivi)
###

#indeg_ttest <- t.test(t_indegree, nt_indegree)
#outdeg_ttest <- t.test(t_outdegree, nt_outdegree)
#tl_ttest <- t.test(t_t_level, nt_t_level)
#closeness_ttest <- t.test(t_closeness, nt_closeness)
#centrality_ttest <- t.test(t_centrality, nt_centrality)
#ivi_ttest <- t.test(t_ivi, nt_ivi)

mask_vect <- terra::vect("C:\\Users\\FMest\\Documents\\github\\red_listed_networks\\europeRaster_poly.shp")

# Add the centroid coordinates
centroids <- terra::centroids(mask_vect, TRUE)

# Add the centroid coordinates to the vector
mask_vect <- terra::cbind2(mask_vect, as.data.frame(terra::crds(centroids)))
longitudes <- seq(min(mask_vect$x), max(mask_vect$x), by = 800000)

#plot(terra::subset(mask_vect, mask_vect$x > longitudes[i] & mask_vect$x  < longitudes[i+1]))

i = 3

mask1 <- terra::subset(mask_vect, mask_vect$x > longitudes[i] & mask_vect$x  < longitudes[i+1])
#plot(mask1)
ext1 <- terra::ext(mask1)

#plot(terra::mask(t_indegree, mask1))
mp1 <- terra::mask(t_indegree, mask1)
mp2 <- terra::mask(nt_indegree, mask1)
#
terra::ext(mp1) <- ext1
terra::ext(mp2) <- ext1
#
plot(mp1)
plot(mp2)
#
tt1 <- t.test(mp1, mp2)
#
rm(mp1, mp2)

################################################################################

#29-08-2023

install.packages("spatialEco")
library(spatialEco)
#
indeg_corr <- rasterCorrelation(nt_indegree, t_indegree, s = 3, type = "spearman")
terra::writeRaster(indeg_corr, filename = "indeg_corr.tif", overwrite=TRUE)
plot(indeg_corr)

#O raster ficou esquisito...

################################################################################

#Cohen's d
library(effsize)
indegree_cohens_d_2 <- effsize::cohen.d(indegree_compare_2[complete.cases(indegree_compare_2),][,4],
                                        indegree_compare_2[complete.cases(indegree_compare_2),][,2], paired = TRUE)
outdegree_cohens_d_2 <- effsize::cohen.d(outdegree_compare_2[complete.cases(outdegree_compare_2),][,4],
                                         outdegree_compare_2[complete.cases(outdegree_compare_2),][,2], paired = TRUE)
tl_cohens_d_2 <- effsize::cohen.d(trophic_level_compare_2[complete.cases(trophic_level_compare_2),][,4],
                                  trophic_level_compare_2[complete.cases(trophic_level_compare_2),][,2], paired = TRUE)
closeness_cohens_d_2 <- effsize::cohen.d(closeness_compare_2[complete.cases(closeness_compare_2),][,4],
                                         closeness_compare_2[complete.cases(closeness_compare_2),][,2], paired = TRUE)
centrality_cohens_d_2 <- effsize::cohen.d(centrality_compare_2[complete.cases(centrality_compare_2),][,4],
                                          centrality_compare_2[complete.cases(centrality_compare_2),][,2], paired = TRUE)
ivi_cohens_d_2 <- effsize::cohen.d(ivi_compare_2[complete.cases(ivi_compare_2),][,4],
                                   ivi_compare_2[complete.cases(ivi_compare_2),][,2], paired = TRUE)

round(indegree_cohens_d_2$estimate, 3)
round(outdegree_cohens_d_2$estimate, 3)
round(tl_cohens_d_2$estimate, 3)
round(closeness_cohens_d_2$estimate, 3)
round(centrality_cohens_d_2$estimate, 3)
round(ivi_cohens_d_2$estimate, 3)


################################################################################
#                       Do I have the updated IUCN data
################################################################################

#FMestre
#10-10-2023

#Having the species in the cheddar community
cheddar_list[[4]]$nodes
igraph_list[[4]]

#Check if the IUCN status is the new one...
new_species_iucn[new_species_iucn$scientificName %in% cheddar_list[[4]]$nodes$node,][,c(1,3)]
#... or the old one.
red_listed_3[red_listed_3$full_name %in% cheddar_list[[4]]$nodes$node,][,2:3]

#Yes!

################################################################################

################################################################################

#red_listed_3
#fw_list

fw_list_with_status <- fw_list

#Add IUCN status
for(i in 1:length(fw_list_with_status)){
  
  fw3 <- fw_list_with_status[[i]]
  fw3$SP_NAME <- stringr::str_replace(fw3$SP_NAME, "_", " ")
  sp_fw3 <- fw3$SP_NAME
  
  if(any(red_listed_3$full_name %in% sp_fw3))
  {
    sp_fw3_redList <- red_listed_3[red_listed_3$full_name %in% sp_fw3,]
    fw4 <- merge(fw3, sp_fw3_redList, by.x = "SP_NAME", by.y = "full_name", all.x = TRUE)
    fw_list_with_status[[i]] <- fw4
  }
  
  
  message(i)
  
}

#table(fw_list_with_status_aggreg_BS[[i]]$SP_NAME == "Natrix natrix")
#identical(as.vector(unlist(lapply(fw_list_with_status, nrow))), as.vector(unlist(lapply(fw_list, nrow))))

#Add aggregated IUCN status
fw_list_with_status_aggreg <- fw_list_with_status

#Add threatened/non-threatened
for(i in 1:length(fw_list_with_status)){
  
  fw5 <- fw_list_with_status[[i]]
  
  if(nrow(fw5)!=0){
    
    categories_fw5 <- fw5$europeanRegionalRedListCategory
    #
    categories_fw5 <- stringr::str_replace(categories_fw5, "VU", "threatened")
    categories_fw5 <- stringr::str_replace(categories_fw5, "EN", "threatened")
    categories_fw5 <- stringr::str_replace(categories_fw5, "CR", "threatened")
    #
    categories_fw5 <- stringr::str_replace(categories_fw5, "LC", "non-threatened")
    categories_fw5 <- stringr::str_replace(categories_fw5, "NT", "non-threatened")
    #
    categories_fw5 <- stringr::str_replace(categories_fw5, "DD", "others")
    categories_fw5 <- stringr::str_replace(categories_fw5, "NE", "others")
    categories_fw5 <- stringr::str_replace(categories_fw5, "RE", "others")
    #
    
    fw5 <- data.frame(fw5, categories_fw5)
    names(fw5)[12] <- "aggreg_IUCN"
    
    fw_list_with_status_aggreg[[i]] <- fw5
  }
  message(i)
  
}



################################################################################
#                    RELATING FW METRICS AND FRAGMENTATION
################################################################################
#FMestre
#24-03-2023

#Load packages
library(ggplot2)
library(terra)
library(igraph)
library(taxize)

#The FW metrics
#metrics_dataset_3
#head(metrics_dataset_3)

#Fragmentation or habitat structure (a proxy of it might be "human footprint")
hfootprint <- terra::rast("D:\\fw_space\\human_footprint_1993-2009\\2009\\wildareas-v3-2009-human-footprint.tif")
#plot(hfootprint)
#crs(hfootprint)

#Extract zonal stats for grids
hfootprint_wgs84 <- terra::project(hfootprint, crs(europe))
#plot(hfootprint_wgs84)
#terra::writeRaster(hfootprint_wgs84, filename = "hfootprint_wgs84.tif")
#hfootprint_wgs84 <- terra::rast("hfootprint_wgs84.tif")

#Getting all the values above 50 to NA
myFunction <- function(x){ x[x >= 50] <- NA; return(x)}
hfootprint_wgs84_v2 <- app(hfootprint_wgs84, fun= myFunction)
#plot(hfootprint_wgs84_v2)
#terra::writeRaster(hfootprint_wgs84_v2, filename = "hfootprint_wgs84_v2.tif")
#hfootprint_wgs84_v2 <- terra::rast("hfootprint_wgs84_v2.tif")

hfootprint_wgs84_v2_grid <- terra::zonal(hfootprint_wgs84_v2, fun=mean,  terra::rasterize(europeRaster_poly_wgs84, hfootprint_wgs84_v2, "PageName"))
#head(hfootprint_wgs84_v2_grid)
names(hfootprint_wgs84_v2_grid)[2] <- "hfootprint"
#head(metrics_dataset_3)
#head(hfootprint_wgs84_v2_grid)

metrics_hfootprint <- merge(x = metrics_dataset_3, all.x = TRUE, y = hfootprint_wgs84_v2_grid, 
                            by.x = "grid_code", by.y = "PageName")

#head(metrics_hfootprint)
#plot(metrics_hfootprint$modularity, metrics_hfootprint$hfootprint)

europeRaster_poly_wgs84_2 <- merge(x = europeRaster_poly_wgs84, all.x = TRUE, y = metrics_hfootprint, 
                                   by.x = "PageName", by.y = "grid_code")
#terra::writeVector(europeRaster_poly_wgs84_2, filename = "europeRaster_poly_wgs84_2.shp")
#europeRaster_poly_wgs84_2 <- terra::vect("europeRaster_poly_wgs84_2.shp")

################################################################################
#                              RELATING WITH BODY SIZE
################################################################################

#BS here is a proxy of dispersal distance and resistance to fragmentation
#(Find references for both in the four groups...)

#The FW metrics
#metrics_dataset_3

#Species names
#red_listed_3 #These are the threatened
#not_threat #threatened
#threatened #non-threatened
all_species_names <- stringr::str_replace(species_names2, "_", " ")#These are the names in the list of the metrics per network
all_species_names <- data.frame(all_species_names,1)
#  red_listed_3$full_name

all_species_status <- merge(x = all_species_names, all.x = TRUE, y = red_listed_3, by.x = "all_species_names", by.y = "full_name")
all_species_status <- all_species_status[,c(1,4)]
all_species_status[is.na(all_species_status$europeanRegionalRedListCategory),][,2] <- "not_listed"

agreg_ts <- c()

for(i in 1:nrow(all_species_status)){
  
  cat1 <- all_species_status[i,]$europeanRegionalRedListCategory
  if(cat1 %in% c("LC", "NT")) agreg_ts[i] <- "not_threatened"
  if(cat1 %in% c("VU", "EN", "CR")) agreg_ts[i] <- "threatened"
  if(!(cat1 %in% c("LC", "NT", "VU", "EN", "CR"))) agreg_ts[i] <- "none"
  
}

all_species_status <- data.frame(all_species_status, agreg_ts)

##AmphiBIO - Amphibian Traits ##################################################
amph <- read.csv("C:\\Users\\FMest\\Documents\\0. Posdoc\\CONTRATO\\species_databases\\AmphiBIO\\AmphiBIO_v1.csv")
amph <- data.frame(amph$Species, amph$Body_mass_g)
#View(amph)

#all_species_status_body_mass_amph <- merge(x = all_species_status_body_mass, all.x = TRUE, all.y = FALSE, y = amph, by.x = "all_species_names", by.y = "amph.Species")
all_species_status_body_mass_amph <- merge(x = all_species_status, 
                                           all.x = TRUE, 
                                           all.y = FALSE, 
                                           y = amph, 
                                           by.x = "all_species_names", 
                                           by.y = "amph.Species"
)

#all_species_status_body_mass_amph_2 <- data.frame(all_species_status_body_mass_amph, all_species_status_body_mass$traits.body.mass*1000)
all_species_status_body_mass_amph_2 <- all_species_status_body_mass_amph #04-05-23 - to correct something... (Don´t know what!)

names(all_species_status_body_mass_amph_2) <- c("species", "status", "agreg_ts", "body_mass_g")
#head(all_species_status_body_mass_amph_2)

#Elton Traits - Mammal and Bird Traits ######################################## 
mamm_elton <- read.delim("C:\\Users\\FMest\\Documents\\0. Posdoc\\CONTRATO\\species_databases\\elton_traits\\MamFuncDat.txt")
bird_elton <- read.delim("C:\\Users\\FMest\\Documents\\0. Posdoc\\CONTRATO\\species_databases\\elton_traits\\BirdFuncDat.txt")
#
mamm_elton <- mamm_elton[, c(2, 24)] 
bird_elton <- bird_elton[, c(8, 36)] 
#
bird_mammal_elton <- rbind(bird_elton, mamm_elton)
#head(bird_mammal_elton)

elton_vector <- c()

#for(i in 1:nrow(all_species_status_body_mass_amph_3)){
for(i in 1:nrow(all_species_status_body_mass_amph_2)){
  
  species_row <- all_species_status_body_mass_amph_2$species[i]
  if(species_row %in% bird_mammal_elton$Scientific) {
    elton_vector[i] <- bird_mammal_elton[bird_mammal_elton$Scientific == species_row,]$BodyMass.Value
  }else elton_vector[i] <- NA  
  
}

all_species_status_body_mass_amph_4 <- data.frame(all_species_status_body_mass_amph_2, elton_vector)
#head(all_species_status_body_mass_amph_4)

unified_bs_2 <- c()

for(i in 1:nrow(all_species_status_body_mass_amph_4)){
  
  linha1 <- all_species_status_body_mass_amph_4[i,]
  if(!is.na(linha1$elton_vector)) unified_bs_2[i] <- linha1$elton_vector 
  if(!is.na(linha1$body_mass_g) && is.na(linha1$elton_vector)) unified_bs_2[i] <- linha1$body_mass_g 
  if(is.na(linha1$elton_vector) && is.na(linha1$body_mass_g)) unified_bs_2[i] <- NA
  
  
}

all_species_status_body_mass_amph_5 <- data.frame(all_species_status_body_mass_amph_4, unified_bs_2)
all_species_status_body_mass_amph_6 <- all_species_status_body_mass_amph_5[,-c(4:5)]

#Meiri et al. ##################################################################
meiri_data <- read.csv("C:\\Users\\FMest\\Documents\\0. Posdoc\\CONTRATO\\species_databases\\Meiri_et_al_2021\\meiri_et_al._2021_appendix.csv", sep = ";")
meiri_data <- data.frame(meiri_data$binomial_2020, meiri_data$binomial_.original.files., meiri_data$adult_body_mass..g.)
names(meiri_data) <- c("binomial_2020", "data_binomial_original", "bmass_g")
#head(meiri_data)

unified_bs_3 <- data.frame(rep(NA, nrow(all_species_status_body_mass_amph_6)), rep(NA, nrow(all_species_status_body_mass_amph_6)))
names(unified_bs_3) <- c("bs_A", "bs_B")

for(i in 1:nrow(all_species_status_body_mass_amph_6)){
  
  species_reptiles <- all_species_status_body_mass_amph_6$species[i]
  
  if(species_reptiles %in%  meiri_data$binomial_2020 | species_reptiles %in%  meiri_data$data_binomial_original){
    
    if(length(meiri_data[which(species_reptiles ==  meiri_data$binomial_2020),]$bmass_g)!=0) {
      unified_bs_3$bs_A[i] <- meiri_data[which(species_reptiles ==  meiri_data$binomial_2020),]$bmass_g
    } else unified_bs_3$bs_A[i] <- NA
    
    if(length(meiri_data[which(species_reptiles ==  meiri_data$data_binomial_original),]$bmass_g)!=0) {    
      unified_bs_3$bs_B[i] <- meiri_data[which(species_reptiles ==  meiri_data$data_binomial_original),]$bmass_g
    } else unified_bs_3$bs_B[i] <- NA
  }
}

all_species_status_body_mass_amph_7 <- data.frame(all_species_status_body_mass_amph_6, unified_bs_3$bs_A)
#head(all_species_status_body_mass_amph_7)

unified_bs_4 <- c()

for(i in 1:nrow(all_species_status_body_mass_amph_7)){
  
  #all_species_status_body_mass_amph_7[i,]
  
  if(!is.na(all_species_status_body_mass_amph_7$unified_bs_2[i])) unified_bs_4[i] <- all_species_status_body_mass_amph_7$unified_bs_2[i]
  if(is.na(all_species_status_body_mass_amph_7$unified_bs_2[i])) unified_bs_4[i] <- all_species_status_body_mass_amph_7$unified_bs_3.bs_A[i]
  
}

#head(all_species_status_body_mass_amph_7)

all_species_status_body_mass_amph_8 <- data.frame(all_species_status_body_mass_amph_7[,c(1:3)], unified_bs_4)
#head(all_species_status_body_mass_amph_8)

missing_species_bs <- all_species_status_body_mass_amph_8[is.na(all_species_status_body_mass_amph_8$unified_bs_4),]$species

#I have to check this with the synonyms, resorting to taxize
#Gather all the BS info from previously used sources

#taxize::use_entrez()
#usethis::edit_r_environ()
#ENTREZ_KEY='fafd2118668fc6bacdf37d11c7c1885f5308'#mykey - have to reload R
#all_species_status_body_mass_amph_6[is.na(all_species_status_body_mass_amph_6$unified_bs_2),][,1][1]

syn_list <- rep(NA, length(missing_species_bs))

for(i in 1:length(missing_species_bs)){
  
  species1 <- missing_species_bs[i]
  
  try(df1 <- get_gbifid(species1, rank="species"),
      silent = TRUE)
  
  syn <- id2name(id = df1[1], db = "gbif")
  syn <- syn[[1]]$name
  syn_list[i] <- syn
  
  #delete
  if(exists("df1"))rm(df1)
  if(exists("species1"))rm(species1)
  if(exists("syn"))rm(syn)
  
  message("########## Did ", i, "! ##########")
  
}
#
#any(syn_list %in% amph$amph.Species)
#any(syn_list %in% meiri_data$binomial_2020)
#any(syn_list %in% meiri_data$data_binomial_original)

synonym_table <- data.frame(missing_species_bs, syn_list, syn_list == missing_species_bs, NA)
names(synonym_table) <- c("original_name", "synonym", "match", "bs")
#head(synonym_table)

for(i in 1:nrow(synonym_table)){
  
  synonym_row <- synonym_table[i,]
  
  synonym2 <- synonym_row$synonym
  
  if (synonym2 %in% amph$amph.Species) bs_syn <- amph[which(synonym2 == amph$amph.Species),]$amph.Body_mass_g
  if (synonym2 %in% meiri_data$binomial_2020) bs_syn <- meiri_data[meiri_data$binomial_2020==synonym2,]$bmass_g
  if (synonym2 %in% meiri_data$data_binomial_original) bs_syn <- meiri_data[meiri_data$data_binomial_original==synonym2,]$bmass_g
  
  if(exists("bs_syn")) synonym_table[i,4] <- bs_syn
  if(!exists("bs_syn")) synonym_table[i,4] <- NA
  
  if(exists("bs_syn"))rm(bs_syn)
  
}

#finally get the body mass of these synonyms into the main table 

all_species_status_body_mass_amph_9 <- data.frame(all_species_status_body_mass_amph_8, NA)
names(all_species_status_body_mass_amph_9)[5] <- "synonym"

for(i in 1:nrow(all_species_status_body_mass_amph_9)){
  
  species_to_evaluate <- all_species_status_body_mass_amph_9[i,]$species
  
  if(species_to_evaluate %in% synonym_table$original_name){
    
    all_species_status_body_mass_amph_9[i,]$synonym <- synonym_table[which(synonym_table$original_name == species_to_evaluate),]$synonym
    all_species_status_body_mass_amph_9[i,]$unified_bs_4 <-synonym_table[which(synonym_table$original_name == species_to_evaluate),]$bs
    
    
  }
  
}

all_species_status_body_mass_amph_9[is.na(all_species_status_body_mass_amph_9$synonym),]$synonym <- "NA"
#nrow(all_species_status_body_mass_amph_9)
#View(all_species_status_body_mass_amph_9)

names(all_species_status_body_mass_amph_9)[4] <- "body_size"
#head(all_species_status_body_mass_amph_9)

################################################################################
#                            CONDUCTING THE ANALYSIS                           
################################################################################
#FMestre
#28-04-2023

#Having the:
#modularity & human footprint (proxy of habitat disturbance)
modularity_hfootprint <- as.data.frame(europeRaster_poly_wgs84_2)
#modularity_hfootprint

#body size (proxy of dispersal) & threat level
all_species_status_body_mass_amph_10 <- all_species_status_body_mass_amph_9
all_species_status_body_mass_amph_10$status <- as.factor(all_species_status_body_mass_amph_10$status)
all_species_status_body_mass_amph_10$agreg_ts <- as.factor(all_species_status_body_mass_amph_10$agreg_ts)

all_species_status_body_mass_amph_12 <- unique(all_species_status_body_mass_amph_10) #to correct an error 

names(all_species_status_body_mass_amph_12)[4] <- "body_size"
#head(all_species_status_body_mass_amph_12)

boxplot(body_size~status,data=all_species_status_body_mass_amph_12, main="Body size per IUCN status",
        xlab="status", ylab="body size")
#
aggregate(x = all_species_status_body_mass_amph_12$body_size,      
          by = list(all_species_status_body_mass_amph_12$status),              
          FUN = mean)  

boxplot(body_size~agreg_ts,data=all_species_status_body_mass_amph_12, main="Body size per aggregated IUCN status",
        xlab="status", ylab="body size")
#
aggregate(x = all_species_status_body_mass_amph_12$body_size,      
          by = list(all_species_status_body_mass_amph_12$agreg_ts),              
          FUN = mean)

#Adding BS information
fw_list_with_status_aggreg_BS <- fw_list_with_status_aggreg

for(i in 1:length(fw_list_with_status_aggreg_BS)){
  
  fw1 <- fw_list_with_status_aggreg_BS[[i]]
  
  #length(fw1$SP_NAME)
  #which(all_species_status_body_mass_amph_11$species %in% fw1$SP_NAME)
  #all_species_status_body_mass_amph_11[fw1$SP_NAME == all_species_status_body_mass_amph_11$species,]
  fw1_bs <- all_species_status_body_mass_amph_12[which(all_species_status_body_mass_amph_12$species %in% fw1$SP_NAME),]
  
  
  #all_species_status_body_mass_amph_12[which(all_species_status_body_mass_amph_12$species %in% fw1$SP_NAME),]
  
  if(nrow(fw1) != nrow(fw1_bs))
  {
    if (any(!fw1$SP_NAME %in% fw1_bs$species)) fw1_bs2 <- data.frame(fw1[!fw1$SP_NAME %in% fw1_bs$species,][,c(1,10,12)], NA, NA)
    if (any(fw1$SP_NAME %in% fw1_bs$species)) fw1_bs2 <- data.frame(matrix(ncol = 5, nrow = 0))
    
    names(fw1_bs2) <- names(fw1_bs)
    fw1_bs <- rbind(fw1_bs, fw1_bs2)
  }
  
  #fw1$SP_NAME[!fw1$SP_NAME %in% fw1_bs$species]
  
  if(nrow(fw1_bs)!=0) 
  {
    fw2 <- merge(x =fw1, y = fw1_bs, by.x = "SP_NAME", by.y = "species", all.x = TRUE)
    fw2 <- fw2[,-c(10,12)]
    names(fw2)[12] <- "aggregated_status"
    names(fw2)[13] <- "body_size"
  }
  
  if(nrow(fw1_bs)==0) fw2 <- fw1
  
  fw_list_with_status_aggreg_BS[[i]] <- fw2
  
  message(i)
  
}

#Save
#save(fw_list_with_status_aggreg_BS, file = "fw_list_with_status_aggreg_BS.RData")

#### Two relevant questions:
#1. Why is there a gradient SW-NE in the centrality?
#2. Why is this gradient more intense in the threatened species than in the non-threatened?

############# Create a species richness vector ############# START
#Get grids
europeRaster_poly_wgs84
#plot(europeRaster_poly_wgs84)
#nrow(europeRaster_poly_wgs84)

species_richness_df <- as.data.frame(matrix(nrow = nrow(europeRaster_poly_wgs84), ncol = 2))
names(species_richness_df) <- c("grid_code", "sp_richness")

#Extract the number of species
for(i in 1:nrow(europeRaster_poly_wgs84)){
  
  grid_name <- as.data.frame(europeRaster_poly_wgs84[i,])[,1]
  grid_richness <- nrow(fw_list_with_status_aggreg_BS[[grid_name]])
  
  species_richness_df[i,1] <- grid_name
  if(!is.null(grid_richness)) species_richness_df[i,2] <- grid_richness  
  
  message(i)
  
}

sp_richness <- terra::merge(x=europeRaster_poly_wgs84, y=species_richness_df, by.x = "PageName", by.y = "grid_code")
#crs(sp_richness)
#writeVector(sp_richness, filename ="sp_richness.shp", overwrite=TRUE, filetype = "ESRI Shapefile")

############# Create a species richness vector ############# END

############# Relating species richness and the FW metrics ############# START

# 1. Getting data

#Centrality - T
centrality_t_spatial <- terra::vect("centrality_t_spatial.shp")

#Centrality - NT
centrality_nt_spatial <- terra::vect("centrality_nt_spatial.shp")

#IVI - T
ivi_t_spatial <- terra::vect("ivi_t_spatial_second_version.shp")

#IVI - NT
ivi_nt_spatial <- terra::vect("ivi_nt_spatial_second_version.shp")

#Closeness - T
closeness_t_spatial <- terra::vect("closeness_t_spatial.shp")

#Closeness - NT
closeness_nt_spatial <- terra::vect("closeness_nt_spatial.shp")

#In-degree - T
indegree_t_spatial <- terra::vect("indegree_t_spatial.shp")

#In-degree - NT
indegree_nt_spatial <- terra::vect("indegree_nt_spatial.shp")

#Out-degree - T
outdegree_t_spatial <- terra::vect("outdegree_t_spatial.shp")

#Out-degree - NT
outdegree_nt_spatial <- terra::vect("outdegree_nt_spatial.shp")

#Trophic level - T
tl_t_spatial <- terra::vect("tl_t_spatial.shp")

#Trophic level - NT
tl_nt_spatial <- terra::vect("tl_nt_spatial.shp")

#Proportion of threatened species
proportion_spatial <- terra::vect("proportion_spatial.shp")

# 2. Test significant relationships

#Centrality - T

species_richness_vs_centrality_T <- merge(x = sp_richness, y = centrality_t_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_centrality_T <- as.data.frame(species_richness_vs_centrality_T)
#
wilcox.test(species_richness_vs_centrality_T$sp_richness, 
            species_richness_vs_centrality_T$centrality, 
            paired = TRUE)

#Centrality - NT
#length(sp_richness$sp_richness)
#length(centrality_nt_spatial$centrality)
#
species_richness_vs_centrality_NT <- merge(x = sp_richness, y = centrality_nt_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_centrality_NT <- as.data.frame(species_richness_vs_centrality_NT)
#
wilcox.test(species_richness_vs_centrality_NT$sp_richness, 
            species_richness_vs_centrality_NT$centrality, 
            paired = TRUE)

#IVI - T
#length(sp_richness$sp_richness)
#length(ivi_t_spatial$ivi)
#
species_richness_vs_ivi_T <- merge(x = sp_richness, y = ivi_t_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_ivi_T <- as.data.frame(species_richness_vs_ivi_T)
#
wilcox.test(species_richness_vs_ivi_T$sp_richness, 
            species_richness_vs_ivi_T$ivi, 
            paired = TRUE)

#IVI - NT
#length(sp_richness$sp_richness)
#length(ivi_nt_spatial$ivi)
#
species_richness_vs_ivi_NT <- merge(x = sp_richness, y = ivi_nt_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_ivi_NT <- as.data.frame(species_richness_vs_ivi_NT)
#
wilcox.test(species_richness_vs_ivi_NT$sp_richness, 
            species_richness_vs_ivi_NT$ivi, 
            paired = TRUE)

#Closeness - T
#length(sp_richness$sp_richness)
#length(closeness_t_spatial$closeness)
#
species_richness_vs_closeness_T <- merge(x = sp_richness, y = closeness_t_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_closeness_T <- as.data.frame(species_richness_vs_closeness_T)
#
wilcox.test(species_richness_vs_closeness_T$sp_richness, 
            species_richness_vs_closeness_T$closeness, 
            paired = TRUE)

#Closeness - NT
#length(sp_richness$sp_richness)
#length(closeness_nt_spatial$closeness)
#
species_richness_vs_closeness_NT <- merge(x = sp_richness, y = closeness_nt_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_closeness_NT <- as.data.frame(species_richness_vs_closeness_NT)
#
wilcox.test(species_richness_vs_closeness_NT$sp_richness, 
            species_richness_vs_closeness_NT$closeness, 
            paired = TRUE)

#In-degree - T
#length(sp_richness$sp_richness)
#length(indegree_t_spatial$indegree)
#
species_richness_vs_indegree_T <- merge(x = sp_richness, y = indegree_t_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_indegree_T <- as.data.frame(species_richness_vs_indegree_T)
#
wilcox.test(species_richness_vs_indegree_T$sp_richness, 
            species_richness_vs_indegree_T$indegree, 
            paired = TRUE)

#In-degree - NT
#length(sp_richness$sp_richness)
#length(indegree_nt_spatial$indegree)
#
species_richness_vs_indegree_NT <- merge(x = sp_richness, y = indegree_nt_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_indegree_NT <- as.data.frame(species_richness_vs_indegree_NT)
#
wilcox.test(species_richness_vs_indegree_NT$sp_richness, 
            species_richness_vs_indegree_NT$indegree, 
            paired = TRUE)

#Out-degree - T
#length(sp_richness$sp_richness)
#length(outdegree_t_spatial$outdegree)
#
species_richness_vs_outdegree_T <- merge(x = sp_richness, y = outdegree_t_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_outdegree_T <- as.data.frame(species_richness_vs_outdegree_T)
#
wilcox.test(species_richness_vs_outdegree_T$sp_richness, 
            species_richness_vs_outdegree_T$outdegree, 
            paired = TRUE)

#Out-degree - NT
#length(sp_richness$sp_richness)
#length(outdegree_nt_spatial$outdegree)
#
species_richness_vs_outdegree_NT <- merge(x = sp_richness, y = outdegree_nt_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_outdegree_NT <- as.data.frame(species_richness_vs_outdegree_NT)
#
wilcox.test(species_richness_vs_outdegree_NT$sp_richness, 
            species_richness_vs_outdegree_NT$outdegree, 
            paired = TRUE)

#Trophic level - T
#length(sp_richness$sp_richness)
#length(tl_t_spatial$trophic_le)
#
species_richness_vs_tl_T <- merge(x = sp_richness, y = tl_t_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_tl_T <- as.data.frame(species_richness_vs_tl_T)
#
wilcox.test(species_richness_vs_tl_T$sp_richness, 
            species_richness_vs_tl_T$trophic_le, 
            paired = TRUE)

#Trophic level - NT
#length(sp_richness$sp_richness)
#length(tl_nt_spatial$trophic_le)
#
species_richness_vs_tl_NT <- merge(x = sp_richness, y = tl_nt_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_tl_NT <- as.data.frame(species_richness_vs_tl_NT)
#
wilcox.test(species_richness_vs_tl_NT$sp_richness, 
            species_richness_vs_tl_NT$trophic_le, 
            paired = TRUE)

#Proportion of threatened species

#length(sp_richness$sp_richness)
#length(proportion_spatial$proportion)

species_richness_vs_proportion <- merge(x = sp_richness, y = proportion_spatial, by.x = "PageName", by.y = "PageName", all = TRUE)
species_richness_vs_proportion <- as.data.frame(species_richness_vs_proportion)

wilcox.test(species_richness_vs_proportion$sp_richness, 
            species_richness_vs_proportion$proportion, 
            paired = TRUE)

############# Relating species richness and the FW metrics ############# END

############# Create a data frame with: ############# START
#species
#Coef. variation of the centrality
#Number of occurrence grids
#Trophic positions

#having...
#all_species_status_body_mass_amph_12

#Had to do this in the cluster... and bring it back here:#######################
#First saving the required files to take to the cluster:
#save(all_species_status_body_mass_amph_12, file = "all_species_status_body_mass_amph_12.RData")
#save(fw_list_with_status_aggreg_BS, file = "fw_list_with_status_aggreg_BS.RData")
#After running the code bellow in the cluster, bring everything here:

#START - Ran in the cluster #####

#To obtain species average centrality, trophic level and presence/absence 

species_names_fw2 <- unique(stringr::str_replace(species_names123, "_", " "))
#identical(species_names_fw, species_names_fw2)

species_names_fw <- stringr::str_replace(species_names, "_", " ")
species_names_fw <- unique(species_names_fw)
#save(species_names_fw, file="species_names_fw.RData")

species_list_centrality <- vector(mode = "list", length = length(species_names_fw))
names(species_list_centrality) <- species_names_fw

#save(species_list_centrality, file = "species_list_centrality.RData")

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

#END - Ran in the cluster #####

#Coming from cluster

#LOAD
load("all_species_status_body_mass_amph_12.RData")
load("fw_list_with_status_aggreg_BS.RData")
load("species_names_fw.RData")
#
load("species_list_centrality.RData")
load("species_list_presence_absence.RData")
load("trophic_level.RData")

###

species_list_centrality
species_list_presence_absence
trophic_level

#DF with the CV of centrality

cv_centrality <- data.frame(species_names_fw, NA)
names(cv_centrality) <- c("species", "cv")  
#View(cv_centrality)

for(i in 1:nrow(cv_centrality)){
  
  spe1 <- cv_centrality[i,1]
  cent <- species_list_centrality[[spe1]]
  
  if(!is.null(cent)) cv_centrality[i,2] <-  sd(cent)/mean(cent)
  
  message(i)
  
}

nr_grids_with_presence <- data.frame(species_names_fw, NA)
names(nr_grids_with_presence) <- c("species", "grids_with_presence")  
#View(nr_grids_with_presence)

for(i in 1:nrow(nr_grids_with_presence)){
  
  spe2 <- nr_grids_with_presence[i,1]
  pres2 <- species_list_presence_absence[[spe2]]
  
  if(!is.null(pres2)) nr_grids_with_presence[i,2] <-  sum(pres2)
  
  message(i)
  
}

trophic_level_vector <- data.frame(species_names_fw, NA)
names(trophic_level_vector) <- c("species", "tl")  
#View(trophic_level_vector)

for(i in 1:nrow(trophic_level_vector)){
  
  spe3 <- trophic_level_vector[i,1]
  pres3 <- trophic_level[[spe3]]
  
  if(!is.null(pres3)) trophic_level_vector[i,2] <-  mean(pres3)
  
  message(i)
  
}

#
cv_presence_tl <-
  data.frame(
    cv_centrality,
    nr_grids_with_presence,
    trophic_level_vector
  )

cv_presence_tl <- cv_presence_tl[,-c(3,5)]
#head(cv_presence_tl)

#View(all_species_status_body_mass_amph_12)

for(i in 1:nrow(cv_presence_tl)){
  
  spe4 <- cv_presence_tl$species[i]
  dt4 <- all_species_status_body_mass_amph_12[all_species_status_body_mass_amph_12$species == spe4,]
  if(nrow(dt4) !=0){
    cv_presence_tl[i,5] <- dt4$status #STATUS
    cv_presence_tl[i,6] <- dt4$agreg_ts #AGGREG STATUS
    cv_presence_tl[i,7] <- dt4$body_size #BS
    cv_presence_tl[i,8] <- dt4$synonym #SYNONYM
  }
}

names(cv_presence_tl)[5:8] <- c("status", "aggreg_status", "body_size", "synonym")
#head(cv_presence_tl)

plot(cv_presence_tl$body_size, cv_presence_tl$cv)
plot(cv_presence_tl$body_size, cv_presence_tl$grids_with_presence)
plot(cv_presence_tl$grids_with_presence, cv_presence_tl$cv)
plot(cv_presence_tl$tl, cv_presence_tl$body_size)

#save(cv_presence_tl, file = "cv_presence_tl.RData")
#load("cv_presence_tl.RData")

cv_presence_tl$status <- as.factor(cv_presence_tl$status)
cv_presence_tl$aggreg_status <- as.factor(cv_presence_tl$aggreg_status)
#
cv_presence_tl_version_2 <- cv_presence_tl[complete.cases(cv_presence_tl$cv),]
#head(cv_presence_tl_version_2)

#save(cv_presence_tl_version_2, file = "cv_presence_tl_version_2.RData")
#load("cv_presence_tl_version_2.RData")

###

ggplot(cv_presence_tl_version_2, aes(x=status, y=cv)) + 
  geom_boxplot(fill="slateblue", alpha=0.2)

ggplot(cv_presence_tl_version_2, aes(x=aggreg_status, y=cv)) + 
  geom_boxplot(fill="slateblue", alpha=0.2)

##

ggplot(cv_presence_tl_version_2, aes(x=status, y=body_size)) + 
  geom_boxplot(fill="slateblue", alpha=0.2)

ggplot(cv_presence_tl_version_2, aes(x=aggreg_status, y=body_size)) + 
  geom_boxplot(fill="slateblue", alpha=0.2)

############# Create a data frame with: ############# END

#head(cv_presence_tl) #Species info
#head(metrics_dataset_3) #Network Metrics

#fw_list[[1]] #species metrics per network
#fw_list_with_status_aggreg_BS[[1]] #species metrics per network with IUCN status and BS

#Get the latitude and longitude 
#head(europeRaster_poly_wgs84_coords)

#Merging network metrics and coordinates
MERGE_metrics_coords <- merge(x=metrics_dataset_3,
                              y=europeRaster_poly_wgs84_coords,
                              by.x="grid_code",
                              by.y="PageName",
                              all = TRUE)

#head(MERGE_metrics_coords)
plot(MERGE_metrics_coords$y, MERGE_metrics_coords$C)
plot(MERGE_metrics_coords$y, MERGE_metrics_coords$indegree)
plot(MERGE_metrics_coords$y, MERGE_metrics_coords$outdegree)
plot(MERGE_metrics_coords$y, MERGE_metrics_coords$top)
plot(MERGE_metrics_coords$x, MERGE_metrics_coords$modularity)

#Merging network's average centrality in threatened species with coordinates

#head(species_richness_vs_centrality_NT)
MERGE_centrality_NT_coords <- merge(x=species_richness_vs_centrality_NT,
                                    y=europeRaster_poly_wgs84_coords,
                                    by.x="PageName",
                                    by.y="PageName",
                                    all = TRUE)


#head(species_richness_vs_centrality_T)
MERGE_centrality_T_coords <- merge(x=species_richness_vs_centrality_T,
                                   y=europeRaster_poly_wgs84_coords,
                                   by.x="PageName",
                                   by.y="PageName",
                                   all = TRUE)

plot(MERGE_centrality_NT_coords$x, MERGE_centrality_NT_coords$centrality)
plot(MERGE_centrality_NT_coords$x, MERGE_centrality_NT_coords$centrality)


par(mfrow=c(1,2))
plot(MERGE_centrality_NT_coords$y, MERGE_centrality_NT_coords$centrality)
plot(MERGE_centrality_T_coords$y, MERGE_centrality_T_coords$centrality)


################################################################################
#   NEWER DATASET ON HABITAT FRAGMENTATION
################################################################################
#FMestre
#17-07-2023

# From the paper "Global forest fragmentation change from 2000 to 2020":
#https://www.nature.com/articles/s41467-023-39221-x
FFI2000 <- terra::rast("D:\\Dados Habitats e Usos Solo\\Global forest fragmentation change from 2000 to 2020\\FFI2000.tif")
#plot(FFI2000)
#
FFI2020 <- terra::rast("D:\\Dados Habitats e Usos Solo\\Global forest fragmentation change from 2000 to 2020\\FFI2020.tif")
#plot(FFI2020)
#
#rm(FFI2000)
#rm(FFI2020)

################################################################################
#                                script 6 clustering
################################################################################


#FMestre
#03-05-2023

#Paper: https://doi.org/10.1111/ele.12002

#library(bipartite)
#library(combinat)
#library(factoextra)
#library(terra)

########Incomplete code... bellow!

#df_S_combined <- as.data.frame(matrix(ncol = length(fw_list_with_status_aggreg_BS), 
#                                      nrow = length(fw_list_with_status_aggreg_BS))
#)


#rownames(df_S_combined) <- names(fw_list_with_status_aggreg_BS)
#colnames(df_S_combined) <- names(fw_list_with_status_aggreg_BS)

#df_OS_combined <- df_S_combined 
#df_WN_combined <- df_S_combined 
#df_ST_combined <- df_S_combined 

#grid_names <- names(iberian_fw_new_comm_collection_new_version_28SET2022)

#matrix_comparisons_combined <- data.frame(t(combn(grid_names, 2, simplify = TRUE)))
#nrow(matrix_comparisons_combined)

##

#Save these i
#save_i_combined <- seq(from = 1, to = nrow(matrix_comparisons_combined), by = round(nrow(matrix_comparisons_combined)/1000, 0))
#save_i_combined <- c(save_i_combined, nrow(matrix_comparisons_combined))

#for(i in 1:nrow(matrix_comparisons_combined)){

#  row1 <-  matrix_comparisons_combined[i,]

#  nt1 <- as.character(row1[1])
#  nt2 <- as.character(row1[2])

#  adj1 <- as.matrix(combined_dataset_globi_maiorano_tg_ADJACENCY_MATRIX[[nt1]])
#  adj2 <- as.matrix(combined_dataset_globi_maiorano_tg_ADJACENCY_MATRIX[[nt2]])

#  try(out1 <- betalinkr(webs2array(list(adj1=adj1, adj2=adj2)), 
#                        partitioning="commondenom", 
#                        partition.st=FALSE
#  ))

#  df_S_combined[nt1,nt2] <- as.numeric(out1["S"])
#  df_OS_combined[nt1,nt2] <- as.numeric(out1["OS"])
#  df_WN_combined[nt1,nt2] <- as.numeric(out1["WN"])
#  df_ST_combined[nt1,nt2] <- as.numeric(out1["ST"])

#  message(round((i*100)/nrow(matrix_comparisons_combined), 2))

#  if(i !=0) rm(out1)

#  if(any(i == save_i_combined)) {
#    save(df_S_combined, file = paste0("matrix_out_combined/df_S_combined_16dez23_", i, ".RData"))
#    save(df_OS_combined, file = paste0("matrix_out_combined/df_OS_combined_16dez23_",i, ".RData"))
#    save(df_WN_combined, file = paste0("matrix_out_combined/df_WN_combined_16dez23_", i, ".RData"))
#    save(df_ST_combined, file = paste0("matrix_out_combined/df_ST_combined_16dez23_", i, ".RData"))

#  }

#}#Ran this at the cluster 

########Incomplete code... up!

################################################################################
################################################################################
################################################################################
#FMestre
#05-04-2023

library(terra)
library(factoextra)

europe <- terra::vect("C:/Users/FMest/Documents/0. Artigos/IUCN_networks/shapefiles/Europe/Europe.shp")
crs(europe)

#Add modularity to the grid ####################################################

#where's the grid? ##########
europeRaster <- terra::rast(x="C:/Users/FMest/Documents/github/red_listed_networks/mask10k-20230214T144658Z-001/mask10k/reference_grid_10km.img")
cells_info <- foreign::read.dbf(file = "C:/Users/FMest/Documents/github/red_listed_networks/mask10k-20230214T144658Z-001/mask10k/reference_grid_10km.img.vat.dbf")
#head(cells_info)
#nrow(cells_info)
#To vector
europeRaster_poly_0 <- terra::as.polygons(europeRaster, values = TRUE, extent=FALSE)
europeRaster_poly <- merge(europeRaster_poly_0, cells_info)
plot(europeRaster_poly)
class(europeRaster_poly)
#
plot(europeRaster_poly_0)
class(europeRaster_poly_0)

#where's the modularity (and other metrics)? ##########

path2 <- "C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\data\\data_nuria_2\\"

files_folder2 <- list.files(path2)
#length(files_folder2)

fw_list2 <- data.frame(rep(NA, length(files_folder2)), 
                       rep(NA, length(files_folder2)),
                       rep(NA, length(files_folder2)),
                       rep(NA, length(files_folder2)), 
                       rep(NA, length(files_folder2)),
                       rep(NA, length(files_folder2))
)
colnames(fw_list2) <- c("network", "modularity", "connectance", "basal", "top", "intermediate")
head(fw_list2)

#
for(i in 1:length(files_folder2)){
  fw_list2[i,2] <- read.table(paste0(path2, files_folder2[i]), sep=",", header = 1)$modularity
  fw_list2[i,1] <- stringr::str_split(files_folder2[i], ".csv")[[1]][1]
  fw_list2[i,3] <- read.table(paste0(path2, files_folder2[i]), sep=",", header = 1)$C
  fw_list2[i,4] <- read.table(paste0(path2, files_folder2[i]), sep=",", header = 1)$basal
  fw_list2[i,5] <- read.table(paste0(path2, files_folder2[i]), sep=",", header = 1)$top
  fw_list2[i,6] <- read.table(paste0(path2, files_folder2[i]), sep=",", header = 1)$intermediate
  
  message(i)
}
#

#Let's bring them together ... ##########
modularity_spatial <- merge(europeRaster_poly_0, fw_list2, by.x = "PageName", by.y = "network")
modularity_spatial_wgs84 <- terra::project(modularity_spatial, europe)
plot(modularity_spatial_wgs84)

#Write vector
#writeVector(modularity_spatial_wgs84, filename ="modularity_spatial_wgs84.shp", overwrite=TRUE, filetype = "ESRI Shapefile")

#... and plot it... ##########
#plot(modularity_spatial_wgs84)

#Clustering ####################################################################

modularity_spatial_2 <- as.data.frame(modularity_spatial)
modularity_spatial_2$modularity

# Compute the WSS for different numbers of clusters
#wss <- numeric(10)
#for (k in 1:10) {
#  kmeans_model <- kmeans(modularity_spatial_2$modularity, centers = k)
#  wss[k] <- kmeans_model$tot.withinss
#}

# Plot the elbow curve

modularity_spatial_2_scaled <- scale(modularity_spatial_2[sample(nrow(modularity_spatial_2), round(nrow(modularity_spatial_2)/5)), 4:8])
#?fviz_nbclust

#Next code takes too long...
#fviz_nbclust(modularity_spatial_2_scaled, kmeans, method = "silhouette", k.max = 10, verbose =TRUE)

# Initialize the k-means model with the number of clusters you want
kmeans_model <- kmeans(modularity_spatial_2[,4:8], centers = 2)

# Get the cluster labels and centroids
labels1 <- kmeans_model$cluster
centroids1 <- kmeans_model$centers
#
labels1
centroids1

length(labels1)
nrow(europeRaster_poly_0)

#Create df
clusters_df <- data.frame(modularity_spatial_2$PageName, labels1)
#View(clusters_df)
names(clusters_df) <- c("grid", "clusters")
#Save
write.csv(clusters_df, file = "clusters_df.csv")


################################################################################
#   script 7 geeting info from cluster
################################################################################

#FMestre
#21-06-2023

library(ggplot2)

load("igraph_node_attrib_df_summarized.RData")
load("igraph_node_attrib_df.RData")

head(igraph_node_attrib_df_summarized)
head(igraph_node_attrib_df)

igraph_node_attrib_df_summarized_2 <- merge(x=igraph_node_attrib_df_summarized, y=igraph_node_attrib_df[,1:4], by.x="name", by.y="name",)
igraph_node_attrib_df_summarized_3 <- unique(igraph_node_attrib_df_summarized_2)

save(igraph_node_attrib_df_summarized_3, file = "igraph_node_attrib_df_summarized_3.RData")

# A basic scatterplot with color depending on Species
names(igraph_node_attrib_df_summarized_3)

ggplot(igraph_node_attrib_df_summarized_3, aes(x=as.factor(status), y=TL, color=body_size)) + 
  geom_point(size=6) +
  theme_ipsum()

