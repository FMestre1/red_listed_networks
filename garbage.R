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


################################################################################
#                                   FIGURES
################################################################################

#FMestre
#10-03-2023

#install.packages("treemap")
library(treemap)

View(red_listed_3)
length(fw_list)

red_listed_3[red_listed_3$group == "Amphibians_Reptiles",]

#Sent to Excel to create a pivot table...
write.csv(red_listed_3, "red_listed_3.csv", row.names=FALSE)
pivot_table1 <- read.csv("pivot_table1.csv", sep = ";")

pivot_table1[pivot_table1$group == "Amphibians_Reptiles",]$group <- "Amphibians and Reptiles"

#Plot
png(filename = "tree.png",width = 2000, height = 1600)

treemap(pivot_table1, index=c("group","status"), 
        fontsize.labels=c(25,20),  
        vSize="count", 
        type="index",
        border.col=c("black","black"),
        border.lwds=c(2,1),
        bg.labels=0,
        overlap.labels = 0,
        palette = "Set2",# Width of colors
        #title= "Species status per group"
        
)

dev.off()

################################################################################
################################################################################
#          Exploratory plotting of these data - Agregated status
################################################################################
#FMestre
#22-02-2023

library(ggplot2)

not_threat <- rbind(lc_table,
                    nt_table)

threatened <- rbind(vu_table,
                    en_table,
                    cr_table)

#in-degree #####################################################################
t_in<-data.frame(rep("threatened", nrow(threatened)), threatened$indegree)
names(t_in) <- c("iucn", "indegree")
head(t_in)
#
nt_in<-data.frame(rep("non-threatened", nrow(not_threat)), not_threat$indegree)
names(nt_in) <- c("iucn", "indegree")
head(nt_in)

INDEGREE_2 <- rbind(t_in,
                    nt_in
)

rm(t_in,
   nt_in
)

#save(INDEGREE_2, file = "INDEGREE_2.RData")

INDEGREE_2 <- INDEGREE_2[complete.cases(INDEGREE_2),]
str(INDEGREE_2)
INDEGREE_2$iucn <- as.factor(INDEGREE_2$iucn)

ggplot(INDEGREE_2, aes(x=fct_reorder(iucn,indegree, .desc=TRUE), y=indegree)) +
  ggtitle("In-degree for threatened and non-threatened species") + 
  xlab("IUCN Categories") + ylab("In-degree") +
  geom_boxplot()

#out-degree #####################################################################
t_out<-data.frame(rep("threatened", nrow(threatened)), threatened$outdegree)
names(t_out) <- c("iucn", "outdegree")
head(t_out)
#
nt_out<-data.frame(rep("non-threatened", nrow(not_threat)), not_threat$outdegree)
names(nt_out) <- c("iucn", "outdegree")
head(nt_out)

OUTDEGREE_2 <- rbind(t_out,
                     nt_out
)

rm(t_out,
   nt_out
)

#save(OUTDEGREE_2, file = "OUTDEGREE_2.RData")

OUTDEGREE_2 <- OUTDEGREE_2[complete.cases(OUTDEGREE_2),]
str(OUTDEGREE_2)
OUTDEGREE_2$iucn <- as.factor(OUTDEGREE_2$iucn)

ggplot(OUTDEGREE_2, aes(x=fct_reorder(iucn,outdegree, .desc=TRUE), y=outdegree)) +
  ggtitle("Out-degree for threatened and non-threatened species") + 
  xlab("IUCN Categories") + ylab("Out-degree") +
  geom_boxplot()

#centrality #####################################################################
t_centrality<-data.frame(rep("threatened", nrow(threatened)), threatened$centrality)
names(t_centrality) <- c("iucn", "centrality")
head(t_centrality)
#
nt_centrality<-data.frame(rep("non-threatened", nrow(not_threat)), not_threat$centrality)
names(nt_centrality) <- c("iucn", "centrality")
head(nt_centrality)

CENTRALITY_2 <- rbind(t_centrality,
                      nt_centrality
)

rm(t_centrality,
   nt_centrality
)

#save(CENTRALITY_2, file = "CENTRALITY_2.RData")

CENTRALITY_2 <- CENTRALITY_2[complete.cases(CENTRALITY_2),]
str(CENTRALITY_2)
CENTRALITY_2$iucn <- as.factor(CENTRALITY_2$iucn)

ggplot(CENTRALITY_2, aes(x=fct_reorder(iucn,centrality, .desc=TRUE), y=centrality)) +
  ggtitle("Centrality for threatened and non-threatened species") + 
  xlab("IUCN Categories") + ylab("Centrality") +
  geom_boxplot()

#centrality #####################################################################
t_ivi<-data.frame(rep("threatened", nrow(threatened)), threatened$ivi)
names(t_ivi) <- c("iucn", "ivi")
head(t_ivi)
#
nt_ivi<-data.frame(rep("non-threatened", nrow(not_threat)), not_threat$ivi)
names(nt_ivi) <- c("iucn", "ivi")
head(nt_ivi)

IVI_2 <- rbind(t_ivi,
               nt_ivi
)

rm(t_ivi,
   nt_ivi
)

#save(IVI_2, file = "IVI_2.RData")

IVI_2 <- IVI_2[complete.cases(IVI_2),]
str(IVI_2)
IVI_2$iucn <- as.factor(IVI_2$iucn)

ggplot(IVI_2, aes(x=fct_reorder(iucn,ivi, .desc=TRUE), y=ivi)) +
  ggtitle("IVI index for threatened and non-threatened species") + 
  xlab("IUCN Categories") + ylab("IVI") +
  geom_boxplot()

#closeness #####################################################################
t_closeness<-data.frame(rep("threatened", nrow(threatened)), threatened$closeness)
names(t_closeness) <- c("iucn", "closeness")
head(t_closeness)
#
nt_closeness<-data.frame(rep("non-threatened", nrow(not_threat)), not_threat$closeness)
names(nt_closeness) <- c("iucn", "closeness")
head(nt_closeness)

CLOSENESS_2 <- rbind(t_closeness,
                     nt_closeness
)

rm(t_closeness,
   nt_closeness
)

#save(CLOSENESS_2, file = "CLOSENESS_2.RData")

CLOSENESS_2 <- CLOSENESS_2[complete.cases(CLOSENESS_2),]
str(CLOSENESS_2)
CLOSENESS_2$iucn <- as.factor(CLOSENESS_2$iucn)

ggplot(CLOSENESS_2, aes(x=fct_reorder(iucn,closeness, .desc=TRUE), y=closeness)) +
  ggtitle("Closeness index for threatened and non-threatened species") + 
  xlab("IUCN Categories") + ylab("IVI") +
  geom_boxplot()


################################################################################
################################################################################
#                               RELATING BOTH
################################################################################
#FMestre
#08-02-2023

#Load packages
library(taxize)

#Red list data - get the names of the species
red_listed_3
red_list_species <- red_listed_3$full_name
#red_listed_3$full_name <- stringr::str_replace(red_listed_3$full_name, "_", " ")

#Node metrics across Europe
#fw_list
species_names2 #the species on the networks across Europe
species_names3 <- stringr::str_replace(species_names2, "_", " ")
species_names4 <- taxize::tax_name(sci=species_names3, get=c("kingdom", "class","order","family","genus"), db="itis") #changed to "itis"
#save(species_names4, file = "species_names4.RData")

#Merge both datasets
merged_tables <- merge(x = red_listed_3, y = species_names4, by.x = "full_name", by.y = "query", all=TRUE)
#View(merged_tables)
#merged_tables[!is.na(merged_tables$europeanRegionalRedListCategory) & !is.na(merged_tables$class),]

#Save as csv file
write.csv(merged_tables, "merged_tables_09_FEV_2023.csv", row.names=FALSE)
#read.csv("merged_tables_09_FEV_2023.csv")

#Only those in the red list data were not merged? Strange?
#table(species_names4$query %in% red_listed_3$full_name) #species on FW that are in Red List
#table(red_listed_3$full_name %in% species_names4$query) #species in red list that are in the FW  
#Ok, we have 190 non-matches!

#Check with synonyms (taxize)
not_matched <- species_names4[!(species_names4$query %in% red_listed_3$full_name),]

list_of_syns <- list()

for(j in 1:nrow(not_matched)){
  sp1 <- not_matched$query[j]
  try(syn_sp1 <- taxize::nbn_synonyms(sp1), silent = TRUE)
  if(exists("syn_sp1")) {
    syn_sp1 <- syn_sp1$nameString
    list_of_syns[[j]] <- syn_sp1
    rm(syn_sp1)
  } else list_of_syns[[j]] <- NA
  
  message(j)
  
}

names(list_of_syns) <- not_matched$query

#Do these synonyms improve things?

additional_matches <- c()

for(i in 1:length(list_of_syns)){
  
  syns2 <- list_of_syns[[i]]
  
  additional_matches[i] <- any(syns2 %in% red_listed_3$full_name)
  
}

table(additional_matches)
#only two! So... the rest are not in the list!

which(additional_matches == TRUE)
#(the 1 and 5 are solved by syns)

list_of_syns[[1]] %in% red_listed_3$full_name
list_of_syns[[5]] %in% red_listed_3$full_name

red_listed_3[red_listed_3$full_name == "Catharacta skua",]
red_listed_3[red_listed_3$full_name == "Eudromias morinellus",]

#Add this to the merged table
View(merged_tables)

merged_tables[merged_tables$full_name == "Catharacta skua",]$kingdom <- species_names4[species_names4$query == names(list_of_syns)[1],]$kingdom
merged_tables[merged_tables$full_name == "Catharacta skua",]$class <- species_names4[species_names4$query == names(list_of_syns)[1],]$class
merged_tables[merged_tables$full_name == "Catharacta skua",]$order <- species_names4[species_names4$query == names(list_of_syns)[1],]$order
merged_tables[merged_tables$full_name == "Catharacta skua",]$family <- species_names4[species_names4$query == names(list_of_syns)[1],]$family
merged_tables[merged_tables$full_name == "Catharacta skua",]$genus <- species_names4[species_names4$query == names(list_of_syns)[1],]$genus

merged_tables[merged_tables$full_name == "Eudromias morinellus",]$kingdom <- species_names4[species_names4$query == names(list_of_syns)[5],]$kingdom
merged_tables[merged_tables$full_name == "Eudromias morinellus",]$class <- species_names4[species_names4$query == names(list_of_syns)[5],]$class
merged_tables[merged_tables$full_name == "Eudromias morinellus",]$order <- species_names4[species_names4$query == names(list_of_syns)[5],]$order
merged_tables[merged_tables$full_name == "Eudromias morinellus",]$family <- species_names4[species_names4$query == names(list_of_syns)[5],]$family
merged_tables[merged_tables$full_name == "Eudromias morinellus",]$genus <- species_names4[species_names4$query == names(list_of_syns)[5],]$genus

#Finally.... keep only those that are in the FW
fw_species_with_red_list_status <- merged_tables[!is.na(merged_tables$genus),]
fw_species_with_red_list_status <- fw_species_with_red_list_status[,-c(2,5)]

fw_species_with_red_list_status$europeanRegionalRedListCategory[is.na(fw_species_with_red_list_status$europeanRegionalRedListCategory)] <- "not_listed"
fw_species_with_red_list_status$endemic_to_europe[is.na(fw_species_with_red_list_status$endemic_to_europe)] <- "not_listed"
#View(fw_species_with_red_list_status)

#Save
#save(fw_species_with_red_list_status, file = "fw_species_with_red_list_status_15_FEV_2023.RData")

unique(fw_species_with_red_list_status$europeanRegionalRedListCategory)

grouped_status <- c()

for(i in 1:nrow(fw_species_with_red_list_status)){
  
  st1 <- fw_species_with_red_list_status$europeanRegionalRedListCategory[i]
  if(st1 == "not_listed") grouped_status[i] <- "not_listed"
  if(st1 == "DD") grouped_status[i] <- "data_deficient"
  if(st1 == "LC") grouped_status[i] <- "not_threatened"
  if(st1 == "NT") grouped_status[i] <- "not_threatened"
  if(st1 == "VU") grouped_status[i] <- "threatened"
  if(st1 == "EN") grouped_status[i] <- "threatened"
  if(st1 == "CR") grouped_status[i] <- "threatened"
  if(st1 == "NE") grouped_status[i] <- "not_evaluated"
  if(st1 == "RE") grouped_status[i] <- "regionally_extinct"
  
}

fw_species_with_red_list_combined_status <- data.frame(fw_species_with_red_list_status, grouped_status)
#View(fw_species_with_red_list_combined_status)

#Save
#save(fw_species_with_red_list_combined_status, file = "fw_species_with_red_list_combined_status_26_JUN_2023.RData")

##########################################

#Get the metrics per status in all trophic structures
#names(fw_list[[1]])[-1]
dd_table <- as.data.frame(matrix(ncol = 7))
lc_table <- as.data.frame(matrix(ncol = 7))
nt_table <- as.data.frame(matrix(ncol = 7))
vu_table <- as.data.frame(matrix(ncol = 7))
en_table <- as.data.frame(matrix(ncol = 7))
cr_table <- as.data.frame(matrix(ncol = 7))
ne_table <- as.data.frame(matrix(ncol = 7))
re_table <- as.data.frame(matrix(ncol = 7))
names(dd_table) <- names(lc_table) <- names(nt_table) <- names(vu_table) <- names(en_table) <- names(cr_table) <- names(ne_table) <- names(re_table) <- names(fw_list[[1]])[-1]  

for(i in 1:length(fw_list)){
  
  net1 <- fw_list[[i]]
  species_net1 <- fw_list[[i]]$SP_NAME
  species_net1 <- stringr::str_replace(species_net1, "_", " ")
  
  if(any(species_net1 %in% fw_species_with_red_list_combined_status$full_name)){
    
    listed_species <- species_net1[species_net1 %in% fw_species_with_red_list_combined_status$full_name]
    n_species <- length(listed_species)
    listed_species_status <- fw_species_with_red_list_combined_status[fw_species_with_red_list_combined_status$full_name %in% listed_species,]
    
    if(any(listed_species_status$europeanRegionalRedListCategory == "DD")){
      species_dd <- listed_species_status[listed_species_status$europeanRegionalRedListCategory == "DD",]
      species_dd <- species_dd$full_name
      species_dd <- stringr::str_replace(species_dd, " ", "_")
      dd_table <- rbind(dd_table, net1[net1$SP_NAME %in% species_dd,][,-1])      
    }
    if(any(listed_species_status$europeanRegionalRedListCategory == "LC")){
      species_lc <- listed_species_status[listed_species_status$europeanRegionalRedListCategory == "LC",]
      species_lc <- species_lc$full_name
      species_lc <- stringr::str_replace(species_lc, " ", "_")
      lc_table <- rbind(lc_table, net1[net1$SP_NAME %in% species_lc,][,-1])
    }
    if(any(listed_species_status$europeanRegionalRedListCategory == "NT")){
      species_nt <- listed_species_status[listed_species_status$europeanRegionalRedListCategory == "NT",]
      species_nt <- species_nt$full_name
      species_nt <- stringr::str_replace(species_nt, " ", "_")
      nt_table <- rbind(nt_table, net1[net1$SP_NAME %in% species_nt,][,-1])
    }
    if(any(listed_species_status$europeanRegionalRedListCategory == "VU")){
      species_vu <- listed_species_status[listed_species_status$europeanRegionalRedListCategory == "VU",]
      species_vu <- species_vu$full_name
      species_vu <- stringr::str_replace(species_vu, " ", "_")
      vu_table <- rbind(vu_table, net1[net1$SP_NAME %in% species_vu,][,-1])
    }
    if(any(listed_species_status$europeanRegionalRedListCategory == "EN")){
      species_en <- listed_species_status[listed_species_status$europeanRegionalRedListCategory == "EN",]
      species_en <- species_en$full_name
      species_en <- stringr::str_replace(species_en, " ", "_")
      en_table <- rbind(en_table, net1[net1$SP_NAME %in% species_en,][,-1])
    }
    if(any(listed_species_status$europeanRegionalRedListCategory == "CR")){
      species_cr <- listed_species_status[listed_species_status$europeanRegionalRedListCategory == "CR",]
      species_cr <- species_cr$full_name
      species_cr <- stringr::str_replace(species_cr, " ", "_")
      cr_table <- rbind(cr_table, net1[net1$SP_NAME %in% species_cr,][,-1])
    }
    if(any(listed_species_status$europeanRegionalRedListCategory == "NE")){
      species_ne <- listed_species_status[listed_species_status$europeanRegionalRedListCategory == "NE",]
      species_ne <- species_ne$full_name
      species_ne <- stringr::str_replace(species_ne, " ", "_")
      ne_table <- rbind(ne_table, net1[net1$SP_NAME %in% species_ne,][,-1])
    }
    if(any(listed_species_status$europeanRegionalRedListCategory == "RE")){
      species_re <- listed_species_status[listed_species_status$europeanRegionalRedListCategory == "RE",]
      species_re <- species_re$full_name
      species_re <- stringr::str_replace(species_re, " ", "_")
      re_table <- rbind(re_table, net1[net1$SP_NAME %in% species_re,][,-1])
    }
  }
  
  message(i)
  
  
}

#Saving of partial tables - Shouldn't these have all the other number, not just 19?
save(dd_table, file = "dd_table_CORRECTEC_19.RData")
save(lc_table, file = "lc_table_CORRECTEC_19.RData")
save(nt_table, file = "nt_table_CORRECTEC_19.RData")
save(vu_table, file = "vu_table_CORRECTEC_19.RData")
save(en_table, file = "en_table_CORRECTEC_19.RData")
save(cr_table, file = "cr_table_CORRECTEC_19.RData")
save(ne_table, file = "ne_table_CORRECTEC_19.RData")
save(re_table, file = "re_table_CORRECTEC_19.RData")


################################################################################
#                         Plotting FW for the figure
################################################################################

#Load datasets
#Load from cluster
load("from_cluster/network_list_cheddar_06jun23.RData")
#load("from_cluster/eur_comm_collection_06jun23.RData")
load("from_cluster/network_list_igraph_2_all_13JUN2023.RData")

#Package
library(cheddar)
library(igraph)

# Get a good example to the fw figure ##########################################

#using igraph list
for(i in 1:length(network_list_igraph_2)){
  
  if(all(c("Lynx pardinus", "Oryctolagus cuniculus", "Vulpes vulpes") %in% igraph::vertex.attributes(network_list_igraph_2[[i]])[1]$name)) print(i)
  #if(any(c("Lynx pardinus", "Oryctolagus cuniculus") %in% igraph::vertex.attributes(network_list_igraph_2[[i]])[1])) print (i)
  #if(any(c("Ursus maritimus") %in% igraph::vertex.attributes(network_list_igraph_2[[i]])[1]$name)) print(i)
  
}

#using cheddar list
for(i in 1:length(network_list_cheddar)){
  
  if(all(c("Lynx pardinus", "Oryctolagus cuniculus", "Vulpes vulpes") %in%  network_list_cheddar[[i]]$nodes$node)) print(i)
  #if(any(c("Lynx pardinus", "Oryctolagus cuniculus") %in% igraph::vertex.attributes(network_list_igraph_2[[i]])[1])) print (i)
  #if(any(c("Ursus maritimus") %in% igraph::vertex.attributes(network_list_igraph_2[[i]])[1]$name)) print(i)
  
}

#cheddar::TrophicLinkPropertyNames (network_list_cheddar[[116183]])
#cheddar::NodePropertyNames (network_list_cheddar[[116183]])

#Which is the most vulnerable prey?
sort(cheddar::TrophicVulnerability(network_list_cheddar[[116183]]))
"Pelophylax perezi" # this is it, but also use...
"Oryctolagus cuniculus"

#Which is the most generalist predator
sort(cheddar::TrophicGenerality(network_list_cheddar[[116183]]))
"Vulpes vulpes" # this is it, but also use...
"Lynx pardinus"

#Lets use the grids SQ189!

#Plot links going and coming from Vulpes vulpes
links <- cbind(TLPS(network_list_cheddar[[116183]]), colour="#c7c7c788")
#
links <- cbind(TLPS(network_list_cheddar[[116183]]), colour="#ffffff00")
links$colour["Vulpes vulpes" == links$resource] <- "red"
links$colour["Vulpes vulpes" == links$consumer] <- "darkgreen"

#cheddar::plot.Community(network_list_cheddar[[116183]], node.labels="node", show.nodes.as="both", link.col=links$colour)
cheddar::plot.Community(network_list_cheddar[[116183]], link.col=links$colour)
#

#Plot links going and coming from Lynx pardinus
links0 <- cbind(TLPS(network_list_cheddar[[116183]]), colour="#c7c7c788")
#
links0 <- cbind(TLPS(network_list_cheddar[[116183]]), colour="#ffffff00")
links0$colour["Lynx pardinus" == links0$resource] <- "red"
links0$colour["Lynx pardinus" == links0$consumer] <- "darkgreen"

#cheddar::plot.Community(network_list_cheddar[[116183]], node.labels="node", show.nodes.as="both", link.col=links$colour)
cheddar::plot.Community(network_list_cheddar[[116183]], link.col=links0$colour)
#

#Plot links coming from Oryctolagus cuniculus
# transparent #ffffff00
links1 <- cbind(TLPS(network_list_cheddar[[116183]]), colour="#c7c7c788")
#
links1 <- cbind(TLPS(network_list_cheddar[[116183]]), colour="#ffffff00")
links1$colour["Oryctolagus cuniculus" == links1$resource] <- "red"
links1$colour["Oryctolagus cuniculus" == links1$consumer] <- "darkgreen"

#cheddar::plot.Community(network_list_cheddar[[116183]], node.labels="node", show.nodes.as="both", link.col=links$colour)
cheddar::plot.Community(network_list_cheddar[[116183]], link.col=links1$colour)
#

############################################################

nodes1 <- cbind(NPS(network_list_cheddar[[116183]]), colour="darkgreen")
#
#nodes1 <- cbind(NPS(network_list_cheddar[[116183]]), colour="#ffffff00")

nodes1$colour["threatened" == nodes1$agreg_ts] <- "red"
nodes1$colour["not_threatened" == nodes1$agreg_ts] <- "darkgreen"

#cheddar::plot.Community(network_list_cheddar[[116183]], node.labels="node", show.nodes.as="both", link.col=links$colour)

pch_vector <- rep(19, length=length(nodes1$colour))
pch_cex <- rep(1, length=length(nodes1$colour))

for(i in 1:length(pch_cex)) if(nodes1$colour[i] == "red") pch_cex[i] <- 2

#The links are setup as transparent: link.col="#ffffff00"
cheddar::plot.Community(network_list_cheddar[[116183]], link.col="#ffffff00", col = nodes1$colour, pch = 19, cex = pch_cex)
#

#DELETE
rm(network_list_igraph_2)
rm(network_list_cheddar)



#######################################################################################
#    Evaluate statistical differences between threatened and non-threatened species
#######################################################################################

proportion <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/proportion_spatial.shp")
richness <- terra::vect("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\shapefiles\\richness_2.shp")

#T-test and Cohen's d

#t test interpretation
#p < 0.05 significance difference

#Coehen's d interpretation
#Small effect size: Cohen's d around 0.2 indicates a small effect size. 
#It suggests that there is a small difference between the groups, which may 
#have limited practical significance.

#Medium effect size: Cohen's d around 0.5 indicates a medium effect size. 
#It suggests a moderate difference between the groups, which is typically 
#considered to have a meaningful practical impact.

#Large effect size: Cohen's d around 0.8 or above indicates a large effect size. 
#It suggests a substantial difference between the groups, which is likely to 
#have significant practical implications.

nt_indegree <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/indegree_nt_spatial.shp")
t_indegree <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/indegree_t_spatial.shp")
nt_indegree <- as.data.frame(nt_indegree)
t_indegree <- as.data.frame(t_indegree)
nt_indegree <- nt_indegree[,c(1,4)]
t_indegree <- t_indegree[,c(1,4)]
names(nt_indegree)[2] <- "NT_indegree"
names(t_indegree)[2] <- "T_indegree"
indegree_compare <- merge(nt_indegree, t_indegree)
indegree_compare <- indegree_compare[complete.cases(indegree_compare),]
#View(indegree_compare)
indegree_ttest <- t.test(indegree_compare[,2], indegree_compare[,3], paired = TRUE)
print(indegree_ttest)
indegree_cohens_d <- effsize::cohen.d(indegree_compare[,2], indegree_compare[,3])
print(indegree_cohens_d)
plot(indegree_compare[,2], indegree_compare[,3])
cor(indegree_compare[,2], indegree_compare[,3], method = "pearson")
#indegree_chisquare <- chisq.test(indegree_compare[,2:3], correct = FALSE)
#indegree_fisher_test <- fisher.test(indegree_compare[,2:3])

#frequency plots
names(indegree_compare)
par(mfrow=c(2, 2))
h_in_nt <- hist(indegree_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "In-degree", breaks = seq(0,65, by=2), ylim = c(0,40000), plot = FALSE)
h_in_t <- hist(indegree_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "In-degree",  breaks = seq(0,65, by=2), ylim = c(0,40000), plot = FALSE)
# Convert the counts to percentages
h_in_nt$counts <- h_in_nt$counts / sum(h_in_nt$counts) * 100
h_in_t$counts <- h_in_t$counts / sum(h_in_t$counts) * 100
#plot(h_in_nt, main = "Non-threatened", col = "darkgreen", xlab = "In-degree")
#plot(h_in_t, main = "Threatened", col = "darkred", xlab = "In-degree")
# Create the line plot
plot(h_in_nt$mids, h_in_nt$counts, type = "n", xlab = "In-degree", ylab = "Frequency (%)", main = "In-degree", ylim = c(0,40))
lines(h_in_nt$mids, h_in_nt$counts, lwd = 3, col = "darkgreen")
lines(h_in_t$mids, h_in_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)
#
nt_outdegree <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/outdegree_nt_spatial.shp")
t_outdegree <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/outdegree_t_spatial.shp")
nt_outdegree <- as.data.frame(nt_outdegree)
t_outdegree <- as.data.frame(t_outdegree)
nt_outdegree <- nt_outdegree[,c(1,4)]
t_outdegree <- t_outdegree[,c(1,4)]
names(nt_outdegree)[2] <- "NT_outdegree"
names(t_outdegree)[2] <- "T_outdegree"
outdegree_compare <- merge(nt_outdegree, t_outdegree)
outdegree_compare <- outdegree_compare[complete.cases(outdegree_compare),]
#View(outdegree_compare)
outdegree_ttest <- t.test(outdegree_compare[,2], outdegree_compare[,3], paired = TRUE)
print(outdegree_ttest)
outdegree_cohens_d <- effsize::cohen.d(outdegree_compare[,2], outdegree_compare[,3])
print(outdegree_cohens_d)
plot(outdegree_compare[,2], outdegree_compare[,3])
cor(outdegree_compare[,2], outdegree_compare[,3], method = "pearson")

#frequency plots
names(outdegree_compare)
par(mfrow=c(1, 2))
h_out_nt <- hist(outdegree_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "Out-degree", breaks = seq(0,65, by=2), ylim = c(0,40000))
h_out_t <- hist(outdegree_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "Out-degree",  breaks = seq(0,65, by=2), ylim = c(0,40000))
# Convert the counts to percentages
h_out_nt$counts <- h_out_nt$counts / sum(h_out_nt$counts) * 100
h_out_t$counts <- h_out_t$counts / sum(h_out_t$counts) * 100
# Create the line plot
plot(h_out_nt$mids, h_out_nt$counts, type = "n", xlab = "Out-degree", ylab = "Frequency (%)", main = "Out-degree", ylim = c(0,40))
lines(h_out_nt$mids, h_out_nt$counts, lwd = 3, col = "darkgreen")
lines(h_out_t$mids, h_out_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)
#
nt_t_level <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/tl_nt_spatial.shp")
t_t_level <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/tl_t_spatial.shp")
nt_t_level <- as.data.frame(nt_t_level)
t_t_level <- as.data.frame(t_t_level)
nt_t_level <- nt_t_level[,c(1,4)]
t_t_level <- t_t_level[,c(1,4)]
names(nt_t_level)[2] <- "NT_trph_level"
names(t_t_level)[2] <- "T_trph_level"
trophic_level_compare <- merge(nt_t_level, t_t_level)
trophic_level_compare <- trophic_level_compare[complete.cases(trophic_level_compare),]
#View(trophic_level_compare)
trophic_level_ttest <- t.test(trophic_level_compare[,2], trophic_level_compare[,3], paired = TRUE)
print(trophic_level_ttest)
trophic_level_cohens_d <- effsize::cohen.d(trophic_level_compare[,2], trophic_level_compare[,3])
print(trophic_level_cohens_d)
plot(trophic_level_compare[,2], trophic_level_compare[,3])
cor(trophic_level_compare[,2], trophic_level_compare[,3], method = "pearson")

#frequency plots
names(trophic_level_compare)
par(mfrow=c(1, 2))
h_tl_nt <- hist(trophic_level_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "Trophic level", breaks = seq(0,4, by=0.25), ylim = c(0,100000))
h_tl_t <- hist(trophic_level_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "Trophic level",  breaks = seq(0,4, by=0.25), ylim = c(0,100000))
# Convert the counts to percentages
h_tl_nt$counts <- h_tl_nt$counts / sum(h_tl_nt$counts) * 100
h_tl_t$counts <- h_tl_t$counts / sum(h_tl_t$counts) * 100
# Create the line plot
plot(h_tl_nt$mids, h_tl_nt$counts, type = "n", xlab = "Trophic level", ylab = "Frequency (%)", main = "Trophic level", ylim = c(0,100))
lines(h_tl_nt$mids, h_tl_nt$counts, lwd = 3, col = "darkgreen")
lines(h_tl_t$mids, h_tl_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)
#
nt_closeness <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/closeness_nt_spatial.shp")
t_closeness <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/closeness_t_spatial.shp")
nt_closeness <- as.data.frame(nt_closeness)
t_closeness <- as.data.frame(t_closeness)
nt_closeness <- nt_closeness[,c(1,4)]
t_closeness <- t_closeness[,c(1,4)]
names(nt_closeness)[2] <- "NT_closeness"
names(t_closeness)[2] <- "T_closeness"
closeness_compare <- merge(nt_closeness, t_closeness)
closeness_compare <- closeness_compare[complete.cases(closeness_compare),]
#View(closeness_compare)
closeness_ttest <- t.test(closeness_compare[,2], closeness_compare[,4], paired = TRUE)
print(closeness_ttest)
outdegree_cohens_d <- effsize::cohen.d(closeness_compare[,2], closeness_compare[,3])
print(outdegree_cohens_d)
plot(closeness_compare[,2], closeness_compare[,3])
cor(closeness_compare[,2], closeness_compare[,3], method = "pearson")

#frequency plots
names(closeness_compare)
par(mfrow=c(1, 2))
h_clos_nt <- hist(closeness_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "Closeness Centrality", breaks = seq(0,460, by=10), ylim = c(0,150000))
h_clos_t <- hist(closeness_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "Closeness Centrality",  breaks = seq(0,460, by=10), ylim = c(0,150000))
# Convert the counts to percentages
h_clos_nt$counts <- h_clos_nt$counts / sum(h_clos_nt$counts) * 100
h_clos_t$counts <- h_clos_t$counts / sum(h_clos_t$counts) * 100
# Create the line plot
plot(h_clos_nt$mids, h_clos_nt$counts, type = "n", xlab = "Closeness Centrality", ylab = "Frequency (%)", main = "Closeness Centrality", ylim = c(0,100))
lines(h_clos_nt$mids, h_clos_nt$counts, lwd = 3, col = "darkgreen")
lines(h_tl_t$mids, h_tl_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)
#
nt_centrality <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/centrality_nt_spatial.shp")
t_centrality <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/centrality_t_spatial.shp")
nt_centrality <- as.data.frame(nt_centrality)
t_centrality <- as.data.frame(t_centrality)
nt_centrality <- nt_centrality[,c(1,4)]
t_centrality <- t_centrality[,c(1,4)]
names(nt_centrality)[2] <- "NT_centrality"
names(t_centrality)[2] <- "T_centrality"
centrality_compare <- merge(nt_centrality, t_centrality)
centrality_compare <- centrality_compare[complete.cases(centrality_compare),]
#View(centrality_compare)
centrality_ttest <- t.test(centrality_compare[,2], centrality_compare[,3], paired = TRUE)
print(centrality_ttest)
centrality_cohens_d <- effsize::cohen.d(centrality_compare[,2], centrality_compare[,3])
print(centrality_cohens_d)
plot(centrality_compare[,2], centrality_compare[,3])
cor(centrality_compare[,2], centrality_compare[,3], method = "pearson")

#frequency plots
names(centrality_compare)
par(mfrow=c(1, 2))
h_bet_nt <- hist(centrality_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "Betwenness Centrality", breaks = seq(0,460, by=10), ylim = c(0,150000))
h_bet_t <- hist(centrality_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "Betwenness Centrality",  breaks = seq(0,460, by=10), ylim = c(0,150000))
# Convert the counts to percentages
h_bet_nt$counts <- h_bet_nt$counts / sum(h_bet_nt$counts) * 100
h_bet_t$counts <- h_bet_t$counts / sum(h_bet_t$counts) * 100
# Create the line plot
plot(h_bet_nt$mids, h_bet_nt$counts, type = "n", xlab = "Betwenness Centrality", ylab = "Frequency (%)", main = "Betwenness Centrality", ylim = c(0,65))
lines(h_bet_nt$mids, h_bet_nt$counts, lwd = 3, col = "darkgreen")
lines(h_tl_t$mids, h_tl_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)
#
nt_ivi <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/ivi_nt_spatial_second_version.shp")
t_ivi <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/ivi_t_spatial_second_version.shp")
nt_ivi <- as.data.frame(nt_ivi)
t_ivi <- as.data.frame(t_ivi)
nt_ivi <- nt_ivi[,c(1,4)]
t_ivi <- t_ivi[,c(1,4)]
names(nt_ivi)[2] <- "NT_ivi"
names(t_ivi)[2] <- "T_ivi"
ivi_compare <- merge(nt_ivi, t_ivi)
ivi_compare <- ivi_compare[complete.cases(ivi_compare),]
#View(centrality_compare)
ivi_ttest <- t.test(ivi_compare[,2], ivi_compare[,3], paired = TRUE)
print(ivi_ttest)
ivi_cohens_d <- effsize::cohen.d(ivi_compare[,2], ivi_compare[,3])
print(ivi_cohens_d)
plot(ivi_compare[,2], ivi_compare[,3])
cor(ivi_compare[,2], ivi_compare[,3], method = "pearson")
#
#frequency plots
names(ivi_compare)
par(mfrow=c(1, 2))
h_ivi_nt <- hist(ivi_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "IVI", breaks = seq(0,100, by=10), ylim = c(0,150000))
h_ivi_t <- hist(ivi_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "IVI",  breaks = seq(0,100, by=10), ylim = c(0,150000))
# Convert the counts to percentages
h_ivi_nt$counts <- h_ivi_nt$counts / sum(h_ivi_nt$counts) * 100
h_ivi_t$counts <- h_ivi_t$counts / sum(h_ivi_t$counts) * 100
# Create the line plot
plot(h_ivi_nt$mids, h_ivi_nt$counts, type = "n", xlab = "IVI", ylab = "Frequency (%)", main = "IVI", ylim = c(0,100))
lines(h_ivi_nt$mids, h_ivi_nt$counts, lwd = 3, col = "darkgreen")
lines(h_ivi_t$mids, h_ivi_t$counts, lwd = 3, col = "darkred")
# Add legend
legend("topright", legend = c("Non-threatened", "Threatened"), col = c("darkgreen", "darkred"), lwd = 3)
#



#Plot with rastervis
##
my.col.regions <- rev(terrain.colors(100))
######################
#max(nt_indegree) #max value: 25.28488 
#max(t_indegree) #max value: 61.000 
ind_min_max <- seq(0, 61, length.out = 100)
ind1 <- rasterVis::levelplot(nt_indegree, col.regions=my.col.regions, at=ind_min_max, main = "Not-threatened")
ind2 <- rasterVis::levelplot(t_indegree, col.regions=my.col.regions, at=ind_min_max, main = "Threatened")
ind_title <- textGrob("Average In-degree", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(ind1, ind2, ncol=2, top=ind_title)
######################
#max(nt_outdegree) #max value: 22.2991447
#max(t_outdegree) #max value: 20.4444447
out_min_max <- seq(0, 23, length.out = 100)
out1 <- rasterVis::levelplot(nt_outdegree, col.regions=my.col.regions, at=out_min_max, main = "Not-threatened")
out2 <- rasterVis::levelplot(t_outdegree, col.regions=my.col.regions, at=out_min_max, main = "Threatened")
out_title <- textGrob("Average Out-degree", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(out1, out2, ncol=2, top=out_title)
######################
#max(nt_t_level) #max value: 2
#max(t_t_level) #max value: 3 
tl_min_max <- seq(0, 3, length.out = 100)
tl1 <- rasterVis::levelplot(nt_t_level, col.regions=my.col.regions, at=tl_min_max, main = "Not-threatened")
tl2 <- rasterVis::levelplot(t_t_level, col.regions=my.col.regions, at=tl_min_max, main = "Threatened")
tl_title <- textGrob("Average Trophic Level", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(tl1, tl2, ncol=2, top=tl_title)
######################
#max(nt_closeness) #max value:
#max(t_closeness) #max value:
closeness_min_max <- seq(0, 1, length.out = 100)
cl1 <- rasterVis::levelplot(nt_closeness, col.regions=my.col.regions, at=closeness_min_max, main = "Not-threatened")
cl2 <- rasterVis::levelplot(t_closeness, col.regions=my.col.regions, at=closeness_min_max, main = "Threatened")
cl_title <- textGrob("Average Closeness Centrality", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(cl1, cl2, ncol=2, top=cl_title)
######################
#max(nt_centrality) #max value: 115.24037170
#max(t_centrality) #max value: 459.1808
centrality_min_max <- seq(0, 460, length.out = 100)
bt1 <- rasterVis::levelplot(nt_centrality, col.regions=my.col.regions, at=centrality_min_max, main = "Not-threatened")
bt2 <- rasterVis::levelplot(t_centrality, col.regions=my.col.regions, at=centrality_min_max, main = "Threatened")
bt_title <- textGrob("Average Betweenness Centrality", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(bt1, bt2, ncol=2, top=bt_title)
######################
#max(nt_ivi) #max value: 82.45447
#max(t_ivi) #max value: 100
ivi_min_max <- seq(0, 100, length.out = 100)
ivi1 <- rasterVis::levelplot(nt_ivi, col.regions=my.col.regions, at=ivi_min_max, main = "Not-threatened")
ivi2 <- rasterVis::levelplot(t_ivi, col.regions=my.col.regions, at=ivi_min_max, main = "Threatened")
ivi_title <- textGrob("Average IVI", gp = gpar(fontsize = 25, fontface = "bold"))
grid.arrange(ivi1, ivi2, ncol=2, top=ivi_title)
######################
rasterVis::levelplot(proportion, par.settings = rasterTheme(viridis_pal()(255)), main = "Proportion of threatened species")
rasterVis::levelplot(richness_raster2, par.settings = rasterTheme(viridis_pal()(255)), main = "Number of nodes in each network")


################################################################################
#                                  Plot maps
################################################################################
#FMestre
#19-07-2023

library(raster)
library(terra)
library(rasterVis)
library(cheddar)
library(effsize)
library(gridExtra)
library(viridis)
library(grid)

#Load all rasters
nt_indegree <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_indegree.tif")
t_indegree <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_indegree.tif")
#
nt_outdegree <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_outdegree.tif")
t_outdegree <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_outdegree.tif")
#
nt_t_level <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_trophic_level.tif")
t_t_level <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_trophic_level.tif")
#
nt_closeness <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_closeness.tif")
t_closeness <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_closeness.tif")
#
nt_centrality <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_centrality.tif")
t_centrality <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_centrality.tif")
#
nt_ivi <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_ivi.tif")
t_ivi <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_ivi.tif")
#
proportion <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\proportion.tif")
richness_raster <- terra::rast("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\richness.tif")
richness_raster2 <- terra::resample(richness_raster, proportion)
#

richness <- terra::vect("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\shapefiles\\richness_2.shp")
proportion_vector <- terra::vect("C:/Users/FMest/Documents/github/red_listed_networks/proportion_spatial.shp")
#
richness_proportion <- data.frame(merge(richness, proportion_vector))
#head(richness_proportion)
richness_proportion <- richness_proportion[,-c(2,3,4)]


####################################################################################

indegree_ttest <- t.test(indegree_compare[,2], indegree_compare[,3], paired = TRUE)
#print(indegree_ttest)
indegree_cohens_d <- effsize::cohen.d(indegree_compare[,2], indegree_compare[,3])
#print(indegree_cohens_d)
#plot(indegree_compare[,2], indegree_compare[,3])
cor(indegree_compare[,2], indegree_compare[,3], method = "pearson")

############

outdegree_ttest <- t.test(outdegree_compare[,2], outdegree_compare[,3], paired = TRUE)
#print(outdegree_ttest)
outdegree_cohens_d <- effsize::cohen.d(outdegree_compare[,2], outdegree_compare[,3])
#print(outdegree_cohens_d)
#plot(outdegree_compare[,2], outdegree_compare[,3])
cor(outdegree_compare[,2], outdegree_compare[,3], method = "pearson")

############

#FMestre
#10-07-2023

################################################################################
#                                  Plot maps
################################################################################

################################################################################
#   Where are areas significantly higher or lower in one metrics in T and NT?
################################################################################

names(indegree_compare)
names(outdegree_compare)
names(trophic_level_difference)
names(closeness_difference)
names(centrality_difference)
names(ivi_difference)
#
indegree_compare_2 <- data.frame(nt_indegree, t_indegree)
outdegree_compare_2 <- data.frame(nt_outdegree, t_outdegree)
trophic_level_compare_2 <- data.frame(nt_t_level, t_t_level)
closeness_compare_2 <- data.frame(nt_closeness, t_closeness)
centrality_compare_2 <- data.frame(nt_centrality, t_centrality)
ivi_compare_2 <- data.frame(nt_ivi, t_ivi)
#
indegree_difference <- indegree_compare_2[,4] - indegree_compare_2[,2]
outdegree_difference <- outdegree_compare_2[,4] - outdegree_compare_2[,2]
trophic_level_difference <- trophic_level_compare_2[,4] - trophic_level_compare_2[,2]
closeness_difference <- closeness_compare_2[,4] - closeness_compare_2[,2]
centrality_difference <- centrality_compare_2[,4] - centrality_compare_2[,2]
ivi_difference <- ivi_compare_2[,4] - ivi_compare_2[,2]

# Perform a t-test to assess the statistical significance
t_test_indeg <- t.test(indegree_compare_2[,4], indegree_compare_2[,2], paired = TRUE, alternative = "two.sided")
t_test_outdeg <- t.test(outdegree_compare_2[,4], outdegree_compare_2[,2], paired = TRUE, alternative = "two.sided")
t_test_tl <- t.test(trophic_level_compare_2[,4], trophic_level_compare_2[,2], paired = TRUE, alternative = "two.sided")
t_test_closeness <- t.test(closeness_compare_2[,4], closeness_compare_2[,2], paired = TRUE, alternative = "two.sided")
t_test_centrality <- t.test(centrality_compare_2[,4], centrality_compare_2[,2], paired = TRUE, alternative = "two.sided")
t_test_ivi <- t.test(ivi_compare_2[,4], ivi_compare_2[,2], paired = TRUE, alternative = "two.sided")

#Use wilcoxon test - to try and avoid the effect of the large sample size leading to significant differences
#relate with richness and proportion

richness_proportion_indegree <- merge(richness_proportion, indegree_compare_2)
richness_proportion_indegree <- richness_proportion_indegree[,-5]
head(richness_proportion_indegree)
wilcox.test(richness_proportion_indegree[,2], richness_proportion_indegree[,4], paired = TRUE)
wilcox.test(richness_proportion_indegree[,2], richness_proportion_indegree[,5], paired = TRUE)
wilcox.test(richness_proportion_indegree[,3], richness_proportion_indegree[,4], paired = TRUE)
wilcox.test(richness_proportion_indegree[,3], richness_proportion_indegree[,5], paired = TRUE)
#
richness_proportion_outdegree <- merge(richness_proportion, outdegree_compare_2)
richness_proportion_outdegree <- richness_proportion_outdegree[,-5]
head(richness_proportion_outdegree)
wilcox.test(richness_proportion_outdegree[,2], richness_proportion_outdegree[,4], paired = TRUE)
wilcox.test(richness_proportion_outdegree[,2], richness_proportion_outdegree[,5], paired = TRUE)
wilcox.test(richness_proportion_outdegree[,3], richness_proportion_outdegree[,4], paired = TRUE)
wilcox.test(richness_proportion_outdegree[,3], richness_proportion_outdegree[,5], paired = TRUE)
#
richness_proportion_tl <- merge(richness_proportion, trophic_level_compare_2)
richness_proportion_tl <- richness_proportion_tl[,-5]
head(richness_proportion_tl)
wilcox.test(richness_proportion_tl[,2], richness_proportion_tl[,4], paired = TRUE)
wilcox.test(richness_proportion_tl[,2], richness_proportion_tl[,5], paired = TRUE)
wilcox.test(richness_proportion_tl[,3], richness_proportion_tl[,4], paired = TRUE)
wilcox.test(richness_proportion_tl[,3], richness_proportion_tl[,5], paired = TRUE)
#
richness_proportion_closeness <- merge(richness_proportion, closeness_compare_2)
richness_proportion_closeness <- richness_proportion_closeness[,-5]
head(richness_proportion_closeness)
wilcox.test(richness_proportion_closeness[,2], richness_proportion_closeness[,4], paired = TRUE)
wilcox.test(richness_proportion_closeness[,2], richness_proportion_closeness[,5], paired = TRUE)
wilcox.test(richness_proportion_closeness[,3], richness_proportion_closeness[,4], paired = TRUE)
wilcox.test(richness_proportion_closeness[,3], richness_proportion_closeness[,5], paired = TRUE)
#
richness_proportion_centrality <- merge(richness_proportion, centrality_compare_2)
richness_proportion_centrality <- richness_proportion_centrality[,-5]
head(richness_proportion_centrality)
wilcox.test(richness_proportion_centrality[,2], richness_proportion_centrality[,4], paired = TRUE)
wilcox.test(richness_proportion_centrality[,2], richness_proportion_centrality[,5], paired = TRUE)
wilcox.test(richness_proportion_centrality[,3], richness_proportion_centrality[,4], paired = TRUE)
wilcox.test(richness_proportion_centrality[,3], richness_proportion_centrality[,5], paired = TRUE)
#
richness_proportion_ivi <- merge(richness_proportion, ivi_compare_2)
richness_proportion_ivi <- richness_proportion_ivi[,-5]
head(richness_proportion_ivi)
wilcox.test(richness_proportion_ivi[,2], richness_proportion_ivi[,4], paired = TRUE)
wilcox.test(richness_proportion_ivi[,2], richness_proportion_ivi[,5], paired = TRUE)
wilcox.test(richness_proportion_ivi[,3], richness_proportion_ivi[,4], paired = TRUE)
wilcox.test(richness_proportion_ivi[,3], richness_proportion_ivi[,5], paired = TRUE)

#######

wilcox_indeg <- wilcox.test(indegree_compare_2[,4], indegree_compare_2[,2], paired = TRUE)
wilcox_outdeg <- wilcox.test(outdegree_compare_2[,4], outdegree_compare_2[,2], paired = TRUE)
wilcox_tl <- wilcox.test(trophic_level_compare_2[,4], trophic_level_compare_2[,2], paired = TRUE)
wilcox_closeness <- wilcox.test(closeness_compare_2[,4], closeness_compare_2[,2], paired = TRUE)
wilcox_centrality <- wilcox.test(centrality_compare_2[,4], centrality_compare_2[,2], paired = TRUE)
wilcox_ivi <- wilcox.test(ivi_compare_2[,4], ivi_compare_2[,2], paired = TRUE)
#

#Wilcoxon's Effects Size ####################################################### START
library("rstatix")
citation("rstatix")

#Indegree
df1_indegree_compare_3 <- data.frame(indegree_compare_2[,4],rep("T_indegree" , length(indegree_compare_2[,4])))
df2_indegree_compare_3 <- data.frame(indegree_compare_2[,2],rep("NT_indegree" , length(indegree_compare_2[,2])))
colnames(df1_indegree_compare_3) <- c("indegree", "group_T_NT")
colnames(df2_indegree_compare_3) <- c("indegree", "group_T_NT")
indegree_compare_3 <- rbind(df1_indegree_compare_3, df2_indegree_compare_3)
indegree_compare_3$group_T_NT <- as.factor(indegree_compare_3$group_T_NT)
#
wilcox_effsize(indegree ~ group_T_NT, paired = TRUE, data=indegree_compare_3)

#Outdegree
df1_outdegree_compare_3 <- data.frame(outdegree_compare_2[,4],rep("T_indegree" , length(outdegree_compare_2[,4])))
df2_outdegree_compare_3 <- data.frame(outdegree_compare_2[,2],rep("NT_indegree" , length(outdegree_compare_2[,2])))
colnames(df1_outdegree_compare_3) <- c("outdegree", "group_T_NT")
colnames(df2_outdegree_compare_3) <- c("outdegree", "group_T_NT")
outdegree_compare_3 <- rbind(df1_outdegree_compare_3, df2_outdegree_compare_3)
outdegree_compare_3$group_T_NT <- as.factor(outdegree_compare_3$group_T_NT)
#
wilcox_effsize(outdegree ~ group_T_NT, paired = TRUE, data=outdegree_compare_3)

#Trophic level
df1_trophic_level_compare_3 <- data.frame(trophic_level_compare_2[,4],rep("T_tlevel" , length(trophic_level_compare_2[,4])))
df2_trophic_level_compare_3 <- data.frame(trophic_level_compare_2[,2],rep("NT_tlevel" , length(trophic_level_compare_2[,2])))
colnames(df1_trophic_level_compare_3) <- c("tlevel", "group_T_NT")
colnames(df2_trophic_level_compare_3) <- c("tlevel", "group_T_NT")
trophic_level_compare_3 <- rbind(df1_trophic_level_compare_3, df2_trophic_level_compare_3)
trophic_level_compare_3$group_T_NT <- as.factor(trophic_level_compare_3$group_T_NT)
#
wilcox_effsize(tlevel ~ group_T_NT, paired = TRUE, data=trophic_level_compare_3)

#Closeness centrality
df1_closeness_compare_3 <- data.frame(closeness_compare_2[,4],rep("T_tlevel" , length(closeness_compare_2[,4])))
df2_closeness_compare_3 <- data.frame(closeness_compare_2[,2],rep("NT_tlevel" , length(closeness_compare_2[,2])))
colnames(df1_closeness_compare_3) <- c("closeness", "group_T_NT")
colnames(df2_closeness_compare_3) <- c("closeness", "group_T_NT")
closeness_compare_3 <- rbind(df1_closeness_compare_3, df2_closeness_compare_3)
closeness_compare_3$group_T_NT <- as.factor(closeness_compare_3$group_T_NT)
#
wilcox_effsize(closeness ~ group_T_NT, paired = TRUE, data=closeness_compare_3)

#Closeness centrality
df1_centrality_compare_3 <- data.frame(centrality_compare_2[,4],rep("T_tlevel" , length(centrality_compare_2[,4])))
df2_centrality_compare_3 <- data.frame(centrality_compare_2[,2],rep("NT_tlevel" , length(centrality_compare_2[,2])))
colnames(df1_centrality_compare_3) <- c("centrality", "group_T_NT")
colnames(df2_centrality_compare_3) <- c("centrality", "group_T_NT")
centrality_compare_3 <- rbind(df1_centrality_compare_3, df2_centrality_compare_3)
centrality_compare_3$group_T_NT <- as.factor(centrality_compare_3$group_T_NT)
#
wilcox_effsize(centrality ~ group_T_NT, paired = TRUE, data=centrality_compare_3)

#IVI
df1_ivi_compare_3 <- data.frame(ivi_compare_2[,4],rep("T_tlevel" , length(ivi_compare_2[,4])))
df2_ivi_compare_3 <- data.frame(ivi_compare_2[,2],rep("NT_tlevel" , length(ivi_compare_2[,2])))
colnames(df1_ivi_compare_3) <- c("ivi", "group_T_NT")
colnames(df2_ivi_compare_3) <- c("ivi", "group_T_NT")
ivi_compare_3 <- rbind(df1_ivi_compare_3, df2_ivi_compare_3)
ivi_compare_3$group_T_NT <- as.factor(ivi_compare_3$group_T_NT)
#
wilcox_effsize(ivi ~ group_T_NT, paired = TRUE, data=ivi_compare_3)

#Wilcoxon's Effects Size ####################################################### END

#Wilcoxon's Significance ####################################################### START

wilcox_indeg$p.value
wilcox_outdeg$p.value
wilcox_tl$p.value
wilcox_closeness$p.value
wilcox_centrality$p.value
wilcox_ivi$p.value

#Wilcoxon's Significance ####################################################### END

#BOXPLOTS
indegree_compare_3_NT <- data.frame(indegree_compare_2[,2], rep("NT", length = length(indegree_compare_2[,2])))
indegree_compare_3_T <- data.frame(indegree_compare_2[,4], rep("T", length = length(indegree_compare_2[,4])))
names(indegree_compare_3_NT) <- c("indeg", "group")
names(indegree_compare_3_T) <- c("indeg", "group")
indegree_compare_3 <- rbind(indegree_compare_3_T, indegree_compare_3_NT)
ggplot(indegree_compare_3, aes(x=group, y=indeg)) + geom_boxplot(notch = TRUE, outlier.colour="grey")
ggplot(indegree_compare_3, aes(x=group, y=indeg)) +
  geom_violin(aes(fill = group)) +
  guides(fill = "none") +
  labs(title = "Average In-degree") +
  scale_fill_manual(values = c("darkgreen", "darkred"))
#####
outdegree_compare_3_NT <- data.frame(outdegree_compare_2[,2], rep("NT", length = length(outdegree_compare_2[,2])))
outdegree_compare_3_T <- data.frame(outdegree_compare_2[,4], rep("T", length = length(outdegree_compare_2[,4])))
names(outdegree_compare_3_NT) <- c("outdeg", "group")
names(outdegree_compare_3_T) <- c("outdeg", "group")
outdegree_compare_3 <- rbind(outdegree_compare_3_T, outdegree_compare_3_NT)
#ggplot(outdegree_compare_3, aes(x=group, y=outdeg)) + geom_boxplot(notch = TRUE, outlier.colour="grey")
ggplot(outdegree_compare_3, aes(x=group, y=outdeg)) +
  geom_violin(aes(fill = group)) +
  guides(fill = "none") +
  labs(title = "Average Out-degree") +
  scale_fill_manual(values = c("darkgreen", "darkred"))
#####
trophic_level_compare_3_NT <- data.frame(trophic_level_compare_2[,2], rep("NT", length = length(trophic_level_compare_2[,2])))
trophic_level_compare_3_T <- data.frame(trophic_level_compare_2[,4], rep("T", length = length(trophic_level_compare_2[,4])))
names(trophic_level_compare_3_NT) <- c("tl", "group")
names(trophic_level_compare_3_T) <- c("tl", "group")
trophic_level_compare_3 <- rbind(trophic_level_compare_3_T, trophic_level_compare_3_NT)
#ggplot(trophic_level_compare_3, aes(x=group, y=tl)) + geom_boxplot(notch = TRUE, outlier.colour="grey")
ggplot(trophic_level_compare_3, aes(x=group, y=tl)) +
  geom_violin(aes(fill = group)) +
  guides(fill = "none") +
  labs(title = "Average Trophic Level") +
  scale_fill_manual(values = c("darkgreen", "darkred"))
#####
closeness_compare_3_NT <- data.frame(closeness_compare_2[,2], rep("NT", length = length(closeness_compare_2[,2])))
closeness_compare_3_T <- data.frame(closeness_compare_2[,4], rep("T", length = length(closeness_compare_2[,4])))
names(closeness_compare_3_NT) <- c("closeness", "group")
names(closeness_compare_3_T) <- c("closeness", "group")
closeness_compare_3 <- rbind(closeness_compare_3_T, closeness_compare_3_NT)
#ggplot(closeness_compare_3, aes(x=group, y=closeness)) + geom_boxplot(notch = TRUE, outlier.colour="grey")
ggplot(closeness_compare_3, aes(x=group, y=closeness)) +
  geom_violin(aes(fill = group)) +
  guides(fill = "none") +
  labs(title = "Average Closeness Centrality") +
  scale_fill_manual(values = c("darkgreen", "darkred"))
#####
centrality_compare_3_NT <- data.frame(centrality_compare_2[,2], rep("NT", length = length(centrality_compare_2[,2])))
centrality_compare_3_T <- data.frame(centrality_compare_2[,4], rep("T", length = length(centrality_compare_2[,4])))
names(centrality_compare_3_NT) <- c("centrality", "group")
names(centrality_compare_3_T) <- c("centrality", "group")
centrality_compare_3 <- rbind(centrality_compare_3_T, centrality_compare_3_NT)
#ggplot(centrality_compare_3, aes(x=group, y=centrality)) + geom_boxplot(notch = TRUE, outlier.colour="grey")
ggplot(centrality_compare_3, aes(x=group, y=centrality)) +
  geom_violin(aes(fill = group)) +
  guides(fill = "none") +
  labs(title = "Average Betweeness Centrality") +
  scale_fill_manual(values = c("darkgreen", "darkred"))
#####
ivi_compare_3_NT <- data.frame(ivi_compare_2[,2], rep("NT", length = length(ivi_compare_2[,2])))
ivi_compare_3_T <- data.frame(ivi_compare_2[,4], rep("T", length = length(ivi_compare_2[,4])))
names(ivi_compare_3_NT) <- c("ivi", "group")
names(ivi_compare_3_T) <- c("ivi", "group")
ivi_compare_3 <- rbind(ivi_compare_3_T, ivi_compare_3_NT)
#ggplot(ivi_compare_3, aes(x=group, y=ivi)) + geom_boxplot(notch = TRUE, outlier.colour="grey")
ggplot(ivi_compare_3, aes(x=group, y=ivi)) +
  geom_violin(aes(fill = group)) +
  guides(fill = "none") +
  labs(title = "Average IVI") +
  scale_fill_manual(values = c("darkgreen", "darkred"))
###

##########################################################################################
# Identify areas where values in shapefile1 are significantly higher than shapefile2
##########################################################################################

# Set the significance level
significance_level <- 0.05

indegree_shape <- terra::vect("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\indegree_nt_spatial.shp")
indegree_shape <- indegree_shape[, -4]
indegree_shape$indegree_difference <- indegree_difference
#terra::plot(indegree_shape, "indegree_difference", col=rainbow(10))
indegree_shape$significantly_higher <- ifelse(t_test_indeg$p.value < significance_level & indegree_shape$indegree_difference > 0, "Yes", "No")
terra::plot(indegree_shape, "significantly_higher")
#######################
outdegree_shape <- terra::vect("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\outdegree_nt_spatial.shp")
outdegree_shape <- outdegree_shape[, -4]
outdegree_shape$outdegree_difference <- outdegree_difference
#terra::plot(outdegree_shape, "outdegree_difference", col=rainbow(10))
outdegree_shape$significantly_higher <- ifelse(t_test_outdeg$p.value < significance_level & outdegree_shape$outdegree_difference > 0, "Yes", "No")
terra::plot(outdegree_shape, "significantly_higher")
#######################
tl_shape <- terra::vect("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\tl_nt_spatial.shp")
tl_shape <- tl_shape[, -4]
tl_shape$trophic_level_difference <- trophic_level_difference
#terra::plot(tl_shape, "trophic_level_difference", col=rainbow(10))
tl_shape$significantly_higher <- ifelse(t_test_tl$p.value < significance_level & tl_shape$trophic_level_difference > 0, "Yes", "No")
terra::plot(tl_shape, "significantly_higher")
#######################
closeness_shape <- terra::vect("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\closeness_nt_spatial.shp")
closeness_shape <- closeness_shape[, -4]
closeness_shape$closeness_difference <- closeness_difference
#terra::plot(closeness_shape, "closeness_difference", col=rainbow(10))
closeness_shape$significantly_higher <- ifelse(t_test_closeness$p.value < significance_level & closeness_shape$closeness_difference > 0, "Yes", "No")
terra::plot(closeness_shape, "significantly_higher")
#######################
centrality_shape <- terra::vect("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\centrality_nt_spatial.shp")
centrality_shape <- closeness_shape[, -4]
centrality_shape$centrality_difference <- centrality_difference
#terra::plot(centrality_shape, "centrality_difference", col=rainbow(10))
centrality_shape$significantly_higher <- ifelse(t_test_centrality$p.value < significance_level & centrality_shape$centrality_difference > 0, "Yes", "No")
terra::plot(centrality_shape, "significantly_higher")
#######################
ivi_shape <- terra::vect("C:\\Users\\asus\\Documents\\github\\red_listed_networks\\ivi_nt_spatial_second_version.shp")
ivi_shape <- ivi_shape[, -4]
ivi_shape$ivi_difference <- ivi_difference
#terra::plot(ivi_shape, "ivi_difference", col=rainbow(10))
ivi_shape$significantly_higher <- ifelse(t_test_ivi$p.value < significance_level & ivi_shape$ivi_difference > 0, "Yes", "No")
terra::plot(ivi_shape, "significantly_higher")

#Save shapefiles
writeVector(indegree_shape, filename = "indegree_shape.shp")
writeVector(outdegree_shape, filename = "outdegree_shape.shp")
writeVector(tl_shape, filename = "tl_shape.shp")
writeVector(closeness_shape, filename = "closeness_shape.shp")
writeVector(centrality_shape, filename = "centrality_shape.shp")
writeVector(ivi_shape, filename = "ivi_shape.shp")

################################################################################
# USING THE RASTERS - Which is much memory consuming...
################################################################################

nt_indegree <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_indegree.tif")
t_indegree <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_indegree.tif")
#
nt_outdegree <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_outdegree.tif")
t_outdegree <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_outdegree.tif")
#
nt_t_level <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_trophic_level.tif")
t_t_level <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_trophic_level.tif")
#
nt_closeness <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_closeness.tif")
t_closeness <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_closeness.tif")
#
nt_centrality <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_centrality.tif")
t_centrality <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_centrality.tif")
#
nt_ivi <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\nt_ivi.tif")
t_ivi <- terra::rast("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\raster_results_IUCN_NETWORKS\\t_ivi.tif")

#####################################################################################
#                           Keystone Species Index
#####################################################################################

res_test <- c("sp2", "sp3", "sp4", "sp5", "sp5", "sp6", "sp7", "sp7", "sp8")
cons_test <- c("sp1", "sp1", "sp1", "sp3", "sp4", "sp2", "sp2", "sp5", "sp5")
data2 <- data.frame(res_test, cons_test, 1)
names(data2) <- c("resource", "consumer", "link")
data1 <- data2

##

nodeID <- levels(factor(c(as.character(data1[,1]),as.character(data1[,2]))))
numnode <- length(nodeID)
mx <- matrix(rep(0,numnode^2),nrow=numnode,ncol=numnode)
rownames(mx) <- nodeID
colnames(mx) <- nodeID
for (i in 1:length(data1[,1])) mx[as.character(data1[i,1]),as.character(data1[i,2])] <- 1

#Nr of preys per species
prey1 <- numeric(numnode)
for (i in 1:numnode) prey1[i] <- sum(mx[,i])
#Nr of predators per species
predator1 <- numeric(numnode)
for (i in 1:numnode) predator1[i] <- sum(mx[i,])

coef1 <- matrix(rep(0,numnode^2), nrow=numnode, ncol=numnode)
for (i in 1:numnode) coef1[i,] <- prey1*mx[i,]
for (i in 1:numnode) for (j in 1:numnode) {if (coef1[i,j]!=0) coef1[i,j]<-1/coef1[i,j]}
vw <- numeric(numnode)
for (i in 1:numnode) vw[i]<- -1*sum(coef1[i,])
for (i in 1:numnode) coef1[i,i] <- coef1[i,i]+(-1)
kbu <- round(as.vector(MASS::ginv(coef1) %*% vw), 3)
#kbu <- round(solve(coef1,vw),3)

for (i in 1:numnode) {coef1[i,] <- predator1*mx[,i]}
for (i in 1:numnode) for (j in 1:numnode) {if (coef1[i,j]!=0) coef1[i,j]<-1/coef1[i,j]}
vw <- numeric(numnode)
for (i in 1:numnode) vw[i] <- -1*sum(coef1[i,])
for (i in 1:numnode) coef1[i,i] <- coef1[i,i]+(-1)
ktd <- round(as.vector(MASS::ginv(coef1) %*% vw),3)
#ktd <- round(solve(coef1,vw),3)

for (i in 1:numnode) {coef1[i,] <- prey1*mx[i,]}
for (i in 1:numnode) for (j in 1:numnode) {if (coef1[i,j]!=0) coef1[i,j] <- 1/coef1[i,j]}
vw <- numeric(numnode)
for (i in 1:numnode) vw [i]<- sum(coef1[i,])
kdir <- vw
for (i in 1:numnode) {coef1[i,] <- kbu*coef1[i,]}
kindir <- numeric(numnode)
for (i in 1:numnode) {kindir[i] <- sum(coef1[i,])}
for (i in 1:numnode) {coef1[,i] <- predator1*mx[,i]}
for (i in 1:numnode) for (j in 1:numnode) {if (coef1[i,j]!=0) coef1[i,j] <- 1/coef1[i,j]}
vw <- numeric(numnode)
for (i in 1:numnode) vw[i]<-sum(coef1[,i])
kdir <- kdir+vw
for (i in 1:numnode) {coef1[,i] <- ktd*coef1[,i]}
for (i in 1:numnode) {kindir[i] <- kindir[i]+sum(coef1[,i])}

#Derive k
k <- kbu + ktd

#prey2 <- prey1
#predator2 <- predator1

#for(i in 1:length(prey2)) prey2[i] <- ifelse(prey2[i] == 0, 0, 1) 
#for(i in 1:length(predator2)) predator2[i] <- ifelse(predator2[i] == 0, 0, 1) 

resu <- data.frame(nodeID,k,kbu,ktd,kdir,kindir)

#resu$kbu <- resu$kbu * predator2
#resu$ktd <- resu$ktd * prey2

rm(data1,
   nodeID,
   numnode,
   mx,
   prey1,
   predator1,
   coef1,
   vw,
   ktd,
   kbu,
   kdir,
   kindir,
   k,
   resu,
   resu2
)




####################################### CHECK ERROR START

is_dag <- c()
for(m in 1:length(cheddar_list)) {
  is_dag[m] <- igraph::is.dag(igraph_list[[m]])
  message(m)
}

head(which(!is_dag),20)

length(keystone_index)

has_index_k <- c()
has_index_kbu <- c()
has_index_ktd <- c()
has_index_kdir <- c()
has_index_kindir <- c()

for(n in 1:length(keystone_index)){
  df_indx <- keystone_index[[n]]
  df_indx_colnames <- colnames(df_indx)
  has_index_k[n] <- "k" %in% df_indx_colnames
  has_index_kbu[n] <-  "kbu" %in% df_indx_colnames
  has_index_ktd[n] <-  "ktd" %in% df_indx_colnames
  has_index_kdir[n] <-  "kdir" %in% df_indx_colnames
  has_index_kindir[n] <-  "kindir" %in% df_indx_colnames
}

table(has_index_k)
table(has_index_kbu)
table(has_index_ktd)
table(has_index_kdir)
table(has_index_kindir)

####################################### CHECK ERROR END




if(!is.na("kbu") & !is.na("ktd")){
  
  k <- kbu+ktd
  
  kdir <- round(kdir, 3)
  kindir <- round(kindir, 3)
  
  resu <- data.frame(nodeID,k,kbu,ktd,kdir,kindir)
  
  prey2 <- prey1
  predator2 <- predator1
  
  for(i in 1:length(prey2)) prey2[i] <- ifelse(prey2[i] == 0, 0, 1) 
  for(i in 1:length(predator2)) predator2[i] <- ifelse(predator2[i] == 0, 0, 1) 
  
  #resu <- data.frame(nodeID,k,kbu,ktd,kdir,kindir)
  resu$kbu <- resu$kbu * predator2
  resu$ktd <- resu$ktd * prey2
  
  nodes1 <- cheddar_list[[m]]$nodes
  
  resu2 <- merge(x = nodes1,
                 y = resu,
                 by.x = "node",
                 by.y = "nodeID"
  )
  
  keystone_index[[m]] <- resu2
  
}

if(!is.na("kbu") & is.na("ktd")){
  
  #k <- kbu+ktd
  
  #resu <- data.frame(nodeID,k,kbu,ktd,kdir,kindir)
  
  predator2 <- predator1
  for(i in 1:length(predator2)) predator2[i] <- ifelse(predator2[i] == 0, 0, 1) 
  resu <- data.frame(nodeID,kbu)
  resu$kbu <- resu$kbu * predator2
  
  nodes1 <- cheddar_list[[m]]$nodes
  
  resu2 <- merge(x = nodes1,
                 y = resu,
                 by.x = "node",
                 by.y = "nodeID"
  )
  
  keystone_index[[m]] <- resu2
  
}

if(is.na("kbu") & !is.na("ktd")){
  
  #k <- kbu+ktd
  
  #resu <- data.frame(nodeID,k,kbu,ktd,kdir,kindir)
  resu <- data.frame(nodeID,ktd)
  prey2 <- prey1
  for(i in 1:length(prey2)) prey2[i] <- ifelse(prey2[i] == 0, 0, 1) 
  resu$ktd <- resu$ktd * prey2
  
  nodes1 <- cheddar_list[[m]]$nodes
  
  resu2 <- merge(x = nodes1,
                 y = resu,
                 by.x = "node",
                 by.y = "nodeID"
  )
  
  keystone_index[[m]] <- resu2
  
}