################################################################################
#                               Red-listed species 
################################################################################

#FMestre
#04-01-2023

#Source: https://www.eea.europa.eu/data-and-maps/data/european-red-lists-7

red_listed <- read.csv("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\European_Red_List_2017_December_csv\\European_Red_List_2017_December.csv", 
                       sep = ",")
                       
#View(red_listed)
#Create vector with genus+species name

full_name <- c()

for(i in 1:nrow(red_listed)){
  
  #BIRD
  #if(red_listed$speciesGroup[i] == "Birds"){
  genus <- red_listed$taxonomicRankGenus[i]
  species1 <- red_listed$taxonomicRankSpecies[i]
  full_name[i] <- paste0(genus, "_", species1)
  #}
  
  #MAMMAL
  if(red_listed$speciesGroup[i] == "Mammals"){
    species2 <- red_listed$scientificName[i]
    full_name[i] <- stringr::str_replace(species2, " ", "_")
    }
  
}  
  
#combine with the previous data frame
red_listed_2 <- data.frame(full_name, red_listed)
rownames(red_listed_2) <- 1:nrow(red_listed_2)
#View(red_listed_2)

#Remove unwanted columns
names(red_listed_2)
#

red_listed_3 <- data.frame(
  red_listed_2$speciesGroup,
  red_listed_2$full_name,
  red_listed_2$europeanRegionalRedListCategory,
  red_listed_2$endemicToEurope
)

names(red_listed_3) <- c("group", "full_name", "europeanRegionalRedListCategory",
                                      "endemic_to_europe")

#View(red_listed_3)
#Save
#save(red_listed_3, file = "red_listed_3_08FEV2023.RData")

#Combine with the new dataset to update
#23-09-2023

red_listed_3$full_name <- stringr::str_replace(red_listed_3$full_name, "_", " ")
#names(red_listed_3)
#names(new_species_iucn)

merged_update_iucn <- merge(x = red_listed_3, y = new_species_iucn, by.x = "full_name", by.y = "scientificName")
merged_update_iucn <- merged_update_iucn[,-3]
#head(merged_update_iucn)
#unique(merged_update_iucn$redlistCategory)

merged_update_iucn[merged_update_iucn$redlistCategory == "Least Concern",][,4] <- "LC"
merged_update_iucn[merged_update_iucn$redlistCategory == "Endangered",][,4] <- "EN"
merged_update_iucn[merged_update_iucn$redlistCategory == "Data Deficient",][,4] <- "DD"
merged_update_iucn[merged_update_iucn$redlistCategory == "Vulnerable",][,4] <- "VU"
merged_update_iucn[merged_update_iucn$redlistCategory == "Near Threatened",][,4] <- "NT"
merged_update_iucn[merged_update_iucn$redlistCategory == "Extinct",][,4] <- "EX"
merged_update_iucn[merged_update_iucn$redlistCategory == "Critically Endangered",][,4] <- "CE"

#head(merged_update_iucn)
