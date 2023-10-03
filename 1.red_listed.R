################################################################################
#                               Red-listed species 
################################################################################

#FMestre
#04-01-2023

#Source: https://www.eea.europa.eu/data-and-maps/data/european-red-lists-7

#red_listed <- read.csv("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\European_Red_List_2017_December_csv\\European_Red_List_2017_December.csv", 
#                       sep = ",")
                       
#View(red_listed)
#Create vector with genus+species name

#full_name <- c()

#for(i in 1:nrow(red_listed)){
  
  #BIRD
  #if(red_listed$speciesGroup[i] == "Birds"){
  ##genus <- red_listed$taxonomicRankGenus[i]
  ##species1 <- red_listed$taxonomicRankSpecies[i]
  ##full_name[i] <- paste0(genus, "_", species1)
  #}
  
  #MAMMAL
  ##if(red_listed$speciesGroup[i] == "Mammals"){
    ##species2 <- red_listed$scientificName[i]
    ##full_name[i] <- stringr::str_replace(species2, " ", "_")
    
#  }
  
#}  
  
#combine with the previous data frame
#red_listed_2 <- data.frame(full_name, red_listed)
#rownames(red_listed_2) <- 1:nrow(red_listed_2)
#View(red_listed_2)

#Remove unwanted columns
#names(red_listed_2)
#

#red_listed_3 <- data.frame(
#  red_listed_2$speciesGroup,
#  red_listed_2$full_name,
#  red_listed_2$europeanRegionalRedListCategory,
#  red_listed_2$endemicToEurope
#)

names(red_listed_3) <- c("group", "full_name", "europeanRegionalRedListCategory",
                         "endemic_to_europe")

#View(red_listed_3)
#Save
#save(red_listed_3, file = "red_listed_3_08FEV2023.RData")

#Combine with the new dataset to update
#23-09-2023

#red_listed_3$full_name <- stringr::str_replace(red_listed_3$full_name, "_", " ")
#names(red_listed_3)
#names(new_species_iucn)

#merged_update_iucn <- merge(x = red_listed_3, y = new_species_iucn, by.x = "full_name", by.y = "scientificName")
#merged_update_iucn <- merged_update_iucn[,-3]
#head(merged_update_iucn)
#unique(merged_update_iucn$redlistCategory)

head(new_species_iucn)

new_species_iucn <- data.frame(new_species_iucn, NA)

new_species_iucn[new_species_iucn$redlistCategory == "Least Concern",][,3] <- "LC"
new_species_iucn[new_species_iucn$redlistCategory == "Endangered",][,3] <- "EN"
new_species_iucn[new_species_iucn$redlistCategory == "Data Deficient",][,3] <- "DD"
new_species_iucn[new_species_iucn$redlistCategory == "Vulnerable",][,3] <- "VU"
new_species_iucn[new_species_iucn$redlistCategory == "Near Threatened",][,3] <- "NT"
new_species_iucn[new_species_iucn$redlistCategory == "Extinct",][,3] <- "EX"
new_species_iucn[new_species_iucn$redlistCategory == "Critically Endangered",][,3] <- "CR"

names(new_species_iucn)[3] <- "code"
#head(new_species_iucn)
#nrow(new_species_iucn)
#
head(all_species_status_body_mass_amph_12)

#all_species_status_body_mass_amph_12$synonym %in% all_species_status_body_mass_amph_13$species
#

all_species_status_body_mass_amph_13 <- all_species_status_body_mass_amph_12

#Where are synonyms?
#which(all_species_status_body_mass_amph_13$synonym != "NA")


all_species_status_body_mass_amph_13$status <- NA
all_species_status_body_mass_amph_13$code <- NA

for(i in 1:nrow(all_species_status_body_mass_amph_13)){

species0 <- all_species_status_body_mass_amph_13$species[i]
species_syn <- all_species_status_body_mass_amph_13$synonym[i]

status_species0 <- new_species_iucn[new_species_iucn$scientificName == species0,]
status_species_syn <- new_species_iucn[new_species_iucn$scientificName == species_syn,]

if(nrow(status_species0)!=0 && nrow(status_species_syn)!=0 && identical(status_species0, status_species_syn)) new_status1 <- status_species0
if(nrow(status_species0) !=0 && !identical(status_species0, status_species_syn)) new_status1 <- status_species0
if(nrow(status_species_syn) !=0 && !identical(status_species0, status_species_syn)) new_status1 <- status_species_syn
    
all_species_status_body_mass_amph_13$status[i] <- new_status1$redlistCategory
all_species_status_body_mass_amph_13$code[i] <- new_status1$code

}

#unique(new_species_iucn$code)

for(i in 1:nrow(all_species_status_body_mass_amph_13)){
  
  cat1 <- all_species_status_body_mass_amph_13[i,]$code
  if(cat1 %in% c("LC", "NT")) all_species_status_body_mass_amph_13$agreg_ts[i] <- "not_threatened"
  if(cat1 %in% c("VU", "EN", "CR")) all_species_status_body_mass_amph_13$agreg_ts[i] <- "threatened"
  if(!(cat1 %in% c("LC", "NT", "VU", "EN", "CR"))) all_species_status_body_mass_amph_13$agreg_ts[i] <- "none"
  message(i)
}

#table(all_species_status_body_mass_amph_13$agreg_ts)
#nrow(all_species_status_body_mass_amph_12)
#nrow(all_species_status_body_mass_amph_13)
