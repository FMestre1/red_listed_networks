################################################################################
#                     Uploading data downloaded from IUCN
################################################################################

#FMestre

#Load packages
library(taxize)

iucn_sept23_aves <- read.csv('C:/Users/asus/Documents/0. Artigos/IUCN_networks/data/iucn_pam_new_data/redlist_species_data_f5ded918-27b1-4050-8a3e-e4d9cbade939/assessments.csv')
iucn_sept23_aves <- iucn_sept23_aves[,3:4]
iucn_sept23_aves_SYN <- read.csv('C:/Users/asus/Documents/0. Artigos/IUCN_networks/data/iucn_pam_new_data/redlist_species_data_f5ded918-27b1-4050-8a3e-e4d9cbade939/synonyms.csv')
iucn_sept23_aves_SYN <- data.frame(iucn_sept23_aves_SYN[,2], paste0(iucn_sept23_aves_SYN[,4], " ", iucn_sept23_aves_SYN[,5]))
names(iucn_sept23_aves_SYN) <- c("speciesname", "synonym")
#
iucn_sept23_rep_mamm_amph <- read.csv("C:/Users/asus/Documents/0. Artigos/IUCN_networks/data/iucn_pam_new_data/redlist_species_data_9e0a257c-d08a-47fa-9bc0-97faa527c323/assessments.csv")
iucn_sept23_rep_mamm_amph <- iucn_sept23_rep_mamm_amph[,3:4]

#remove empty spaces at the end of string
iucn_sept23_aves_SYN[,2] <- trimws(iucn_sept23_aves_SYN[,2])

#Join both
new_species_iucn <- rbind(iucn_sept23_rep_mamm_amph, iucn_sept23_aves)
#write.csv(new_species_iucn, file = "new_species_iucn.csv")

class_species <- c(rep(NA, nrow(iucn_sept23_rep_mamm_amph)), rep("aves", nrow(iucn_sept23_aves)))

for(i in 1:nrow(iucn_sept23_rep_mamm_amph)){
  
tryCatch(
    expr = {zz <- taxize::classification(iucn_sept23_rep_mamm_amph[i,1], get = "class", db = 'itis')
    },
    error = function(e) NULL
  )
  
  
  if(exists("zz")) {if(any(!is.na(zz[[1]]))){ 
  zz <- zz[[1]]
  class_species[i] <- as.character(zz[zz$rank == "class",][1])  
 } else class_species[i] <- NA
   } else class_species[i] <- NA
  
  if(exists("zz")) rm(zz)
   
  message(i)

}

new_species_iucn <- data.frame(new_species_iucn, NA)

new_species_iucn[new_species_iucn$redlistCategory == "Least Concern",][,3] <- "LC"
new_species_iucn[new_species_iucn$redlistCategory == "Endangered",][,3] <- "EN"
new_species_iucn[new_species_iucn$redlistCategory == "Data Deficient",][,3] <- "DD"
new_species_iucn[new_species_iucn$redlistCategory == "Vulnerable",][,3] <- "VU"
new_species_iucn[new_species_iucn$redlistCategory == "Near Threatened",][,3] <- "NT"
new_species_iucn[new_species_iucn$redlistCategory == "Extinct",][,3] <- "EX"
new_species_iucn[new_species_iucn$redlistCategory == "Critically Endangered",][,3] <- "CR"

names(new_species_iucn)[3] <- "code"

head(all_species_status_body_mass_amph_12)

all_species_status_body_mass_amph_13 <- all_species_status_body_mass_amph_12

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

all_species_status_body_mass_amph_13$agreg_ts <- NA

for(i in 1:nrow(all_species_status_body_mass_amph_13)){
  
  cat1 <- all_species_status_body_mass_amph_13[i,]$code
  if(cat1 %in% c("LC", "NT")) all_species_status_body_mass_amph_13$agreg_ts[i] <- "not_threatened"
  if(cat1 %in% c("VU", "EN", "CR")) all_species_status_body_mass_amph_13$agreg_ts[i] <- "threatened"
  if(!(cat1 %in% c("LC", "NT", "VU", "EN", "CR"))) all_species_status_body_mass_amph_13$agreg_ts[i] <- "none"
  message(i)
}

#Save & Load
#save(all_species_status_body_mass_amph_13, file = "all_species_status_body_mass_amph_13.RData")
#load("all_species_status_body_mass_amph_13.RData")
