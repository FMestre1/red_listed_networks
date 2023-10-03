################################################################################
#                     Newer data downloaded from IUCN...
################################################################################
#FMestre
#22-09-2023

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
#head(iucn_sept23_rep_mamm_amph_SYN, 20)

#Join both
#head(iucn_sept23_rep_mamm_amph)
#head(iucn_sept23_aves)
new_species_iucn <- rbind(iucn_sept23_rep_mamm_amph, iucn_sept23_aves)
#head(new_species_iucn)
#View(new_species_iucn)

class_species <- c(rep(NA, nrow(iucn_sept23_rep_mamm_amph)), rep("aves", nrow(iucn_sept23_aves)))

for(i in 1:nrow(iucn_sept23_rep_mamm_amph)){
  
  #zz <- taxize::classification(iucn_sept23_rep_mamm_amph[i,1], get = "class", db = 'itis')
  
  
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

#iucn_sept23_rep_mamm_amph[[i-1]]

#AQUI
#new_species_iucn_2 <- data.frame(new_species_iucn, class_species)
#View(new_species_iucn_2)

