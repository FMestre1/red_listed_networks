################################################################################
#                     Newer data downloaded from IUCN...
################################################################################
#FMestre
#22-09-2023


iucn_sept23_rep_mamm_amph <- read.csv('C:/Users/asus/Documents/0. Artigos/IUCN_networks/data/iucn_pam_new_data/redlist_species_data_f5ded918-27b1-4050-8a3e-e4d9cbade939/assessments.csv')
#View(iucn_sept23_rep_mamm_amph)
iucn_sept23_rep_mamm_amph <- iucn_sept23_rep_mamm_amph[,3:4]
#head(iucn_sept23_rep_mamm_amph)
iucn_sept23_rep_mamm_amph_SYN <- read.csv('C:/Users/asus/Documents/0. Artigos/IUCN_networks/data/iucn_pam_new_data/redlist_species_data_f5ded918-27b1-4050-8a3e-e4d9cbade939/synonyms.csv')
#View(iucn_sept23_rep_mamm_amph_SYN)
iucn_sept23_rep_mamm_amph_SYN <- data.frame(iucn_sept23_rep_mamm_amph_SYN[,2], paste0(iucn_sept23_rep_mamm_amph_SYN[,4], " ", iucn_sept23_rep_mamm_amph_SYN[,5]))
names(iucn_sept23_rep_mamm_amph_SYN) <- c("speciesname", "synonym")
#head(iucn_sept23_rep_mamm_amph_SYN, 20)
#
iucn_sept23_aves <- read.csv("C:/Users/asus/Documents/0. Artigos/IUCN_networks/data/iucn_pam_new_data/redlist_species_data_9e0a257c-d08a-47fa-9bc0-97faa527c323/assessments.csv")
#View(iucn_sept23_aves)
iucn_sept23_aves <- iucn_sept23_aves[,3:4]
#head(iucn_sept23_aves)

#remove empty spaces at the end of string
iucn_sept23_rep_mamm_amph_SYN[,2] <- trimws(iucn_sept23_rep_mamm_amph_SYN[,2])
#head(iucn_sept23_rep_mamm_amph_SYN, 20)

#Join both
#head(iucn_sept23_rep_mamm_amph)
#head(iucn_sept23_aves)
new_species_iucn <- rbind(iucn_sept23_rep_mamm_amph, iucn_sept23_aves)
#head(new_species_iucn)
#View(new_species_iucn)
