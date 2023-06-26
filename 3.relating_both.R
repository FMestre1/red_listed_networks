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
red_listed_3$full_name <- stringr::str_replace(red_listed_3$full_name, "_", " ")

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
