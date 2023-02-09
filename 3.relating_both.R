#FMestre
#08-02-2023

library(taxize)

#Red list data
red_listed_3
red_list_species <- red_listed_3$full_name
red_listed_3$full_name <- stringr::str_replace(red_listed_3$full_name, "_", " ")


#Node metrics across Europe
fw_list
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

#Only those in the red list data were not merged? Strange?
table(species_names4$query %in% red_listed_3$full_name) #species on FW that are in Red List
table(red_listed_3$full_name %in% species_names4$query) #species in red list that are in the FW  

#Ok, we have 190 non-matches!

#Check with synonyms
library(taxize)

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
#the 1 and 5 are solved by syns

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

save(fw_species_with_red_list_status, file = "fw_species_with_red_list_status_09_FEV_2023.RData")


