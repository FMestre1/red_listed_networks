################################################################################
#                               Red-listed species 
################################################################################

#FMestre
#04-01-2023

#Source: https://www.eea.europa.eu/data-and-maps/data/european-red-lists-7

red_listed <- read.csv("C:\\Users\\FMest\\Documents\\0. Artigos\\IUCN_networks\\data\\European_Red_List_2017_December_csv\\European_Red_List_2017_December.csv", 
                       sep = ","
                       )
#head(red_listed)
#View(red_listed)
#unique(red_listed$speciesGroup)

#red_listed_mammals_birds <- red_listed[red_listed$speciesGroup == c("Birds", "Mammals"),]
#unique(red_listed_mammals_birds$speciesGroup)
#nrow(red_listed_mammals_birds)

#View(red_listed_mammals_birds)
#nrow(red_listed_mammals_birds)

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
View(red_listed_2)

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
#save(red_listed_3, file = "red_listed_3_08FEV2023.RData")



################################################################################