################################################################################
#                    RELATING FW METRICS AND FRAGMENTATION
################################################################################

#FMestre
#24-03-2023

#The FW metrics
metrics_dataset_3

#Fragmentation or habitat structure

################################################################################
#                              RELATING WITH BODY SIZE
################################################################################

#BS here is a proxy of dispersal distance and resistance to fragmentation
#Find references for both in the four groups...


#The FW metrics
metrics_dataset_3

#Species names
  #red_listed_3 #These are the threatened
  #not_threat #threatened
  #threatened #non-threatened
  all_species_names <- stringr::str_replace(species_names2, "_", " ")#These are the names in the list of the metrics per network
  all_species_names <- data.frame(all_species_names,1)
  
  red_listed_3$full_name
  
  all_species_status <- merge(x = all_species_names, all.x = TRUE, y = red_listed_3, by.x = "all_species_names", by.y = "full_name")
  View(all_species_status)
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
  View(all_species_status)
  

#Body size of all the species - sources (in kg)
traits <- read.csv("C:\\Users\\FMest\\Documents\\0. Posdoc\\CONTRATO\\species_databases\\AnimalTraits - a curated animal trait database for body mass, metabolic rate and brain size\\observations.csv")
traits <- traits[traits$phylum == "Chordata",]
#unique(traits$class)
#names(traits)
traits <- data.frame(traits$species, traits$body.mass)
#View(traits)
traits <- traits[!is.na(traits$traits.body.mass),]
traits$traits.body.mass*1000


#Now... the painful task of merging both...

#Há um erro qq neste merge.... dá tão poucos matches...

#all_species_status$all_species_names
#traits$traits.species
#all_species_status$all_species_names[!(all_species_status$all_species_names %in% traits$traits.species)]

all_species_status_body_mass <- merge(x = all_species_status, all.x = TRUE, y = traits, by.x = "all_species_names", by.y = "traits.species")
View(all_species_status_body_mass)

#Were there some missing? Which?
all_species_status_body_mass[is.na(all_species_status_body_mass$traits.body.mass),]
nrow(all_species_status_body_mass[is.na(all_species_status_body_mass$traits.body.mass),])
#of
nrow(all_species_status)


##AMPHIBIANS
amph <- read.csv("C:\\Users\\FMest\\Documents\\0. Posdoc\\CONTRATO\\species_databases\\AmphiBIO\\AmphiBIO_v1.csv")
amph <- data.frame(amph$Species, amph$Body_mass_g)
View(amph)

all_species_status_body_mass_amph <- merge(x = all_species_status_body_mass, all.x = TRUE, all.y = FALSE, y = amph, by.x = "all_species_names", by.y = "amph.Species")
#nrow(all_species_status_body_mass)
#nrow(all_species_status_body_mass_amph)

all_species_status_body_mass_amph_2 <- data.frame(all_species_status_body_mass_amph, all_species_status_body_mass$traits.body.mass*1000)

names(all_species_status_body_mass_amph_2) <- c("species", "status", "agreg_ts", "traits.body.mass_kg", "amph_bs_g", "traits.body.mass_g")
all_species_status_body_mass_amph_2 <- all_species_status_body_mass_amph_2[,-4]#remove the animaltraits in kg
#head(all_species_status_body_mass_amph_2)
#table(is.na(all_species_status_body_mass_amph$traits.body.mass))
#table(is.na(all_species_status_body_mass_amph$amph.Body_mass_g))
#all_species_status_body_mass_amph[!is.na(all_species_status_body_mass_amph$amph.Body_mass_g),]
#View(all_species_status_body_mass_amph_2[!is.na(all_species_status_body_mass_amph_2$amph_bs_g) == TRUE & !is.na(all_species_status_body_mass_amph_2$traits.body.mass_g) == TRUE,])
#View(all_species_status_body_mass_amph_2[!is.na(all_species_status_body_mass_amph_2$amph_bs_g),])
#View(all_species_status_body_mass_amph_2[!is.na(all_species_status_body_mass_amph_2$traits.body.mass_g),])

unified_bs <- c()

for(i in 1:nrow(all_species_status_body_mass_amph_2)){
  
  ln1 <- all_species_status_body_mass_amph_2[i,]
  if(!is.na(ln1$traits.body.mass_g)) unified_bs[i] <- ln1$traits.body.mass_g 
  if(!is.na(ln1$amph_bs_g)) unified_bs[i] <- ln1$amph_bs_g 
  if(is.na(ln1$traits.body.mass_g) && is.na(ln1$amph_bs_g)) unified_bs[i] <- NA
  
}

all_species_status_body_mass_amph_3 <- data.frame(all_species_status_body_mass_amph_2[,1:3], unified_bs)#gathering information from animaltraits and amphibio
#table(!is.na(all_species_status_body_mass_amph_3$unified_bs))

#Elton Traits

mamm_elton <- read.delim("C:\\Users\\FMest\\Documents\\0. Posdoc\\CONTRATO\\species_databases\\elton_traits\\MamFuncDat.txt")
bird_elton <- read.delim("C:\\Users\\FMest\\Documents\\0. Posdoc\\CONTRATO\\species_databases\\elton_traits\\BirdFuncDat.txt")
#
mamm_elton <- mamm_elton[, c(2, 24)] 
bird_elton <- bird_elton[, c(8, 36)] 
#
bird_mammal_elton <- rbind(bird_elton, mamm_elton)
nrow(bird_mammal_elton)

elton_vector <- c()

for(i in 1:nrow(all_species_status_body_mass_amph_3)){
  
  species_row <- all_species_status_body_mass_amph_2$species[i]
  if(species_row %in% bird_mammal_elton$Scientific) {
    elton_vector[i] <- bird_mammal_elton[bird_mammal_elton$Scientific == species_row,]$BodyMass.Value
      }else elton_vector[i] <- NA  
  
}

all_species_status_body_mass_amph_4 <- data.frame(all_species_status_body_mass_amph_3, elton_vector)
View(all_species_status_body_mass_amph_4)

unified_bs_2 <- c()

for(i in 1:nrow(all_species_status_body_mass_amph_4)){
  
  linha1 <- all_species_status_body_mass_amph_4[i,]
  if(!is.na(linha1$elton_vector)) unified_bs_2[i] <- linha1$elton_vector 
  if(!is.na(linha1$unified_bs) && is.na(linha1$elton_vector)) unified_bs_2[i] <- linha1$unified_bs 
  if(is.na(linha1$elton_vector) && is.na(linha1$unified_bs)) unified_bs_2[i] <- NA
  
  
}

all_species_status_body_mass_amph_5 <- data.frame(all_species_status_body_mass_amph_4, unified_bs_2)
all_species_status_body_mass_amph_6 <- all_species_status_body_mass_amph_5[,-c(4:5)]
#View(all_species_status_body_mass_amph_6)
#table(is.na(all_species_status_body_mass_amph_6$unified_bs_2)) #TRUE são as que não têm dados de body size!
#save(all_species_status_body_mass_amph_6, file = "all_species_status_body_mass_amph_6.RData")
#Which don'tn have BS info?
#View(all_species_status_body_mass_amph_6[is.na(all_species_status_body_mass_amph_6$unified_bs_2),])


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

meiri_data <- read.csv("C:\\Users\\FMest\\Documents\\0. Posdoc\\CONTRATO\\species_databases\\Meiri_et_al_2021\\meiri_et_al._2021_appendix.csv", sep = ";")
View(meiri_data)
names(meiri_data)
meiri_data <- data.frame(meiri_data$binomial_2020, meiri_data$binomial_.original.files., meiri_data$adult_body_mass..g.)
names(meiri_data) <- c("binomial_2020", "data_binomial_original", "bmass_g")
head(meiri_data)

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

#table(unified_bs_3 == unified_bs_3)

all_species_status_body_mass_amph_7 <- data.frame(all_species_status_body_mass_amph_6, unified_bs_3$bs_A)
head(all_species_status_body_mass_amph_7)

unified_bs_4 <- c()

for(i in 1:nrow(all_species_status_body_mass_amph_7)){
  
  all_species_status_body_mass_amph_7[i,]

  if(!is.na(all_species_status_body_mass_amph_7$unified_bs_2[i])) unified_bs_4[i] <- all_species_status_body_mass_amph_7$unified_bs_2[i]
  if(is.na(all_species_status_body_mass_amph_7$unified_bs_2[i])) unified_bs_4[i] <- all_species_status_body_mass_amph_7$unified_bs_3.bs_A[i]
  
}

all_species_status_body_mass_amph_8 <- data.frame(all_species_status_body_mass_amph_7[,1:3], unified_bs_4)

#sum(is.na(all_species_status_body_mass_amph_8$unified_bs_4))
#Which?

missing_species_bs <- all_species_status_body_mass_amph_8[is.na(all_species_status_body_mass_amph_8$unified_bs_4),]$species


#I have to check this with the synonyms, resorting to taxize
library(taxize)

#Gather all the BS info from previously used sources

#taxize::use_entrez()
#usethis::edit_r_environ()
#ENTREZ_KEY='fafd2118668fc6bacdf37d11c7c1885f5308'#mykey - have to reload R
#all_species_status_body_mass_amph_6[is.na(all_species_status_body_mass_amph_6$unified_bs_2),][,1][1]

syn_list <- rep(NA, length(missing_species_bs))

for(i in 1:length(missing_species_bs)){

species1 <- missing_species_bs[i]

#try(df1 <- get_nbnid(
#  sci_com = species1,
#  ask = TRUE,
#  rank = "species",
#  #rows = 1,
#  rec_only = TRUE
#  ), 
#  silent = TRUE)

#try(df1 <- get_ids(
#  sci = species1,
#  db= "gbif",
#  ask = TRUE,
#  rank = "species",
#  #rows = 1,
#  rec_only = TRUE
#), 
#silent = TRUE)

try(df1 <- get_gbifid(species1, rank="species"),
    silent = TRUE)

syn <- id2name(id = df1[1], db = "gbif")
syn <- syn[[1]]$name
syn_list[i] <- syn

#if(exists("df1")){
#syn1 <- taxize::synonyms(id = df1[1])
#syn1 <- syn1$nameString
#if(is.null(syn1)) syn_list[i] <- "no synonym"
#if(!is.null(syn1)) syn_list[i] <- syn1
#}

#delete
if(exists("df1"))rm(df1)
if(exists("species1"))rm(species1)
if(exists("syn"))rm(syn)

message("########## Did ", i, "! ##########")

}

names(syn_list) <- missing_species_bs

#length(syn_list)
#length(missing_species_bs)
#head(syn_list)
#syn_list == missing_species_bs

syn_list %in% traits$traits.species
syn_list %in% amph$amph.Species
syn_list %in% meiri_data$binomial_2020
syn_list %in% meiri_data$data_binomial_original

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

#sum(is.na(all_species_status_body_mass_amph_8$unified_bs_4))
#sum(is.na(all_species_status_body_mass_amph_9$unified_bs_4))

all_species_status_body_mass_amph_9[is.na(all_species_status_body_mass_amph_9$synonym),]$synonym <- "NA"
#nrow(all_species_status_body_mass_amph_9)
#View(all_species_status_body_mass_amph_9)

################################################################################
## Fragmentation or habitat structure ##########################################
################################################################################

#FMestre
#05-04-2023

library(terra)

grid_europe <- terra::vect("C:/Users/FMest/Documents/0. Artigos/IUCN_networks/shapefiles/grid_10_EUROPE.shp")
crs(grid_europe)

europe <- terra::vect("C:/Users/FMest/Documents/0. Artigos/IUCN_networks/shapefiles/Europe/Europe.shp")
crs(europe)

grid_europe_wgs84 <- terra::project(grid_europe, europe)
crs(grid_europe_wgs84)

europe <- terra::vect("C:/Users/FMest/Documents/0. Artigos/IUCN_networks/shapefiles/Europe/Europe.shp")
crs(europe)

################################################################################
#                                      FGA
################################################################################

FGA2_2009 <- terra::rast("C:\\fw_space\\fragmentation\\eea_r_3035_100_m_fga2-s-2009_p_2009-2016_v01_r00\\FGA2_S_2009_v3.tif")
#plot(FGA2_2009)
#terra::crs(FGA2_2009)

FGA2_2018 <- terra::rast("C:\\fw_space\\fragmentation\\eea_r_3035_100_m_fga2-s-2018_p_2017-2019_v01_r00\\SEFF_2018_10m.tif")
#plot(FGA2_2018)
#terra::crs(FGA2_2018)

FGA2_2009_2 <- terra::crop(FGA2_2009, FGA2_2018)
#terra::ext(FGA2_2009_2)
#terra::ext(FGA2_2018)

FGA2_mean <- terra::mean(FGA2_2009_2, FGA2_2018)
plot(FGA2_mean)
plot(grid_europe, add = TRUE)
crs(FGA2_mean)
#
#FGA2_mean_wgs84 <- terra::project(FGA2_mean, europe)

crs(FGA2_mean_wgs84)

################################################################################
#                                  CORINE
################################################################################

#corine_wgs84 <- terra::rast("C:\\fw_space\\Corine\\corine_5classes_wgs84.tif")
#plot(corine_wgs84)
