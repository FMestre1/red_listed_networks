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

#Getting Elton traits in...

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
#table(is.na(all_species_status_body_mass_amph_6$unified_bs_2))
#save(all_species_status_body_mass_amph_6, file = "all_species_status_body_mass_amph_6.RData")


  #Fragmentation or habitat structure
