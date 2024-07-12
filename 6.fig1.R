################################################################################
#               Fig1 - Proportion of species per threat status             
################################################################################

#FMestre

#Load packages
library(treemap)

load("all_species_status_body_mass_amph_13.RData")
head(all_species_status_body_mass_amph_13)
nrow(all_species_status_body_mass_amph_13)
table(all_species_status_body_mass_amph_13$status)
table(all_species_status_body_mass_amph_13$agreg_ts)

for(i in 1:nrow(all_species_status_body_mass_amph_13)){
  
  tryCatch(
    expr = {zz <- taxize::classification(all_species_status_body_mass_amph_13[1,i], get = "class", db = 'itis')
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

table(class_species)

#################################

#write.csv(red_listed_3, "red_listed_3.csv", row.names=FALSE)
#red_listed_3 <- read.csv("old_results/red_listed_3.csv", sep = ";")
pivot_table1 <- read.csv("old_results/pivot_table1.csv", sep = ";")
#sum(pivot_table1$count)

pivot_table1[pivot_table1$group == "Amphibians_Reptiles",]$group <- "Amphibians and Reptiles"

#Plot
png(filename = "tree.png",width = 2000, height = 1600)

treemap(pivot_table1, index=c("group","status"), 
        fontsize.labels=c(25,20),  
        vSize="count", 
        type="index",
        border.col=c("black","black"),
        border.lwds=c(2,1),
        bg.labels=0,
        overlap.labels = 0,
        palette = "Set2",# Width of colors
        #title= "Species status per group"
        
)

dev.off()



