################################################################################
#               Fig1 - Proportion of species per threat status             
################################################################################

#FMestre

#Load packages
library(treemap)
library(grDevices)

all_species_status_body_mass_amph_15 <- read.csv("all_species_status_body_mass_amph_15.csv", sep = ";")
head(all_species_status_body_mass_amph_15)
nrow(all_species_status_body_mass_amph_15)
table(all_species_status_body_mass_amph_15$corrected_status)
table(all_species_status_body_mass_amph_15$corrected_agreg_ts)

pivot_table2 <- read.csv("pivot_table_2.csv", sep = ";")

#Plot
png(filename = "tree_2.png",width = 2000, height = 1600)

treemap(pivot_table2, index=c("class","status"), 
        fontsize.labels=c(45,40),  
        vSize="number", 
        type="index",
        border.col=c("black","black"),
        border.lwds=c(2,1),
        bg.labels=0,
        overlap.labels = 0,
        palette = "Set2",# Width of colors
        #title= "Species status per group"
        
)

dev.off()

