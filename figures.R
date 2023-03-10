################################################################################
# FIGURES
################################################################################

#FMestre
#10-03-2023

#install.packages("treemap")
library(treemap)

View(red_listed_3)
length(fw_list)

red_listed_3[red_listed_3$group == "Amphibians_Reptiles",]

#Sent to Excel to create a pivot table...
write.csv(red_listed_3, "red_listed_3.csv", row.names=FALSE)
pivot_table1 <- read.csv("pivot_table1.csv", sep = ";")

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
