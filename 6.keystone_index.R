#FMestre
#31-10-23

#Load packages
library(igraph)
library(cheddar)
library(terra)

#Load geographic information
europeRaster <- terra::rast(x="C:\\Users\\asus\\Documents\\github\\red_listed_networks\\mask10k-20230214T144658Z-001\\mask10k\\reference_grid_10km.img")
cells_info <- foreign::read.dbf(file = "C:/Users/asus/Documents/github/red_listed_networks/mask10k-20230214T144658Z-001/mask10k/reference_grid_10km.img.vat.dbf")
europeRaster_poly <- terra::as.polygons(europeRaster, values = TRUE, extent=FALSE)
europeRaster_poly <- terra::merge(europeRaster_poly, cells_info)

#Loading the newly created lists
#load("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\networks_SET23\\igraph_list_02SET23.RData")
#load("C:\\Users\\asus\\Documents\\0. Artigos\\IUCN_networks\\data\\networks_SET23\\cheddar_list_02SET23.RData")

#load("fw_list_with_status_20OUT.RData")

#Having...
igraph_list
cheddar_list

#Jordan's Keystone Index

#From:
#Jordán, F., Takács-Sánta, A., & Molnar, I. (1999). 
#A reliability theoretical quest for keystones. Oikos, 453-462.

#Adapted from code gently provided by F. Jordán 

igraph1 <- igraph_list[[10]]
cheddar1 <- cheddar_list[[10]]

#as.data.frame(igraph::E(igraph1))

#Plot example networks
#par(mfrow = c(1, 2))    
#plot(igraph1)
#plot(cheddar1)

keystone_index <- list()

for(m in 1:length(cheddar_list)){

data1 <- data.frame(cheddar_list[[m]]$trophic.links,1)
nodeID <- levels(factor(c(as.character(data1[,1]),as.character(data1[,2]))))
numnode <- length(nodeID)
mx <- matrix(rep(0,numnode^2),nrow=numnode,ncol=numnode)
rownames(mx) <- nodeID
colnames(mx) <- nodeID

for (i in 1:length(data1[,1])) mx[as.character(data1[i,1]),as.character(data1[i,2])] <- 1

prey1 <- numeric(numnode)
for (i in 1:numnode) prey1[i] <- sum(mx[,i])

predator1 <- numeric(numnode)
for (i in 1:numnode) predator1[i] <- sum(mx[i,])

coef1 <- matrix(rep(0,numnode^2), nrow=numnode, ncol=numnode)

for (i in 1:numnode) coef1[i,] <- prey1*mx[i,]
for (i in 1:numnode) for (j in 1:numnode) {if (coef1[i,j]!=0) coef1[i,j]<-1/coef1[i,j]}
vw <- numeric(numnode)
for (i in 1:numnode) vw[i]<- -1*sum(coef1[i,])
for (i in 1:numnode) coef1[i,i] <- coef1[i,i]+(-1)

#kbu <- solve(coef1,vw)
kbu <- tryCatch(solve(coef1,vw), error=function(e) NA)
for (i in 1:numnode) {coef1[i,] <- predator1*mx[,i]}
for (i in 1:numnode) for (j in 1:numnode) {if (coef1[i,j]!=0) coef1[i,j]<-1/coef1[i,j]}
vw <- numeric(numnode)
for (i in 1:numnode) vw[i] <- -1*sum(coef1[i,])
for (i in 1:numnode) coef1[i,i] <- coef1[i,i]+(-1)

#ktd <- solve(coef1,vw)
ktd <- tryCatch(solve(coef1,vw), error=function(e) NA)

for (i in 1:numnode) {coef1[i,] <- prey1*mx[i,]}
for (i in 1:numnode) for (j in 1:numnode) {if (coef1[i,j]!=0) coef1[i,j] <- 1/coef1[i,j]}
vw <- numeric(numnode)
for (i in 1:numnode) vw [i]<- sum(coef1[i,])
kdir <- vw
for (i in 1:numnode) {coef1[i,] <- kbu*coef1[i,]}
kindir <- numeric(numnode)
for (i in 1:numnode) {kindir[i] <- sum(coef1[i,])}
for (i in 1:numnode) {coef1[,i] <- predator1*mx[,i]}
for (i in 1:numnode) for (j in 1:numnode) {if (coef1[i,j]!=0) coef1[i,j] <- 1/coef1[i,j]}
vw <- numeric(numnode)
for (i in 1:numnode) vw[i]<-sum(coef1[,i])
kdir <- kdir+vw
for (i in 1:numnode) {coef1[,i] <- ktd*coef1[,i]}
for (i in 1:numnode) {kindir[i] <- kindir[i]+sum(coef1[,i])}

if(!is.na("kbu") & !is.na("ktd")){
  
  k <- kbu+ktd
  
  resu <- data.frame(nodeID,k,kbu,ktd,kdir,kindir)
  nodes1 <- cheddar_list[[m]]$nodes
  
  resu2 <- merge(x = nodes1,
                 y = resu,
                 by.x = "node",
                 by.y = "nodeID"
  )
  
  keystone_index[[m]] <- resu2
  
}

if(!is.na("kbu") & is.na("ktd")){
  
  #k <- kbu+ktd
  
  #resu <- data.frame(nodeID,k,kbu,ktd,kdir,kindir)
  resu <- data.frame(nodeID,kbu)
  
  nodes1 <- cheddar_list[[m]]$nodes
  
  resu2 <- merge(x = nodes1,
                 y = resu,
                 by.x = "node",
                 by.y = "nodeID"
  )
  
  keystone_index[[m]] <- resu2
  
}

if(is.na("kbu") & !is.na("ktd")){
  
  #k <- kbu+ktd
  
  #resu <- data.frame(nodeID,k,kbu,ktd,kdir,kindir)
  resu <- data.frame(nodeID,ktd)
  
  nodes1 <- cheddar_list[[m]]$nodes
  
  resu2 <- merge(x = nodes1,
                 y = resu,
                 by.x = "node",
                 by.y = "nodeID"
  )
  
  keystone_index[[m]] <- resu2
  
}

rm(data1,
   nodeID,
   numnode,
   mx,
   prey1,
   predator1,
   coef1,
   vw,
   ktd,
   kbu,
   kdir,
   kindir,
   k,
   resu,
   resu2
)

message(m)

}

names(keystone_index) <- names(cheddar_list)

keystone_index_2 <- list()
names_keystone_2 <- c()

for(i in 1:length(keystone_index)){

name1 <- names(keystone_index)[i]
names_keystone_2[i] <- name1
fw_list_with_status_net <- fw_list_with_status[[name1]]
  
merge1 <- merge(x = keystone_index[[i]],
       y = fw_list_with_status_net[,c(1,3,4,5,6,7,8,10)],
       by.x = "node",
       by.y = "SP_NAME"
       )

names(merge1)[3] <- "status_code"
names(merge1)[17] <- "agreg_ts"

keystone_index_2[[i]] <- merge1

message(i)

}

names(keystone_index_2) <- names_keystone_2

#Save
#save(keystone_index, file = "keystone_index_02NOV23.RData")
#save(keystone_index_2, file = "keystone_index_2_02NOV23.RData")

#Average indexes for threatened and non-threatened species

keystone_k_nt <- rep(NA, length(keystone_index_2))
keystone_kbu_nt <- rep(NA, length(keystone_index_2))
keystone_ktd_nt <- rep(NA, length(keystone_index_2))
keystone_kdir_nt <- rep(NA, length(keystone_index_2))
keystone_kindir_nt <- rep(NA, length(keystone_index_2))
#
keystone_k_t <- rep(NA, length(keystone_index_2))
keystone_kbu_t <- rep(NA, length(keystone_index_2))
keystone_ktd_t <- rep(NA, length(keystone_index_2))
keystone_kdir_t <- rep(NA, length(keystone_index_2))
keystone_kindir_t <- rep(NA, length(keystone_index_2))

for(i in 1:length(keystone_index_2)){
  
  #grid_name <- names(keystone_index_2)[i]
  tab1 <- keystone_index_2[[i]]
  tab1_nt <- tab1[tab1$agreg_ts == "not_threatened",]
  tab1_t <- tab1[tab1$agreg_ts == "threatened",]
  
  keystone_k_nt[i] <- mean(tab1_nt$k)
  keystone_kbu_nt[i] <- mean(tab1_nt$kbu)
  keystone_ktd_nt[i] <- mean(tab1_nt$ktd)
  keystone_kdir_nt[i] <- mean(tab1_nt$kdir)
  keystone_kindir_nt[i] <- mean(tab1_nt$kindir)
  #
  keystone_k_t[i] <- mean(tab1_t$k)
  keystone_kbu_t[i] <- mean(tab1_t$kbu)
  keystone_ktd_t[i] <- mean(tab1_t$ktd)
  keystone_kdir_t[i] <- mean(tab1_t$kdir)
  keystone_kindir_t[i] <- mean(tab1_t$kindir)
  
  message(i)
  
}

#Combine in a frame
keystone_indexes_df <- data.frame(names_keystone_2,
           keystone_k_nt,
           keystone_kbu_nt,
           keystone_ktd_nt,
           keystone_kdir_nt,
           keystone_kindir_nt,
           keystone_k_t,
           keystone_kbu_t,
           keystone_ktd_t,
           keystone_kdir_t,
           keystone_kindir_t
           )


colnames(keystone_indexes_df) <- c("grid", "k_nt", "kbu_nt", "ktd_nt", "kdir_nt", "kindir_nt", "k_t",
"kbu_t", "ktd_t", "kdir_t", "kindir_t") 

#Create maps

keystone_index_geo <- merge(europeRaster_poly, keystone_indexes_df, by.x = "PageName", by.y = "grid")

writeVector(keystone_index_geo, 
            filename = "shapes_20OUT23\\keystone_indexes_03NOV.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)







