################################################################################
#                         Ploting FW for the figure
################################################################################

# Get a good example to the fw figure

for(i in 1:length(network_list_igraph_2)){
  
  if(all(c("Lynx pardinus", "Oryctolagus cuniculus", "Vulpes vulpes") %in% igraph::vertex.attributes(network_list_igraph_2[[i]])[1]$name)) print(i)
  #if(any(c("Lynx pardinus", "Oryctolagus cuniculus") %in% igraph::vertex.attributes(network_list_igraph_2[[i]])[1])) print (i)
  #if(any(c("Ursus maritimus") %in% igraph::vertex.attributes(network_list_igraph_2[[i]])[1]$name)) print(i)
  
}

for(i in 1:length(network_list_cheddar)){
  
  if(all(c("Lynx pardinus", "Oryctolagus cuniculus", "Vulpes vulpes") %in%  network_list_cheddar[[i]]$nodes$node)) print(i)
  #if(any(c("Lynx pardinus", "Oryctolagus cuniculus") %in% igraph::vertex.attributes(network_list_igraph_2[[i]])[1])) print (i)
  #if(any(c("Ursus maritimus") %in% igraph::vertex.attributes(network_list_igraph_2[[i]])[1]$name)) print(i)
  
}


#cheddar::TrophicLinkPropertyNames (network_list_cheddar[[116183]])
#cheddar::NodePropertyNames (network_list_cheddar[[116183]])

#Which is the most vulnerable prey
sort(cheddar::TrophicVulnerability(network_list_cheddar[[116183]]))
"Pelophylax perezi" # this is it, but also use...
"Oryctolagus cuniculus"

#Which is the most generalist predator
sort(cheddar::TrophicGenerality(network_list_cheddar[[116183]]))
"Vulpes vulpes" # this is it, but also use...
"Lynx pardinus"

#Plot links going and coming from Vulpes vulpes
links <- cbind(TLPS(network_list_cheddar[[116183]]), colour="#c7c7c788")
#
links <- cbind(TLPS(network_list_cheddar[[116183]]), colour="#ffffff00")
links$colour["Vulpes vulpes" == links$resource] <- "red"
links$colour["Vulpes vulpes" == links$consumer] <- "darkgreen"

#cheddar::plot.Community(network_list_cheddar[[116183]], node.labels="node", show.nodes.as="both", link.col=links$colour)
cheddar::plot.Community(network_list_cheddar[[116183]], link.col=links$colour)
#

#Plot links going and coming from Lynx pardinus
links0 <- cbind(TLPS(network_list_cheddar[[116183]]), colour="#c7c7c788")
#
links0 <- cbind(TLPS(network_list_cheddar[[116183]]), colour="#ffffff00")
links0$colour["Lynx pardinus" == links0$resource] <- "red"
links0$colour["Lynx pardinus" == links0$consumer] <- "darkgreen"

#cheddar::plot.Community(network_list_cheddar[[116183]], node.labels="node", show.nodes.as="both", link.col=links$colour)
cheddar::plot.Community(network_list_cheddar[[116183]], link.col=links0$colour)
#

#Plot links coming from Oryctolagus cuniculus
# transparent #ffffff00
links1 <- cbind(TLPS(network_list_cheddar[[116183]]), colour="#c7c7c788")
#
links1 <- cbind(TLPS(network_list_cheddar[[116183]]), colour="#ffffff00")
links1$colour["Oryctolagus cuniculus" == links1$resource] <- "red"
links1$colour["Oryctolagus cuniculus" == links1$consumer] <- "darkgreen"

#cheddar::plot.Community(network_list_cheddar[[116183]], node.labels="node", show.nodes.as="both", link.col=links$colour)
cheddar::plot.Community(network_list_cheddar[[116183]], link.col=links1$colour)
#

############################################################

nodes1 <- cbind(NPS(network_list_cheddar[[116183]]), colour="darkgreen")
#
#nodes1 <- cbind(NPS(network_list_cheddar[[116183]]), colour="#ffffff00")

nodes1$colour["threatened" == nodes1$agreg_ts] <- "red"
nodes1$colour["not_threatened" == nodes1$agreg_ts] <- "darkgreen"

#cheddar::plot.Community(network_list_cheddar[[116183]], node.labels="node", show.nodes.as="both", link.col=links$colour)

pch_vector <- rep(19, length=length(nodes1$colour))
pch_cex <- rep(1, length=length(nodes1$colour))

for(i in 1:length(pch_cex)) if(nodes1$colour[i] == "red") pch_cex[i] <- 2

#The links are setup as transparent: link.col="#ffffff00"
cheddar::plot.Community(network_list_cheddar[[116183]], link.col="#ffffff00", col = nodes1$colour, pch = 19, cex = pch_cex)
#
