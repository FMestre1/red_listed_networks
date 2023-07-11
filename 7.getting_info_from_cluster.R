################################################################################
#
################################################################################

#FMestre
#21-06-2023

library(ggplot2)

load("igraph_node_attrib_df_summarized.RData")
load("igraph_node_attrib_df.RData")

head(igraph_node_attrib_df_summarized)
head(igraph_node_attrib_df)

igraph_node_attrib_df_summarized_2 <- merge(x=igraph_node_attrib_df_summarized, y=igraph_node_attrib_df[,1:4], by.x="name", by.y="name",)
igraph_node_attrib_df_summarized_3 <- unique(igraph_node_attrib_df_summarized_2)

save(igraph_node_attrib_df_summarized_3, file = "igraph_node_attrib_df_summarized_3.RData")

# A basic scatterplot with color depending on Species
names(igraph_node_attrib_df_summarized_3)

ggplot(igraph_node_attrib_df_summarized_3, aes(x=as.factor(status), y=TL, color=body_size)) + 
  geom_point(size=6) +
  theme_ipsum()
