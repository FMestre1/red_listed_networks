
################################################################################
#                              #STANDARDIZE MAPS
################################################################################

#FMestre

#Load packages
library(vegan)
library(terra)

###########
#   IVI
###########

#NT

ivi_nt_std_vector <- as.vector(vegan::decostand(ivi_nt_spatial$ivi, method = "standardize", na.rm = TRUE))
ivi_nt_spatial_STD <- ivi_nt_spatial
ivi_nt_spatial_STD$ivi_STD <- ivi_nt_std_vector
#
writeVector(ivi_nt_spatial_STD, 
            filename = "shape_15JUL24\\ivi_nt_spatial_STD_15JUL24.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#T

ivi_t_std_vector <- as.vector(vegan::decostand(ivi_t_spatial$ivi, method = "standardize", na.rm = TRUE))
ivi_t_spatial_STD <- ivi_t_spatial
ivi_t_spatial_STD$ivi_STD <- ivi_t_std_vector
#
writeVector(ivi_t_spatial_STD, 
            filename = "shape_15JUL24\\ivi_t_spatial_STD_15JUL24.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

###########
#CENTRALITY
###########

#NT

centrality_nt_std_vector <- as.vector(vegan::decostand(centrality_nt_spatial$centrality, method = "standardize", na.rm = TRUE))
centrality_nt_std_vector_STD <- centrality_nt_spatial
centrality_nt_std_vector_STD$centrality_STD <- centrality_nt_std_vector
#
writeVector(centrality_nt_std_vector_STD, 
            filename = "shape_15JUL24\\centrality_nt_std_vector_STD_15JUL24.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#T

centrality_t_std_vector <- as.vector(vegan::decostand(centrality_t_spatial$centrality, method = "standardize", na.rm = TRUE))
centrality_t_std_vector_STD <- centrality_t_spatial
centrality_t_std_vector_STD$centrality_STD <- centrality_t_std_vector
#
writeVector(centrality_t_std_vector_STD, 
            filename = "shape_15JUL24\\centrality_t_std_vector_STD_15JUL24.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


###########
#IN-DEGREE
###########

#NT

indegree_nt_std_vector <- as.vector(vegan::decostand(indegree_nt_spatial$indegree, method = "standardize", na.rm = TRUE))
indegree_nt_std_vector_STD <- indegree_nt_spatial
indegree_nt_std_vector_STD$indegree_STD <- indegree_nt_std_vector
#
writeVector(indegree_nt_std_vector_STD, 
            filename = "shape_15JUL24\\indegree_nt_std_vector_STD_15JUL24.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#T

indegree_t_std_vector <- as.vector(vegan::decostand(indegree_t_spatial$indegree, method = "standardize", na.rm = TRUE))
indegree_t_std_vector_STD <- indegree_t_spatial
indegree_t_std_vector_STD$indegree_STD <- indegree_t_std_vector
#
writeVector(indegree_t_std_vector_STD, 
            filename = "shape_15JUL24\\indegree_t_std_vector_STD_15JUL24.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


###########
#OUT-DEGREE
###########

#NT

outdegree_nt_std_vector <- as.vector(vegan::decostand(outdegree_nt_spatial$outdegree, method = "standardize", na.rm = TRUE))
outdegree_nt_std_vector_STD <- outdegree_nt_spatial
outdegree_nt_std_vector_STD$outdegree_STD <- outdegree_nt_std_vector
#
writeVector(outdegree_nt_std_vector_STD, 
            filename = "shape_15JUL24\\outdegree_nt_std_vector_STD_15JUL24.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#T

outdegree_t_std_vector <- as.vector(vegan::decostand(outdegree_t_spatial$outdegree, method = "standardize", na.rm = TRUE))
outdegree_t_std_vector_STD <- outdegree_t_spatial
outdegree_t_std_vector_STD$outdegree_STD <- outdegree_t_std_vector
#
writeVector(outdegree_t_std_vector_STD, 
            filename = "shape_15JUL24\\outdegree_t_std_vector_STD_15JUL24.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)


###########
#CLOSENESS
###########

#NT

closeness_nt_std_vector <- as.vector(vegan::decostand(closeness_nt_spatial$closeness, method = "standardize", na.rm = TRUE))
closeness_nt_std_vector_STD <- closeness_nt_spatial
closeness_nt_std_vector_STD$closeness_STD <- closeness_nt_std_vector
#
writeVector(closeness_nt_std_vector_STD, 
            filename = "shape_15JUL24\\closeness_nt_std_vector_STD_15JUL24.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

#T

closeness_t_std_vector <- as.vector(vegan::decostand(closeness_t_spatial$closeness, method = "standardize", na.rm = TRUE))
closeness_t_std_vector_STD <- closeness_t_spatial
closeness_t_std_vector_STD$closeness_STD <- closeness_t_std_vector
#
writeVector(closeness_t_std_vector_STD, 
            filename = "shape_15JUL24\\closeness_t_std_vector_STD_15JUL24.shp",
            filetype=NULL, 
            layer=NULL, 
            insert=FALSE,
            overwrite=TRUE, 
            options="ENCODING=UTF-8"
)

