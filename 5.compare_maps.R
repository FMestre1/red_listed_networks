
################################################################################
#                            COMPARE MAPS - WILCOXON
################################################################################

#FMestre
#15-03-2023

#IVI
ivi_t_spatial$ivi
ivi_nt_spatial$ivi

cor.test(ivi_t_spatial$ivi, ivi_nt_spatial$ivi, method = "spearman", use = "complete.obs")
#
ivi_test <- data.frame(ivi_t_spatial$ivi, ivi_nt_spatial$ivi)
ivi_test <- ivi_test[complete.cases(ivi_test),]
wilcx_ivi <- wilcox.test(ivi_test[,1], ivi_test[,2], paired = TRUE)

#CLOSENESS
closeness_t_spatial$closeness
closeness_nt_spatial$closeness

cor.test(closeness_t_spatial$closeness, closeness_nt_spatial$closeness, method = "spearman", use = "complete.obs")
#
close_test <- data.frame(closeness_t_spatial$closeness, closeness_nt_spatial$closeness)
close_test <- close_test[complete.cases(close_test),]
wilcx_close <- wilcox.test(close_test[,1], close_test[,2], paired = TRUE)

#CENTRALITY
centrality_t_spatial$centrality
centrality_nt_spatial$centrality

cor.test(centrality_t_spatial$centrality, centrality_nt_spatial$centrality, method = "spearman", use = "complete.obs")
#
ctrl_test <- data.frame(centrality_t_spatial$centrality, centrality_nt_spatial$centrality)
ctrl_test <- ctrl_test[complete.cases(ctrl_test),]
ctrl_ivi <- wilcox.test(ctrl_test[,1], ctrl_test[,2], paired = TRUE)


################################################################################
#                                 Effect Size
################################################################################

#Delete everything from the environment
rm(list = ls())

#Load vectors
ivi_nt_spatial <- terra::vect("shape_15JUL24\\ivi_nt_spatial_second_version_15JUL.shp")
ivi_t_spatial <- terra::vect("shape_15JUL24\\ivi_t_spatial_second_version_15JUL.shp")
centrality_nt_spatial <- terra::vect("shape_15JUL24\\centrality_nt_spatial_15JUL.shp")
centrality_t_spatial <- terra::vect("shape_15JUL24\\centrality_t_spatial_15JUL.shp")
indegree_t_spatial <- terra::vect("shape_15JUL24\\indegree_t_spatial_15JUL.shp")
indegree_nt_spatial <- terra::vect("shape_15JUL24\\indegree_nt_spatial_15JUL.shp")
outdegree_t_spatial <- terra::vect("shape_15JUL24\\outdegree_t_spatial_15JUL.shp")
outdegree_nt_spatial <- terra::vect("shape_15JUL24\\outdegree_nt_spatial_15JUL.shp")
closeness_t_spatial <- terra::vect("shape_15JUL24\\closeness_t_spatial_15JUL.shp")
closeness_nt_spatial <- terra::vect("shape_15JUL24\\closeness_nt_spatial_15JUL.shp")
tl_nt_spatial <- terra::vect("shape_15JUL24\\tl_nt_spatial_15JUL.shp")
tl_t_spatial <- terra::vect("shape_15JUL24\\tl_t_spatial_15JUL.shp")

#Coehen's d interpretation
#Small effect size: Cohen's d around 0.2 indicates a small effect size. 
#It suggests that there is a small difference between the groups, which may 
#have limited practical significance.

#Medium effect size: Cohen's d around 0.5 indicates a medium effect size. 
#It suggests a moderate difference between the groups, which is typically 
#considered to have a meaningful practical impact.

#Large effect size: Cohen's d around 0.8 or above indicates a large effect size. 
#It suggests a substantial difference between the groups, which is likely to 
#have significant practical implications.

##### IN-DEGREE #####

indegree_t_spatial_df <- as.data.frame(indegree_t_spatial)
indegree_nt_spatial_df <- as.data.frame(indegree_nt_spatial)
indegree_t_spatial_df <- indegree_t_spatial_df[,c(1,4)]
indegree_nt_spatial_df <- indegree_nt_spatial_df[,c(1,4)]
names(indegree_t_spatial_df)[2] <- "T_indegree"
names(indegree_nt_spatial_df)[2] <- "NT_indegree"
indegree_compare <- merge(indegree_nt_spatial_df, indegree_t_spatial_df)
indegree_compare <- indegree_compare[complete.cases(indegree_compare),]
#View(indegree_compare)
#indegree_ttest <- t.test(indegree_compare[,2], indegree_compare[,3], paired = TRUE)
indegree_cohens_d <- effsize::cohen.d(indegree_compare[,2], indegree_compare[,3])
indegree_wicox <- wilcox.test(indegree_compare[,2], indegree_compare[,3], paired = TRUE, alternative = "two.sided")

##### OUT-DEGREE #####

outdegree_t_spatial_df <- as.data.frame(outdegree_t_spatial)
outdegree_nt_spatial_df <- as.data.frame(outdegree_nt_spatial)
outdegree_t_spatial_df <- outdegree_t_spatial_df[,c(1,4)]
outdegree_nt_spatial_df <- outdegree_nt_spatial_df[,c(1,4)]
names(outdegree_t_spatial_df)[2] <- "T_outdegree"
names(outdegree_nt_spatial_df)[2] <- "NT_outdegree"
outdegree_compare <- merge(outdegree_nt_spatial_df, outdegree_t_spatial_df)
outdegree_compare <- outdegree_compare[complete.cases(outdegree_compare),]
#View(outdegree_compare)
#outdegree_ttest <- t.test(outdegree_compare[,2], outdegree_compare[,3], paired = TRUE)
outdegree_cohens_d <- effsize::cohen.d(outdegree_compare[,2], outdegree_compare[,3])
outdegree_wicox <- wilcox.test(outdegree_compare[,2], outdegree_compare[,3], paired = TRUE, alternative = "two.sided")

##### TROPHIC LEVEL #####

tl_t_spatial_df <- as.data.frame(tl_t_spatial)
tl_nt_spatial_df <- as.data.frame(tl_nt_spatial)
tl_t_spatial_df <- tl_t_spatial_df[,c(1,4)]
tl_nt_spatial_df <- tl_nt_spatial_df[,c(1,4)]
names(tl_t_spatial_df)[2] <- "T_tl"
names(tl_nt_spatial_df)[2] <- "NT_tl"
tl_compare <- merge(tl_nt_spatial_df, tl_t_spatial_df)
tl_compare <- tl_compare[complete.cases(tl_compare),]
#View(tl_compare)
#tl_ttest <- t.test(tl_compare[,2], tl_compare[,3], paired = TRUE)
tl_cohens_d <- effsize::cohen.d(tl_compare[,2], tl_compare[,3])
tl_wicox <- wilcox.test(tl_compare[,2], tl_compare[,3], paired = TRUE, alternative = "two.sided")

##### BETWEENNESS CENTRALITY #####

beet_centrality_t_spatial_df <- as.data.frame(centrality_t_spatial)
beet_centrality_nt_spatial_df <- as.data.frame(centrality_nt_spatial)
beet_centrality_t_spatial_df <- beet_centrality_t_spatial_df[,c(1,4)]
beet_centrality_nt_spatial_df <- beet_centrality_nt_spatial_df[,c(1,4)]
names(beet_centrality_t_spatial_df)[2] <- "T_beet_centrality"
names(beet_centrality_nt_spatial_df)[2] <- "NT_beet_centrality"
beet_centrality_compare <- merge(beet_centrality_nt_spatial_df, beet_centrality_t_spatial_df)
beet_centrality_compare <- beet_centrality_compare[complete.cases(beet_centrality_compare),]
#View(beet_centrality_compare)
#beet_centrality_ttest <- t.test(beet_centrality_compare[,2], beet_centrality_compare[,3], paired = TRUE)
beet_centrality_cohens_d <- effsize::cohen.d(beet_centrality_compare[,2], beet_centrality_compare[,3])
beet_centrality_wicox <- wilcox.test(beet_centrality_compare[,2], beet_centrality_compare[,3], paired = TRUE, alternative = "two.sided")

##### CLOSENESS CENTRALITY #####

closeness_t_spatial_df <- as.data.frame(closeness_t_spatial)
closeness_nt_spatial_df <- as.data.frame(closeness_nt_spatial)
closeness_t_spatial_df <- closeness_t_spatial_df[,c(1,4)]
closeness_nt_spatial_df <- closeness_nt_spatial_df[,c(1,4)]
names(closeness_t_spatial_df)[2] <- "T_closeness"
names(closeness_nt_spatial_df)[2] <- "NT_closeness"
closeness_compare <- merge(closeness_nt_spatial_df, closeness_t_spatial_df)
closeness_compare <- closeness_compare[complete.cases(closeness_compare),]
#View(closeness_compare)
#closeness_ttest <- t.test(closeness_compare[,2], closeness_compare[,3], paired = TRUE)
closeness_cohens_d <- effsize::cohen.d(closeness_compare[,2], closeness_compare[,3])
closeness_wicox <- wilcox.test(closeness_compare[,2], closeness_compare[,3], paired = TRUE, alternative = "two.sided")

##### IVI CENTRALITY #####

ivi_t_spatial_df <- as.data.frame(ivi_t_spatial)
ivi_nt_spatial_df <- as.data.frame(ivi_nt_spatial)
ivi_t_spatial_df <- ivi_t_spatial_df[,c(1,4)]
ivi_nt_spatial_df <- ivi_nt_spatial_df[,c(1,4)]
names(ivi_t_spatial_df)[2] <- "T_ivi"
names(ivi_nt_spatial_df)[2] <- "NT_ivi"
ivi_compare <- merge(ivi_nt_spatial_df, ivi_t_spatial_df)
ivi_compare <- ivi_compare[complete.cases(ivi_compare),]
#View(ivi_compare)
#ivi_ttest <- t.test(ivi_compare[,2], ivi_compare[,3], paired = TRUE)
ivi_cohens_d <- effsize::cohen.d(ivi_compare[,2], ivi_compare[,3])
ivi_wicox <- wilcox.test(ivi_compare[,2], ivi_compare[,3], paired = TRUE, alternative = "two.sided")
