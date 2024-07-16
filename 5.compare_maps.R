
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

