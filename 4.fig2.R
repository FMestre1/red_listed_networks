################################################################################
#                       Fig2 - Maps of average metrics           
################################################################################

#FMestre
#23-10-2023

#Load packages
library(terra)
library(viridis)
library(rasterVis)
library(gridExtra)
library(grid)

#my.theme <- rasterTheme(YlOrRdTheme(option = "D")(10))

### In-degree ###
nt_indegree_df <- as.data.frame(indegree_nt_spatial)
t_indegree_df <- as.data.frame(indegree_t_spatial)
nt_indegree_df <- nt_indegree_df[,c(1,4)]
t_indegree_df <- t_indegree_df[,c(1,4)]
#
names(nt_indegree_df)[2] <- "NT_indegree"
names(t_indegree_df)[2] <- "T_indegree"
#
indegree_compare <- merge(nt_indegree_df, t_indegree_df)
indegree_compare <- indegree_compare[complete.cases(indegree_compare),]
#View(indegree_compare)

#frequency plots
max(indegree_compare[,2])
max(indegree_compare[,3])

h_in_nt <- hist(indegree_compare[,2])
h_in_t <- hist(indegree_compare[,3])
h_in_nt$counts <- h_in_nt$counts / sum(h_in_nt$counts) * 100
h_in_t$counts <- h_in_t$counts / sum(h_in_t$counts) * 100
#
df_indegreee_plot_nt <- data.frame(h_in_nt$mids, h_in_nt$counts)
names(df_indegreee_plot_nt) <- c("mids", "counts")
#
df_indegreee_plot_t <- data.frame(h_in_t$mids, h_in_t$counts)
names(df_indegreee_plot_t) <- c("mids", "counts")
#
ind_title <- textGrob("Average In-degree", gp = gpar(fontsize = 25, fontface = "bold"))
ind_min_max <- seq(0, 65, length.out = 100)
ind1 <- rasterVis::levelplot(nt_indegree, par.settings=BuRdTheme, at=ind_min_max, main = "Not-threatened", scales = list(draw = FALSE))
ind2 <- rasterVis::levelplot(t_indegree, par.settings=BuRdTheme, at=ind_min_max, main = "Threatened", scales = list(draw = FALSE))
#
indeg_plot <- ggplot() +
  geom_line(data=df_indegreee_plot_t, aes(x=mids, y=counts), color = "darkred", lwd=1) +
  geom_line(data=df_indegreee_plot_nt, aes(x=mids, y=counts), color = "darkgreen", lwd=1) + 
  labs(x = "In-degree", y = "Frequency", title = "In-degree")


### Out-degree ###
nt_outdegree_df <- as.data.frame(outdegree_nt_spatial)
t_outegree_df <- as.data.frame(outdegree_t_spatial)
nt_outdegree_df <- nt_outdegree_df[,c(1,4)]
t_outdegree_df <- t_outegree_df[,c(1,4)]
#
names(nt_outdegree_df)[2] <- "NT_outdegree"
names(t_outdegree_df)[2] <- "T_outdegree"
#
outdegree_compare <- merge(nt_outdegree_df, t_outdegree_df)
outdegree_compare <- outdegree_compare[complete.cases(outdegree_compare),]

#frequency plots
max(outdegree_compare[,2])
max(outdegree_compare[,3])

#names(indegree_compare)
h_out_nt <- hist(outdegree_compare[,2])
h_out_t <- hist(outdegree_compare[,3])
# Convert the counts to percentages
h_out_nt$counts <- h_out_nt$counts / sum(h_out_nt$counts) * 100
h_out_t$counts <- h_out_t$counts / sum(h_out_t$counts) * 100
#
df_outdegreee_plot_nt <- data.frame(h_out_nt$mids, h_out_nt$counts)
names(df_outdegreee_plot_nt) <- c("mids", "counts")
#
df_outdegreee_plot_t <- data.frame(h_out_t$mids, h_out_t$counts)
names(df_outdegreee_plot_t) <- c("mids", "counts")
#
out_title <- textGrob("Average Out-degree", gp = gpar(fontsize = 25, fontface = "bold"))
outd_min_max <- seq(0, 35, length.out = 100)
out1 <- rasterVis::levelplot(nt_outdegree, par.settings=BuRdTheme, at=outd_min_max, main = "Not-threatened", scales = list(draw = FALSE))
out2 <- rasterVis::levelplot(t_outdegree, par.settings=BuRdTheme, at=outd_min_max, main = "Threatened", scales = list(draw = FALSE))
#
outdeg_plot <- ggplot() +
  geom_line(data=df_outdegreee_plot_t, aes(x=mids, y=counts), color = "darkred", lwd=1) +
  geom_line(data=df_outdegreee_plot_nt, aes(x=mids, y=counts), color = "darkgreen", lwd=1) + 
  labs(x = "Out-degree", y = "Frequency", title = "Out-degree")


### Trophic level ###
nt_t_level <- as.data.frame(tl_nt_spatial)
t_t_level <- as.data.frame(tl_t_spatial)
nt_t_level <- nt_t_level[,c(1,4)]
t_t_level <- t_t_level[,c(1,4)]
names(nt_t_level)[2] <- "NT_trph_level"
names(t_t_level)[2] <- "T_trph_level"
trophic_level_compare <- merge(nt_t_level, t_t_level)
trophic_level_compare <- trophic_level_compare[complete.cases(trophic_level_compare),]
names(trophic_level_compare)
#frequency plots
max(trophic_level_compare[,2])
max(trophic_level_compare[,3])
h_tl_nt <- hist(trophic_level_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "Trophic level", breaks = seq(0,4, by=0.25), ylim = c(0,100000))
h_tl_t <- hist(trophic_level_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "Trophic level",  breaks = seq(0,4, by=0.25), ylim = c(0,100000))
# Convert the counts to percentages
h_tl_nt$counts <- h_tl_nt$counts / sum(h_tl_nt$counts) * 100
h_tl_t$counts <- h_tl_t$counts / sum(h_tl_t$counts) * 100
# Create the line plot
df_tl_plot_nt <- data.frame(h_tl_nt$mids, h_tl_nt$counts)
names(df_tl_plot_nt) <- c("mids", "counts")
#
df_tl_plot_t <- data.frame(h_tl_t$mids, h_tl_t$counts)
names(df_tl_plot_t) <- c("mids", "counts")
#
tl_title <- textGrob("Average Trophic level", gp = gpar(fontsize = 25, fontface = "bold"))
tl_min_max <- seq(0, 3, length.out = 100)
tl1 <- rasterVis::levelplot(nt_tl, par.settings=BuRdTheme, at=tl_min_max, main = "Not-threatened", scales = list(draw = FALSE))
tl2 <- rasterVis::levelplot(t_tl, par.settings=BuRdTheme, at=tl_min_max, main = "Threatened", scales = list(draw = FALSE))
#
tl_plot <- ggplot() +
  geom_line(data=df_tl_plot_t, aes(x=mids, y=counts), color = "darkred", lwd=1) +
  geom_line(data=df_tl_plot_nt, aes(x=mids, y=counts), color = "darkgreen", lwd=1) + 
  labs(x = "Trophic level", y = "Frequency", title = "Trophic level")

##### PLOT #####
title_plot1 <- textGrob("Average in-degree", gp = gpar(fontsize = 20, fontface = "bold"))
grid.arrange(ind1, ind2, ncol=2, top=title_plot1)

title_plot2 <- textGrob("Average out-degree", gp = gpar(fontsize = 20, fontface = "bold"))
grid.arrange(out1, out2, ncol=2, top=title_plot2)

title_plot3 <- textGrob("Average trophic level", gp = gpar(fontsize = 20, fontface = "bold"))
grid.arrange(tl1, tl2, ncol=2, top=title_plot3)

grid.arrange(indeg_plot, outdeg_plot, tl_plot, nrow=3)
##### PLOT #####

################################################################################
################################################################################
################################################################################

#IVI
nt_ivi2 <- as.data.frame(ivi_nt_spatial)
t_ivi2 <- as.data.frame(ivi_t_spatial)
nt_ivi2 <- nt_ivi2[,c(1,4)]
t_ivi2 <- t_ivi2[,c(1,4)]
names(nt_ivi2)[2] <- "NT_ivi"
names(t_ivi2)[2] <- "T_ivi"
ivi_compare <- merge(nt_ivi2, t_ivi2)
ivi_compare <- ivi_compare[complete.cases(ivi_compare),]
#
#frequency plots
names(ivi_compare)
max(ivi_compare[,2])
max(ivi_compare[,3])
#
h_ivi_nt <- hist(ivi_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "IVI", breaks = seq(0,100, by=10), ylim = c(0,150000))
h_ivi_t <- hist(ivi_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "IVI",  breaks = seq(0,100, by=10), ylim = c(0,150000))
# Convert the counts to percentages
h_ivi_nt$counts <- h_ivi_nt$counts / sum(h_ivi_nt$counts) * 100
h_ivi_t$counts <- h_ivi_t$counts / sum(h_ivi_t$counts) * 100
#
# Create the line plot
df_ivi_plot_nt <- data.frame(h_ivi_nt$mids, h_ivi_nt$counts)
names(df_ivi_plot_nt) <- c("mids", "counts")
#
df_ivi_plot_t <- data.frame(h_ivi_t$mids, h_ivi_t$counts)
names(df_ivi_plot_t) <- c("mids", "counts")
#
ivi_title <- textGrob("IVI", gp = gpar(fontsize = 25, fontface = "bold"))
ivi_min_max <- seq(0, 100, length.out = 100)
ivi1 <- rasterVis::levelplot(nt_ivi, par.settings=BuRdTheme, at=ivi_min_max, main = "Not-threatened", scales = list(draw = FALSE))
ivi2 <- rasterVis::levelplot(t_ivi, par.settings=BuRdTheme, at=ivi_min_max, main = "Threatened", scales = list(draw = FALSE))
#
ivi_plot <- ggplot() +
  geom_line(data=df_ivi_plot_t, aes(x=mids, y=counts), color = "darkred", lwd=1) +
  geom_line(data=df_ivi_plot_nt, aes(x=mids, y=counts), color = "darkgreen", lwd=1) + 
  labs(x = "IVI", y = "Frequency", title = "IVI")

#Betweenness centrality
nt_centrality <- as.data.frame(centrality_nt_spatial)
t_centrality <- as.data.frame(centrality_t_spatial)
nt_centrality <- nt_centrality[,c(1,4)]
t_centrality <- t_centrality[,c(1,4)]
names(nt_centrality)[2] <- "NT_centrality"
names(t_centrality)[2] <- "T_centrality"
centrality_compare <- merge(nt_centrality, t_centrality)
centrality_compare <- centrality_compare[complete.cases(centrality_compare),]
#frequency plots
names(centrality_compare)
max(centrality_compare[,2])
max(centrality_compare[,3])
#
h_bet_nt <- hist(centrality_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "Betwenness Centrality", breaks = seq(0,460, by=10), ylim = c(0,150000))
h_bet_t <- hist(centrality_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "Betwenness Centrality",  breaks = seq(0,460, by=10), ylim = c(0,150000))
# Convert the counts to percentages
h_bet_nt$counts <- h_bet_nt$counts / sum(h_bet_nt$counts) * 100
h_bet_t$counts <- h_bet_t$counts / sum(h_bet_t$counts) * 100


# Create the line plot
df_btcent_plot_nt <- data.frame(h_bet_nt$mids, h_bet_nt$counts)
names(df_btcent_plot_nt) <- c("mids", "counts")
#
df_btcent_plot_t <- data.frame(h_bet_t$mids, h_bet_t$counts)
names(df_btcent_plot_t) <- c("mids", "counts")
#
btcent_title <- textGrob("Betweenness centrality", gp = gpar(fontsize = 25, fontface = "bold"))
btcent_min_max <- seq(0, 460, length.out = 100)
bt_cent1 <- rasterVis::levelplot(nt_centrality, par.settings=BuRdTheme, at=btcent_min_max, main = "Not-threatened", scales = list(draw = FALSE))
bt_cent2 <- rasterVis::levelplot(t_centrality, par.settings=BuRdTheme, at=btcent_min_max, main = "Threatened", scales = list(draw = FALSE))
#
btcent_plot <- ggplot() +
  geom_line(data=df_btcent_plot_t, aes(x=mids, y=counts), color = "darkred", lwd=1) +
  geom_line(data=df_btcent_plot_nt, aes(x=mids, y=counts), color = "darkgreen", lwd=1) + 
  labs(x = "Betweenness centrality", y = "Frequency", title = "Betweenness centrality")

#Closeness centrality
nt_closeness <- as.data.frame(closeness_nt_spatial)
t_closeness <- as.data.frame(closeness_t_spatial)
nt_closeness <- nt_closeness[,c(1,4)]
t_closeness <- t_closeness[,c(1,4)]
names(nt_closeness)[2] <- "NT_closeness"
names(t_closeness)[2] <- "T_closeness"
closeness_compare <- merge(nt_closeness, t_closeness)
closeness_compare <- closeness_compare[complete.cases(closeness_compare),]
#frequency plots
names(closeness_compare)
max(closeness_compare[,2])
max(closeness_compare[,3])
#
h_clos_nt <- hist(closeness_compare[,2], freq = TRUE, main = "Non-threatened", col = "darkgreen", xlab = "Closeness Centrality", breaks = seq(0,460, by=10), ylim = c(0,150000))
h_clos_t <- hist(closeness_compare[,3], freq = TRUE, main = "Threatened", col = "darkred", xlab = "Closeness Centrality",  breaks = seq(0,460, by=10), ylim = c(0,150000))
# Convert the counts to percentages
h_clos_nt$counts <- h_clos_nt$counts / sum(h_clos_nt$counts) * 100
h_clos_t$counts <- h_clos_t$counts / sum(h_clos_t$counts) * 100
# Create the line plot
df_close_plot_nt <- data.frame(h_clos_nt$mids, h_clos_nt$counts)
names(df_close_plot_nt) <- c("mids", "counts")
#
df_close_plot_t <- data.frame(h_clos_t$mids, h_clos_t$counts)
names(df_close_plot_t) <- c("mids", "counts")
#
close_title <- textGrob("Closeness centrality", gp = gpar(fontsize = 25, fontface = "bold"))
close_min_max <- seq(0, 35, length.out = 100)
close_cent1 <- rasterVis::levelplot(nt_outdegree, par.settings=BuRdTheme, at=close_min_max, main = "Not-threatened", scales = list(draw = FALSE))
close_cent2 <- rasterVis::levelplot(t_outdegree,  par.settings=BuRdTheme, at=close_min_max, main = "Threatened", scales = list(draw = FALSE))
#
close_plot <- ggplot() +
  geom_line(data=df_close_plot_t, aes(x=mids, y=counts), color = "darkred", lwd=1) +
  geom_line(data=df_close_plot_nt, aes(x=mids, y=counts), color = "darkgreen", lwd=1) + 
  labs(x = "Closeness centrality", y = "Frequency", title = "Closeness centrality")

##### PLOT #####
title_plot4 <- textGrob("IVI", gp = gpar(fontsize = 20, fontface = "bold"))
grid.arrange(ivi1, ivi2, ncol=2, top=title_plot4)

title_plot5 <- textGrob("Betweenness centrality", gp = gpar(fontsize = 20, fontface = "bold"))
grid.arrange(bt_cent1, bt_cent2, ncol=2, top=title_plot5)

title_plot6 <- textGrob("Closeness centrality", gp = gpar(fontsize = 20, fontface = "bold"))
grid.arrange(close_cent1, close_cent2, ncol=2, top=title_plot6)

grid.arrange(ivi_plot, btcent_plot, close_plot, nrow=3)
##### PLOT #####
