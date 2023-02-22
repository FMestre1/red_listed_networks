#FMestre
#22-02-2023

library(ggplot2)
library(ggplot2)


not_threat <- rbind(lc_table,
                    nt_table)

threatened <- rbind(vu_table,
                    en_table,
                    cr_table)


#in-degree #####################################################################

t_in<-data.frame(rep("threatened", nrow(threatened)), threatened$indegree)
names(t_in) <- c("iucn", "indegree")
head(t_in)
#
nt_in<-data.frame(rep("non-threatened", nrow(not_threat)), not_threat$indegree)
names(nt_in) <- c("iucn", "indegree")
head(nt_in)


INDEGREE_2 <- rbind(t_in,
                    nt_in
                    )

rm(t_in,
   nt_in
)

#save(INDEGREE_2, file = "INDEGREE_2.RData")

INDEGREE_2 <- INDEGREE_2[complete.cases(INDEGREE_2),]
str(INDEGREE_2)
INDEGREE_2$iucn <- as.factor(INDEGREE_2$iucn)


ggplot(INDEGREE_2, aes(x=fct_reorder(iucn,indegree, .desc=TRUE), y=indegree)) +
  ggtitle("In-degree for threatened and non-threatened species") + 
  xlab("IUCN Categories") + ylab("In-degree") +
  geom_boxplot()


#out-degree #####################################################################

t_out<-data.frame(rep("threatened", nrow(threatened)), threatened$outdegree)
names(t_out) <- c("iucn", "outdegree")
head(t_out)
#
nt_out<-data.frame(rep("non-threatened", nrow(not_threat)), not_threat$outdegree)
names(nt_out) <- c("iucn", "outdegree")
head(nt_out)


OUTDEGREE_2 <- rbind(t_out,
                    nt_out
)

rm(t_out,
   nt_out
)

#save(OUTDEGREE_2, file = "OUTDEGREE_2.RData")

OUTDEGREE_2 <- OUTDEGREE_2[complete.cases(OUTDEGREE_2),]
str(OUTDEGREE_2)
OUTDEGREE_2$iucn <- as.factor(OUTDEGREE_2$iucn)


ggplot(OUTDEGREE_2, aes(x=fct_reorder(iucn,outdegree, .desc=TRUE), y=outdegree)) +
  ggtitle("Out-degree for threatened and non-threatened species") + 
  xlab("IUCN Categories") + ylab("Out-degree") +
  geom_boxplot()


#centrality #####################################################################

t_centrality<-data.frame(rep("threatened", nrow(threatened)), threatened$centrality)
names(t_centrality) <- c("iucn", "centrality")
head(t_centrality)
#
nt_centrality<-data.frame(rep("non-threatened", nrow(not_threat)), not_threat$centrality)
names(nt_centrality) <- c("iucn", "centrality")
head(nt_centrality)


CENTRALITY_2 <- rbind(t_centrality,
                      nt_centrality
)

rm(t_centrality,
   nt_centrality
)

#save(CENTRALITY_2, file = "CENTRALITY_2.RData")

CENTRALITY_2 <- CENTRALITY_2[complete.cases(CENTRALITY_2),]
str(CENTRALITY_2)
CENTRALITY_2$iucn <- as.factor(CENTRALITY_2$iucn)


ggplot(CENTRALITY_2, aes(x=fct_reorder(iucn,centrality, .desc=TRUE), y=centrality)) +
  ggtitle("Centrality for threatened and non-threatened species") + 
  xlab("IUCN Categories") + ylab("Centrality") +
  geom_boxplot()


#centrality #####################################################################

t_ivi<-data.frame(rep("threatened", nrow(threatened)), threatened$ivi)
names(t_ivi) <- c("iucn", "ivi")
head(t_ivi)
#
nt_ivi<-data.frame(rep("non-threatened", nrow(not_threat)), not_threat$ivi)
names(nt_ivi) <- c("iucn", "ivi")
head(nt_ivi)


IVI_2 <- rbind(t_ivi,
                      nt_ivi
)

rm(t_ivi,
   nt_ivi
)

#save(IVI_2, file = "IVI_2.RData")

IVI_2 <- IVI_2[complete.cases(IVI_2),]
str(IVI_2)
IVI_2$iucn <- as.factor(IVI_2$iucn)


ggplot(IVI_2, aes(x=fct_reorder(iucn,ivi, .desc=TRUE), y=ivi)) +
  ggtitle("IVI index for threatened and non-threatened species") + 
  xlab("IUCN Categories") + ylab("IVI") +
  geom_boxplot()

#closeness #####################################################################

t_closeness<-data.frame(rep("threatened", nrow(threatened)), threatened$closeness)
names(t_closeness) <- c("iucn", "closeness")
head(t_closeness)
#
nt_closeness<-data.frame(rep("non-threatened", nrow(not_threat)), not_threat$closeness)
names(nt_closeness) <- c("iucn", "closeness")
head(nt_closeness)


CLOSENESS_2 <- rbind(t_closeness,
               nt_closeness
)

rm(t_closeness,
   nt_closeness
)

#save(CLOSENESS_2, file = "CLOSENESS_2.RData")

CLOSENESS_2 <- CLOSENESS_2[complete.cases(CLOSENESS_2),]
str(CLOSENESS_2)
CLOSENESS_2$iucn <- as.factor(CLOSENESS_2$iucn)


ggplot(CLOSENESS_2, aes(x=fct_reorder(iucn,closeness, .desc=TRUE), y=closeness)) +
  ggtitle("Closeness index for threatened and non-threatened species") + 
  xlab("IUCN Categories") + ylab("IVI") +
  geom_boxplot()











