################################################################################
#             Exploratory plotting of these data - IUCN status
################################################################################
#FMestre
#18-02-2023

#Load packages
library(ggplot2)
library(forcats)


#Using these tables:
dd_table
lc_table
nt_table
vu_table
en_table
cr_table
ne_table
re_table

# IVI ##########################################################################
dd_ivi<-data.frame(rep("DD", nrow(dd_table)), dd_table$ivi)
names(dd_ivi) <- c("iucn", "ivi")
head(dd_ivi)
#
lc_ivi<-data.frame(rep("LC", nrow(lc_table)), lc_table$ivi)
names(lc_ivi) <- c("iucn", "ivi")
head(lc_ivi)
#
nt_ivi<-data.frame(rep("NT", nrow(nt_table)), nt_table$ivi)
names(nt_ivi) <- c("iucn", "ivi")
head(nt_ivi)
#
vu_ivi<-data.frame(rep("VU", nrow(vu_table)), vu_table$ivi)
names(vu_ivi) <- c("iucn", "ivi")
head(vu_ivi)
#
en_ivi<-data.frame(rep("EN", nrow(en_table)), en_table$ivi)
names(en_ivi) <- c("iucn", "ivi")
head(en_ivi)
#
cr_ivi<-data.frame(rep("CR", nrow(cr_table)), cr_table$ivi)
names(cr_ivi) <- c("iucn", "ivi")
head(cr_ivi)
#
ne_ivi<-data.frame(rep("NE", nrow(ne_table)), ne_table$ivi)
names(ne_ivi) <- c("iucn", "ivi")
head(ne_ivi)
#
re_ivi<-data.frame(rep("RE", nrow(re_table)), re_table$ivi)
names(re_ivi) <- c("iucn", "ivi")
head(re_ivi)

IVI <- rbind(dd_ivi,
      lc_ivi,
      nt_ivi,
      vu_ivi,
      en_ivi,
      cr_ivi,
      ne_ivi,
      re_ivi)

rm(dd_ivi,
   lc_ivi,
   nt_ivi,
   vu_ivi,
   en_ivi,
   cr_ivi,
   ne_ivi,
   re_ivi)

IVI <- IVI[complete.cases(IVI),]
str(IVI)
IVI$iucn <- as.factor(IVI$iucn)

#save(IVI, file = "IVI.RData")

ggplot(IVI, aes(x=fct_reorder(iucn,ivi, .desc=TRUE), y=ivi)) +
  ggtitle("IVI index by IUCN status") + 
  xlab("IUCN Categories") + ylab("IVI") +
  geom_boxplot()

#in-degree #####################################################################
dd_in<-data.frame(rep("DD", nrow(dd_table)), dd_table$indegree)
names(dd_in) <- c("iucn", "indegree")
head(dd_in)
#
lc_in<-data.frame(rep("LC", nrow(lc_table)), lc_table$indegree)
names(lc_in) <- c("iucn", "indegree")
head(lc_in)
#
nt_in<-data.frame(rep("NT", nrow(nt_table)), nt_table$indegree)
names(nt_in) <- c("iucn", "indegree")
head(nt_in)
#
vu_in<-data.frame(rep("VU", nrow(vu_table)), vu_table$indegree)
names(vu_in) <- c("iucn", "indegree")
head(vu_in)
#
en_in<-data.frame(rep("EN", nrow(en_table)), en_table$indegree)
names(en_in) <- c("iucn", "indegree")
head(en_in)
#
cr_in<-data.frame(rep("CR", nrow(cr_table)), cr_table$indegree)
names(cr_in) <- c("iucn", "indegree")
head(cr_in)
#
ne_in<-data.frame(rep("NE", nrow(ne_table)), ne_table$indegree)
names(ne_in) <- c("iucn", "indegree")
head(ne_in)
#
re_in<-data.frame(rep("RE", nrow(re_table)), re_table$indegree)
names(re_in) <- c("iucn", "indegree")
head(re_in)

INDEGREE <- rbind(dd_in,
             lc_in,
             nt_in,
             vu_in,
             en_in,
             cr_in,
             ne_in,
             re_in)

rm(dd_in,
      lc_in,
      nt_in,
      vu_in,
      en_in,
      cr_in,
      ne_in,
      re_in)


INDEGREE <- INDEGREE[complete.cases(INDEGREE),]
str(INDEGREE)
INDEGREE$iucn <- as.factor(INDEGREE$iucn)

#save(INDEGREE, file = "INDEGREE.RData")

ggplot(INDEGREE, aes(x=fct_reorder(iucn,indegree, .desc=TRUE), y=indegree)) +
  ggtitle("In-degree by IUCN status") + 
  xlab("IUCN Categories") + ylab("In-degree") +
  geom_boxplot()


#out-degree #####################################################################
dd_out<-data.frame(rep("DD", nrow(dd_table)), dd_table$outdegree)
names(dd_out) <- c("iucn", "outdegree")
head(dd_out)
#
lc_out<-data.frame(rep("LC", nrow(lc_table)), lc_table$outdegree)
names(lc_out) <- c("iucn", "outdegree")
head(lc_out)
#
nt_out<-data.frame(rep("NT", nrow(nt_table)), nt_table$outdegree)
names(nt_out) <- c("iucn", "outdegree")
head(nt_out)
#
vu_out<-data.frame(rep("VU", nrow(vu_table)), vu_table$outdegree)
names(vu_out) <- c("iucn", "outdegree")
head(vu_out)
#
en_out<-data.frame(rep("EN", nrow(en_table)), en_table$outdegree)
names(en_out) <- c("iucn", "outdegree")
head(en_out)
#
cr_out<-data.frame(rep("CR", nrow(cr_table)), cr_table$outdegree)
names(cr_out) <- c("iucn", "outdegree")
head(cr_out)
#
ne_out<-data.frame(rep("NE", nrow(ne_table)), ne_table$outdegree)
names(ne_out) <- c("iucn", "outdegree")
head(ne_out)
#
re_out<-data.frame(rep("RE", nrow(re_table)), re_table$outdegree)
names(re_out) <- c("iucn", "outdegree")
head(re_out)

OUTDEGREE <- rbind(dd_out,
                  lc_out,
                  nt_out,
                  vu_out,
                  en_out,
                  cr_out,
                  ne_out,
                  re_out)

rm(dd_out,
   lc_out,
   nt_out,
   vu_out,
   en_out,
   cr_out,
   ne_out,
   re_out)

OUTDEGREE <- OUTDEGREE[complete.cases(OUTDEGREE),]
str(OUTDEGREE)
OUTDEGREE$iucn <- as.factor(OUTDEGREE$iucn)

#save(OUTDEGREE, file = "OUTDEGREE.RData")

ggplot(OUTDEGREE, aes(x=fct_reorder(iucn,outdegree, .desc=TRUE), y=outdegree)) +
  ggtitle("Out-degree by IUCN status") + 
  xlab("IUCN Categories") + ylab("Out-degree") +
  geom_boxplot()

#centrality #####################################################################
dd_centrality<-data.frame(rep("DD", nrow(dd_table)), dd_table$centrality)
names(dd_centrality) <- c("iucn", "centrality")
head(dd_centrality)
#
lc_centrality<-data.frame(rep("LC", nrow(lc_table)), lc_table$centrality)
names(lc_centrality) <- c("iucn", "centrality")
head(lc_centrality)
#
nt_centrality<-data.frame(rep("NT", nrow(nt_table)), nt_table$centrality)
names(nt_centrality) <- c("iucn", "centrality")
head(nt_centrality)
#
vu_centrality<-data.frame(rep("VU", nrow(vu_table)), vu_table$centrality)
names(vu_centrality) <- c("iucn", "centrality")
head(vu_centrality)
#
en_centrality<-data.frame(rep("EN", nrow(en_table)), en_table$centrality)
names(en_centrality) <- c("iucn", "centrality")
head(en_centrality)
#
cr_centrality<-data.frame(rep("CR", nrow(cr_table)), cr_table$centrality)
names(cr_centrality) <- c("iucn", "centrality")
head(cr_centrality)
#
ne_centrality<-data.frame(rep("NE", nrow(ne_table)), ne_table$centrality)
names(ne_centrality) <- c("iucn", "centrality")
head(ne_centrality)
#
re_centrality<-data.frame(rep("RE", nrow(re_table)), re_table$centrality)
names(re_centrality) <- c("iucn", "centrality")
head(re_centrality)

CENTRALITY <- rbind(dd_centrality,
                   lc_centrality,
                   nt_centrality,
                   vu_centrality,
                   en_centrality,
                   cr_centrality,
                   ne_centrality,
                   re_centrality)

rm(dd_centrality,
   lc_centrality,
   nt_centrality,
   vu_centrality,
   en_centrality,
   cr_centrality,
   ne_centrality,
   re_centrality)

CENTRALITY <- CENTRALITY[complete.cases(CENTRALITY),]
str(CENTRALITY)
CENTRALITY$iucn <- as.factor(CENTRALITY$iucn)

#save(CENTRALITY, file = "CENTRALITY.RData")

ggplot(CENTRALITY, aes(x=fct_reorder(iucn,centrality, .desc=TRUE), y=centrality)) +
  ggtitle("Centrality by IUCN status") + 
  xlab("IUCN Categories") + ylab("Centrality") +
  geom_boxplot()


#closeness #####################################################################
dd_closeness<-data.frame(rep("DD", nrow(dd_table)), dd_table$closeness)
names(dd_closeness) <- c("iucn", "closeness")
head(dd_closeness)
#
lc_closeness<-data.frame(rep("LC", nrow(lc_table)), lc_table$closeness)
names(lc_closeness) <- c("iucn", "closeness")
head(lc_closeness)
#
nt_closeness<-data.frame(rep("NT", nrow(nt_table)), nt_table$closeness)
names(nt_closeness) <- c("iucn", "closeness")
head(nt_closeness)
#
vu_closeness<-data.frame(rep("VU", nrow(vu_table)), vu_table$closeness)
names(vu_closeness) <- c("iucn", "closeness")
head(vu_closeness)
#
en_closeness<-data.frame(rep("EN", nrow(en_table)), en_table$closeness)
names(en_closeness) <- c("iucn", "closeness")
head(en_closeness)
#
cr_closeness<-data.frame(rep("CR", nrow(cr_table)), cr_table$closeness)
names(cr_closeness) <- c("iucn", "closeness")
head(cr_closeness)
#
ne_closeness<-data.frame(rep("NE", nrow(ne_table)), ne_table$closeness)
names(ne_closeness) <- c("iucn", "closeness")
head(ne_closeness)
#
re_closeness<-data.frame(rep("RE", nrow(re_table)), re_table$closeness)
names(re_closeness) <- c("iucn", "closeness")
head(re_closeness)

CLOSENESS <- rbind( dd_closeness,
                    lc_closeness,
                    nt_closeness,
                    vu_closeness,
                    en_closeness,
                    cr_closeness,
                    ne_closeness,
                    re_closeness)

rm(dd_closeness,
   lc_closeness,
   nt_closeness,
   vu_closeness,
   en_closeness,
   cr_closeness,
   ne_closeness,
   re_closeness)

CLOSENESS <- CLOSENESS[complete.cases(CLOSENESS),]
str(CLOSENESS)
CLOSENESS$iucn <- as.factor(CLOSENESS$iucn)

#save(CLOSENESS, file = "CLOSENESS.RData")

ggplot(CLOSENESS, aes(x=fct_reorder(iucn,closeness, .desc=TRUE), y=closeness)) +
  ggtitle("Closeness by IUCN status") + 
  xlab("IUCN Categories") + ylab("Closeness") +
  geom_boxplot()

