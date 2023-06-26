################################################################################
#               Combining the partial tables saved previously
################################################################################
#FMestre
#17-02-2023

#These tables are now in the folder "rdata tables"

#DD table
dd_table_1 <- get(load("dd_table_CORRECTEC_1.RData"))
dd_table_2 <- get(load("dd_table_CORRECTEC_2.RData"))
dd_table_3 <- get(load("dd_table_CORRECTEC_3.RData"))
dd_table_4 <- get(load("dd_table_CORRECTEC_4.RData"))
dd_table_5 <- get(load("dd_table_CORRECTEC_5.RData"))
dd_table_6 <- get(load("dd_table_CORRECTEC_6.RData"))
dd_table_7 <- get(load("dd_table_CORRECTEC_7.RData"))
dd_table_8 <- get(load("dd_table_CORRECTEC_8.RData"))
dd_table_9 <- get(load("dd_table_CORRECTEC_9.RData"))
dd_table_10 <- get(load("dd_table_CORRECTEC_10.RData"))
dd_table_11 <- get(load("dd_table_CORRECTEC_11.RData"))
dd_table_12 <- get(load("dd_table_CORRECTEC_12.RData"))
dd_table_13 <- get(load("dd_table_CORRECTEC_13.RData"))
dd_table_14 <- get(load("dd_table_CORRECTEC_14.RData"))
dd_table_15 <- get(load("dd_table_CORRECTEC_15.RData"))
dd_table_16 <- get(load("dd_table_CORRECTEC_16.RData"))
dd_table_17 <- get(load("dd_table_CORRECTEC_17.RData"))
dd_table_18 <- get(load("dd_table_CORRECTEC_18.RData"))
dd_table_19 <- get(load("dd_table_CORRECTEC_19.RData"))

dd_table <- rbind(dd_table_1, dd_table_2, dd_table_3, dd_table_4, dd_table_5,  
                  dd_table_6, dd_table_7, dd_table_8, dd_table_9, dd_table_10,  
                  dd_table_11, dd_table_12, dd_table_13, dd_table_14, dd_table_15,  
                  dd_table_16, dd_table_17, dd_table_18, dd_table_19)

rm(dd_table_1, dd_table_2, dd_table_3, dd_table_4, dd_table_5,  
   dd_table_6, dd_table_7, dd_table_8, dd_table_9, dd_table_10,  
   dd_table_11, dd_table_12, dd_table_13, dd_table_14, dd_table_15,  
   dd_table_16, dd_table_17, dd_table_18, dd_table_18, dd_table_18)

dd_table <- dd_table[-1,]

View(dd_table)


#LC table
lc_table_1 <- get(load("lc_table_CORRECTEC_1.RData"))
lc_table_2 <- get(load("lc_table_CORRECTEC_2.RData"))
lc_table_3 <- get(load("lc_table_CORRECTEC_3.RData"))
lc_table_4 <- get(load("lc_table_CORRECTEC_4.RData"))
lc_table_5 <- get(load("lc_table_CORRECTEC_5.RData"))
lc_table_6 <- get(load("lc_table_CORRECTEC_6.RData"))
lc_table_7 <- get(load("lc_table_CORRECTEC_7.RData"))
lc_table_8 <- get(load("lc_table_CORRECTEC_8.RData"))
lc_table_9 <- get(load("lc_table_CORRECTEC_9.RData"))
lc_table_10 <- get(load("lc_table_CORRECTEC_10.RData"))
lc_table_11 <- get(load("lc_table_CORRECTEC_11.RData"))
lc_table_12 <- get(load("lc_table_CORRECTEC_12.RData"))
lc_table_13 <- get(load("lc_table_CORRECTEC_13.RData"))
lc_table_14 <- get(load("lc_table_CORRECTEC_14.RData"))
lc_table_15 <- get(load("lc_table_CORRECTEC_15.RData"))
lc_table_16 <- get(load("lc_table_CORRECTEC_16.RData"))
lc_table_17 <- get(load("lc_table_CORRECTEC_17.RData"))
lc_table_18 <- get(load("lc_table_CORRECTEC_18.RData"))
lc_table_19 <- get(load("lc_table_CORRECTEC_19.RData"))

lc_table <- rbind(lc_table_1, lc_table_2, lc_table_3, lc_table_4, lc_table_5,  
                  lc_table_6, lc_table_7, lc_table_8, lc_table_9, lc_table_10,  
                  lc_table_11, lc_table_12, lc_table_13, lc_table_14, lc_table_15,  
                  lc_table_16, lc_table_17, lc_table_18, lc_table_19)

rm(lc_table_1, lc_table_2, lc_table_3, lc_table_4, lc_table_5,  
   lc_table_6, lc_table_7, lc_table_8, lc_table_9, lc_table_10,  
   lc_table_11, lc_table_12, lc_table_13, lc_table_14, lc_table_15,  
   lc_table_16, lc_table_17, lc_table_18, lc_table_19)

lc_table <- lc_table[-1,]

View(lc_table)

#NT table
nt_table_1 <- get(load("nt_table_CORRECTEC_1.RData"))
nt_table_2 <- get(load("nt_table_CORRECTEC_2.RData"))
nt_table_3 <- get(load("nt_table_CORRECTEC_3.RData"))
nt_table_4 <- get(load("nt_table_CORRECTEC_4.RData"))
nt_table_5 <- get(load("nt_table_CORRECTEC_5.RData"))
nt_table_6 <- get(load("nt_table_CORRECTEC_6.RData"))
nt_table_7 <- get(load("nt_table_CORRECTEC_7.RData"))
nt_table_8 <- get(load("nt_table_CORRECTEC_8.RData"))
nt_table_9 <- get(load("nt_table_CORRECTEC_9.RData"))
nt_table_10 <- get(load("nt_table_CORRECTEC_10.RData"))
nt_table_11 <- get(load("nt_table_CORRECTEC_11.RData"))
nt_table_12 <- get(load("nt_table_CORRECTEC_12.RData"))
nt_table_13 <- get(load("nt_table_CORRECTEC_13.RData"))
nt_table_14 <- get(load("nt_table_CORRECTEC_14.RData"))
nt_table_15 <- get(load("nt_table_CORRECTEC_15.RData"))
nt_table_16 <- get(load("nt_table_CORRECTEC_16.RData"))
nt_table_17 <- get(load("nt_table_CORRECTEC_17.RData"))
nt_table_18 <- get(load("nt_table_CORRECTEC_18.RData"))
nt_table_19 <- get(load("nt_table_CORRECTEC_19.RData"))

nt_table <- rbind(nt_table_1, nt_table_2, nt_table_3, nt_table_4, nt_table_5,  
                  nt_table_6, nt_table_7, nt_table_8, nt_table_9, nt_table_10,  
                  nt_table_11, nt_table_12, nt_table_13, nt_table_14, nt_table_15,  
                  nt_table_16, nt_table_17, nt_table_18, nt_table_19)

rm(nt_table_1, nt_table_2, nt_table_3, nt_table_4, nt_table_5,  
   nt_table_6, nt_table_7, nt_table_8, nt_table_9, nt_table_10,  
   nt_table_11, nt_table_12, nt_table_13, nt_table_14, nt_table_15,  
   nt_table_16, nt_table_17, nt_table_18, nt_table_19)

nt_table <- nt_table[-1,]

View(nt_table)

#VU table
vu_table_1 <- get(load("vu_table_CORRECTEC_1.RData"))
vu_table_2 <- get(load("vu_table_CORRECTEC_2.RData"))
vu_table_3 <- get(load("vu_table_CORRECTEC_3.RData"))
vu_table_4 <- get(load("vu_table_CORRECTEC_4.RData"))
vu_table_5 <- get(load("vu_table_CORRECTEC_5.RData"))
vu_table_6 <- get(load("vu_table_CORRECTEC_6.RData"))
vu_table_7 <- get(load("vu_table_CORRECTEC_7.RData"))
vu_table_8 <- get(load("vu_table_CORRECTEC_8.RData"))
vu_table_9 <- get(load("vu_table_CORRECTEC_9.RData"))
vu_table_10 <- get(load("vu_table_CORRECTEC_10.RData"))
vu_table_11 <- get(load("vu_table_CORRECTEC_11.RData"))
vu_table_12 <- get(load("vu_table_CORRECTEC_12.RData"))
vu_table_13 <- get(load("vu_table_CORRECTEC_13.RData"))
vu_table_14 <- get(load("vu_table_CORRECTEC_14.RData"))
vu_table_15 <- get(load("vu_table_CORRECTEC_15.RData"))
vu_table_16 <- get(load("vu_table_CORRECTEC_16.RData"))
vu_table_17 <- get(load("vu_table_CORRECTEC_17.RData"))
vu_table_18 <- get(load("vu_table_CORRECTEC_18.RData"))
vu_table_19 <- get(load("vu_table_CORRECTEC_19.RData"))

vu_table <- rbind(vu_table_1, vu_table_2, vu_table_3, vu_table_4, vu_table_5,  
                  vu_table_6, vu_table_7, vu_table_8, vu_table_9, vu_table_10,  
                  vu_table_11, vu_table_12, vu_table_13, vu_table_14, vu_table_15,  
                  vu_table_16, vu_table_17, vu_table_18, vu_table_19)

rm(vu_table_1, vu_table_2, vu_table_3, vu_table_4, vu_table_5,  
   vu_table_6, vu_table_7, vu_table_8, vu_table_9, vu_table_10,  
   vu_table_11, vu_table_12, vu_table_13, vu_table_14, vu_table_15,  
   vu_table_16, vu_table_17, vu_table_18, vu_table_19)

vu_table <- vu_table[-1,]

View(vu_table)

#EN table
en_table_1 <- get(load("en_table_CORRECTEC_1.RData"))
en_table_2 <- get(load("en_table_CORRECTEC_2.RData"))
en_table_3 <- get(load("en_table_CORRECTEC_3.RData"))
en_table_4 <- get(load("en_table_CORRECTEC_4.RData"))
en_table_5 <- get(load("en_table_CORRECTEC_5.RData"))
en_table_6 <- get(load("en_table_CORRECTEC_6.RData"))
en_table_7 <- get(load("en_table_CORRECTEC_7.RData"))
en_table_8 <- get(load("en_table_CORRECTEC_8.RData"))
en_table_9 <- get(load("en_table_CORRECTEC_9.RData"))
en_table_10 <- get(load("en_table_CORRECTEC_10.RData"))
en_table_11 <- get(load("en_table_CORRECTEC_11.RData"))
en_table_12 <- get(load("en_table_CORRECTEC_12.RData"))
en_table_13 <- get(load("en_table_CORRECTEC_13.RData"))
en_table_14 <- get(load("en_table_CORRECTEC_14.RData"))
en_table_15 <- get(load("en_table_CORRECTEC_15.RData"))
en_table_16 <- get(load("en_table_CORRECTEC_16.RData"))
en_table_17 <- get(load("en_table_CORRECTEC_17.RData"))
en_table_18 <- get(load("en_table_CORRECTEC_18.RData"))
en_table_19 <- get(load("en_table_CORRECTEC_19.RData"))

en_table <- rbind(en_table_1, en_table_2, en_table_3, en_table_4, en_table_5,  
                  en_table_6, en_table_7, en_table_8, en_table_9, en_table_10,  
                  en_table_11, en_table_12, en_table_13, en_table_14, en_table_15,  
                  en_table_16, en_table_17, en_table_18, en_table_19)

rm(en_table_1, en_table_2, en_table_3, en_table_4, en_table_5,  
   en_table_6, en_table_7, en_table_8, en_table_9, en_table_10,  
   en_table_11, en_table_12, en_table_13, en_table_14, en_table_15,  
   en_table_16, en_table_17, en_table_18, en_table_19)

en_table <- en_table[-1,]

View(en_table)

#CR table
cr_table_1 <- get(load("cr_table_CORRECTEC_1.RData"))
cr_table_2 <- get(load("cr_table_CORRECTEC_2.RData"))
cr_table_3 <- get(load("cr_table_CORRECTEC_3.RData"))
cr_table_4 <- get(load("cr_table_CORRECTEC_4.RData"))
cr_table_5 <- get(load("cr_table_CORRECTEC_5.RData"))
cr_table_6 <- get(load("cr_table_CORRECTEC_6.RData"))
cr_table_7 <- get(load("cr_table_CORRECTEC_7.RData"))
cr_table_8 <- get(load("cr_table_CORRECTEC_8.RData"))
cr_table_9 <- get(load("cr_table_CORRECTEC_9.RData"))
cr_table_10 <- get(load("cr_table_CORRECTEC_10.RData"))
cr_table_11 <- get(load("cr_table_CORRECTEC_11.RData"))
cr_table_12 <- get(load("cr_table_CORRECTEC_12.RData"))
cr_table_13 <- get(load("cr_table_CORRECTEC_13.RData"))
cr_table_14 <- get(load("cr_table_CORRECTEC_14.RData"))
cr_table_15 <- get(load("cr_table_CORRECTEC_15.RData"))
cr_table_16 <- get(load("cr_table_CORRECTEC_16.RData"))
cr_table_17 <- get(load("cr_table_CORRECTEC_17.RData"))
cr_table_18 <- get(load("cr_table_CORRECTEC_18.RData"))
cr_table_19 <- get(load("cr_table_CORRECTEC_19.RData"))

cr_table <- rbind(cr_table_1, cr_table_2, cr_table_3, cr_table_4, cr_table_5,  
                  cr_table_6, cr_table_7, cr_table_8, cr_table_9, cr_table_10,  
                  cr_table_11, cr_table_12, cr_table_13, cr_table_14, cr_table_15,  
                  cr_table_16, cr_table_17, cr_table_18, cr_table_19)

rm(cr_table_1, cr_table_2, cr_table_3, cr_table_4, cr_table_5,  
   cr_table_6, cr_table_7, cr_table_8, cr_table_9, cr_table_10,  
   cr_table_11, cr_table_12, cr_table_13, cr_table_14, cr_table_15,  
   cr_table_16, cr_table_17, cr_table_18, cr_table_19)

cr_table <- cr_table[-1,]

View(cr_table)

#NE table
ne_table_1 <- get(load("ne_table_CORRECTEC_1.RData"))
ne_table_2 <- get(load("ne_table_CORRECTEC_2.RData"))
ne_table_3 <- get(load("ne_table_CORRECTEC_3.RData"))
ne_table_4 <- get(load("ne_table_CORRECTEC_4.RData"))
ne_table_5 <- get(load("ne_table_CORRECTEC_5.RData"))
ne_table_6 <- get(load("ne_table_CORRECTEC_6.RData"))
ne_table_7 <- get(load("ne_table_CORRECTEC_7.RData"))
ne_table_8 <- get(load("ne_table_CORRECTEC_8.RData"))
ne_table_9 <- get(load("ne_table_CORRECTEC_9.RData"))
ne_table_10 <- get(load("ne_table_CORRECTEC_10.RData"))
ne_table_11 <- get(load("ne_table_CORRECTEC_11.RData"))
ne_table_12 <- get(load("ne_table_CORRECTEC_12.RData"))
ne_table_13 <- get(load("ne_table_CORRECTEC_13.RData"))
ne_table_14 <- get(load("ne_table_CORRECTEC_14.RData"))
ne_table_15 <- get(load("ne_table_CORRECTEC_15.RData"))
ne_table_16 <- get(load("ne_table_CORRECTEC_16.RData"))
ne_table_17 <- get(load("ne_table_CORRECTEC_17.RData"))
ne_table_18 <- get(load("ne_table_CORRECTEC_18.RData"))
ne_table_19 <- get(load("ne_table_CORRECTEC_19.RData"))

ne_table <- rbind(ne_table_1, ne_table_2, ne_table_3, ne_table_4, ne_table_5,  
                  ne_table_6, ne_table_7, ne_table_8, ne_table_9, ne_table_10,  
                  ne_table_11, ne_table_12, ne_table_13, ne_table_14, ne_table_15,  
                  ne_table_16, ne_table_17, ne_table_18, ne_table_19)

rm(ne_table_1, ne_table_2, ne_table_3, ne_table_4, ne_table_5,  
   ne_table_6, ne_table_7, ne_table_8, ne_table_9, ne_table_10,  
   ne_table_11, ne_table_12, ne_table_13, ne_table_14, ne_table_15,  
   ne_table_16, ne_table_17, ne_table_18, ne_table_19)

ne_table <- ne_table[-1,]

View(ne_table)


#RE tables
re_table_1 <- get(load("re_table_CORRECTEC_1.RData"))
re_table_2 <- get(load("re_table_CORRECTEC_2.RData"))
re_table_3 <- get(load("re_table_CORRECTEC_3.RData"))
re_table_4 <- get(load("re_table_CORRECTEC_4.RData"))
re_table_5 <- get(load("re_table_CORRECTEC_5.RData"))
re_table_6 <- get(load("re_table_CORRECTEC_6.RData"))
re_table_7 <- get(load("re_table_CORRECTEC_7.RData"))
re_table_8 <- get(load("re_table_CORRECTEC_8.RData"))
re_table_9 <- get(load("re_table_CORRECTEC_9.RData"))
re_table_10 <- get(load("re_table_CORRECTEC_10.RData"))
re_table_11 <- get(load("re_table_CORRECTEC_11.RData"))
re_table_12 <- get(load("re_table_CORRECTEC_12.RData"))
re_table_13 <- get(load("re_table_CORRECTEC_13.RData"))
re_table_14 <- get(load("re_table_CORRECTEC_14.RData"))
re_table_15 <- get(load("re_table_CORRECTEC_15.RData"))
re_table_16 <- get(load("re_table_CORRECTEC_16.RData"))
re_table_17 <- get(load("re_table_CORRECTEC_17.RData"))
re_table_18 <- get(load("re_table_CORRECTEC_18.RData"))
re_table_19 <- get(load("re_table_CORRECTEC_19.RData"))

re_table <- rbind(re_table_1, re_table_2, re_table_3, re_table_4, re_table_5,  
                  re_table_6, re_table_7, re_table_8, re_table_9, re_table_10,  
                  re_table_11, re_table_12, re_table_13, re_table_14, re_table_15,  
                  re_table_16, re_table_17, re_table_18, re_table_19)

rm(re_table_1, re_table_2, re_table_3, re_table_4, re_table_5,  
   re_table_6, re_table_7, re_table_8, re_table_9, re_table_10,  
   re_table_11, re_table_12, re_table_13, re_table_14, re_table_15,  
   re_table_16, re_table_17, re_table_18, re_table_19)

re_table <- re_table[-1,]
View(re_table)


#####

dd_table
lc_table
nt_table
vu_table
en_table
cr_table
ne_table
re_table


