#-----------------------------------------------------------------------------------#
#Need to understand how many chemicals are positive at each hormone
#Katie Paul Friedman, paul-friedman.katie@epa.gov
#Original 12 June 2017
#Updated  1 Aug 2017
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
#loading libraries and setting directories
#-----------------------------------------------------------------------------------#

rm(list=ls()) 
library(data.table)
library(plyr)
library(reshape2)
library(gplots)
library(RColorBrewer)
library(NMF)
library(ggplot2)
library(tcpl)
library(stringr)
tcplConf()
tcplConfList()
getwd()
setwd('L:/Lab/NCCT_ToxCast/Derik Haggard/OECD Analysis/R Project/ANOVA/Global ANOVA')
anova.result<-fread('Global_H295R_ANOVA_OECD_filtered_wide_output_2017-07-31.txt')
load('Global_ANOVA_06122017.RData')

#-----------------------------------------------------------------------------------#
#understand the number of hits by chemical - Table 4 in paper
#-----------------------------------------------------------------------------------#

#first get the anova results and chemical information lined up in same table
colnames(anova.result)
head(anova.result)
anova.result <- data.table(anova.result)
anova.result[, date_chnm_plate := str_replace(date_chnm_plate, "PharmaGSID_", "PharmaGSID")]
sample.info<-strsplit(anova.result$date_chnm_plate, split="_", fixed=TRUE)
anova.result$chnm<-paste(sapply(sample.info, "[", 2))
chemdb<-tcplLoadChem()
anova.result$casn<-chemdb$casn[match(anova.result$chnm, chemdb$chnm)]

unique(anova.result[,chnm])
#654 unique chemicals by chnm

#generate Table of positive ANOVA # and % by steroid hormone analyte (Table 4)
anova.sums<-numcolwise(sum)(anova.result)
anova.sums.long<-melt(anova.sums, na.rm=TRUE, value.name='value')
colnames(anova.sums.long)
anova.sums.long<-as.data.table(anova.sums.long)
anova.sums.long[,percent.library:= round((value/655)*100, 1)]

#-----------------------------------------------------------------------------------#
#recall sensitivity assessment
#-----------------------------------------------------------------------------------#

aeid.table<-tcplLoadAeid(fld='asid', val=8)
aeid.table<-aeid.table[!grep("noMTC", aenm)]
aeid.table<-aeid.table[!grep("MTT", aenm)]
aeids<-aeid.table$aeid
aeids

acid.table<-tcplLoadAcid(fld='asid', val=8)
acid.table<-acid.table[!grep("noMTC", acnm)]
acid.table<-acid.table[!grep("MTT", acnm)]
acids<-acid.table$acid
acids

sc2<-tcplLoadData(lvl=2, fld='aeid', val=aeids, type='sc')

sc2[aeid %in% c(890,891),Hormone := '11-deoxycortisol']
sc2[aeid %in% c(892,893),Hormone := 'OH-Pregnenolone']
sc2[aeid %in% c(894,895),Hormone := 'OH-Progesterone']
sc2[aeid %in% c(896,897),Hormone := 'Androstenedione']
sc2[aeid %in% c(898,899),Hormone := 'Corticosterone']
sc2[aeid %in% c(900,901),Hormone := 'Cortisol']
sc2[aeid %in% c(902,903),Hormone := 'DHEA']
sc2[aeid %in% c(904,905),Hormone := 'DOC']
sc2[aeid %in% c(906,907),Hormone := 'Estradiol']
sc2[aeid %in% c(908,909),Hormone := 'Estrone']
sc2[aeid %in% c(910,911),Hormone := 'Pregnenolone']
sc2[aeid %in% c(912,913),Hormone := 'Progesterone']
sc2[aeid %in% c(914,915),Hormone := 'Testosterone']
head(sc2)
sc2$casn<-chemdb$casn[match(sc2$spid, chemdb$spid)]
sc2$chnm<-chemdb$chnm[match(sc2$spid, chemdb$spid)]
sc2$chid<-chemdb$chid[match(sc2$spid, chemdb$spid)]
sc2[,horm.hitc := max(hitc), by=list(casn, Hormone)]
colnames(sc2)

#now create unique data.frame by chemical and make it wide
sc2<-unique(sc2[, list(chid, casn, chnm, Hormone, horm.hitc)])
sc2<-sc2[order(chnm)]
unique(sc2$chnm)
#2012 entries by chnm (substract 3 for FOR, PRO, DMSO)
colnames(anova.result)
sc2<-sc2[chnm %in% anova.result$chnm]
unique(sc2$chnm)
#only 638 of the 655 unique chems in anova.result were tested in single conc

sc2.wide<-dcast(sc2, chnm + chid + casn ~ Hormone, value.var='horm.hitc')
sc2.wide<-as.data.frame(sc2.wide)

sc2.wide<-sc2.wide[!sc2.wide$chnm %in% c('PRO',
                     'FOR',
                     'DMSO',
                     'NA'),]
sc2.wide<-as.data.frame(sc2.wide)
colnames(sc2.wide)
sc2.wide[is.na(sc2.wide)]<-0
sc2.wide<-as.data.table(sc2.wide)
sc2.wide[, hitcsum := rowSums(.SD), .SDcols = colnames(sc2.wide)[4:13]]
colnames(anova.result)
anova.result[,hitcsum.anova := rowSums(.SD), .SDcols = colnames(anova.result)[2:12]]

anova.result$hitcsum.sc<-sc2.wide$hitcsum[match(anova.result$casn, sc2.wide$casn)]

#this sort of starts to fall apart because we have different numbers of hormones
#not sure how to calculate the recall sensitivity

#want to know how many hit 3 or more hormones   
anova.result.unique<-anova.result[hitcsum.sc>2] 
colnames(anova.result.unique)
anova.result.unique[,c(1:12, 15):=NULL]
setkey(anova.result.unique, 'chnm')
anova.result.unique<-subset(unique(anova.result.unique))
anova.result.unique
#518 chemicals out of 638
518/638

#this wouldn't include flag logic
oecd.anova.result<-anova.result[casn %in% casns]
colnames(oecd.anova.result)
oecd.anova.result<-oecd.anova.result[,c('Androstenedione',
                     'Testosterone',
                     'Estrone',
                     'Estradiol',
                     'date_chnm_plate',
                     'chnm',
                     'casn')]
oecd.anova.result<-oecd.anova.result[order(chnm)]
