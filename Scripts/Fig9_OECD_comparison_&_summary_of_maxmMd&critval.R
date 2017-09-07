#-----------------------------------------------------------------------------------#
# Haggard et al. "Using a high-throughput steroidogenesis assay to predict effects on 
# ...estrogen and androgen production or overall perturbation of the steroid biosynthesis pathway"

# Summary of maximum mean Mahalanobis distance (maxmMd) values.

# Figure 9: Geometric tiling to compare the OECD validation and HT H295R Results.

#Original 09 June 2017
#Updated  10 Aug 2017
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
#loading libraries and sources
#-----------------------------------------------------------------------------------#

rm(list=ls()) 
library(data.table)
library(plyr)
library(reshape2)
library(gplots)
library(RColorBrewer)
library(NMF)
library(ggplot2)
library(grid)
library(gridExtra)
library(tcpl)
tcplConf()
tcplConfList()

# set directory to appropriate location
setwd('./ANOVA/Global ANOVA')
# this RData object contains all of the Mahalanobis distance information

load('../../RData/AllResps2017-08-09.RData')


# this csv contains the manual comparison of the OECD reference chemical performance with the HT H295R performance (using ANOVA-based logic)
OECD_hmap <- fread("LT_HT_OECD_comparisonTable.csv")

#-----------------------------------------------------------------------------------#
# Summary information for the paper on the maxmMd and critical value
#-----------------------------------------------------------------------------------#

# Describe the set of maximum mean Mahalanobis distances (maxmMds)

colnames(Mahalanobis_dat)
#maximum mean Mahalanobis distance = 'maxD11P'

hist(Mahalanobis_dat$maxD11P)

unique(Mahalanobis_dat$casn) #653 unique casns

unique(Mahalanobis_dat$chnm) #654 unique chnms (triadimenol was tested twice in two different concentration ranges)

range(Mahalanobis_dat$maxD11P, na.rm=TRUE)
#0.7333849 50.0743003

median(Mahalanobis_dat$maxD11P, na.rm=TRUE)
#5.128606

mean(Mahalanobis_dat$maxD11P, na.rm=TRUE)
#7.499737

# Describe the set of the critical values (Scrit01 = critical value)

range(Mahalanobis_dat$Scrit01, na.rm=TRUE)
#1.152012 1.807976

median(Mahalanobis_dat$Scrit01, na.rm=TRUE)
#1.642736

mean(Mahalanobis_dat$Scrit01, na.rm=TRUE)
#1.584002

Mahalanobis_dat[,adj.maxmMd := (maxD11P - Scrit01)] # calculate adjusted maxmMd = maxmMd - critical value

range(Mahalanobis_dat$adj.maxmMd, na.rm=TRUE)
#-0.6363836 51.7823943

median(Mahalanobis_dat$adj.maxmMd, na.rm=TRUE)
#3.521399

mean(Mahalanobis_dat$adj.maxmMd, na.rm=TRUE)
#5.915735

#-----------------------------------------------------------------------------------#
# Prepare pathway responses for comaprison in Figure 9 (geometric tiling) 
# using the maxmMd
#-----------------------------------------------------------------------------------#

# align chemical information for the OECD reference chemicals so it can be merged with maxmMd information

casns <- OECD_hmap$CASRN

chem <- tcplLoadChem('casn', casns)

chnms <- unique(chem$chnm)

# extract the mMd information for the OECD reference chemicals

OECD_md<-Mahalanobis_dat[casn %in% casns,c('date_chnm_plate', 
                                           'chnm', 
                                           'casn',
                                           'maxD11P',
                                           'Scrit01')
                         ]
# calculate the adjusted maxmMd and take the log10maxmMd for plotting

OECD_md[,adj.maxmMd := (maxD11P - Scrit01)] # adjusted maxmMd = maxmMd - critical value
is.numeric(OECD_md$adj.maxmMd)
OECD_md[adj.maxmMd<0, adj.maxmMd := NA] # adjusted maxmMd values less than zero are negative, NA
OECD_md[,log.adj.maxmMd := log10(maxD11P - Scrit01)] # take the log10 of the adjusted maxmMd for visualization
OECD_md[log.adj.maxmMd == 'NaN', log.adj.maxmMd:=NA]
OECD_md<-OECD_md[order(-adj.maxmMd)]

# BPA and Danazol were screened more than once and want to take the median of their log10(maxmMd) values for visualization

OECD_md[casn=='80-05-7'] # bisphenol A (BPA)
OECD_md[casn=='17230-88-5'] # Danazol
OECD_md[, median.log.adj.maxmMd := median(log.adj.maxmMd), by=casn]
OECD_md<-OECD_md[!(OECD_md$date_chnm_plate =='20140326_Bisphenol A_Plate20' & casn=='80-05-7'),]
OECD_md<-OECD_md[!(OECD_md$date_chnm_plate =='20140326_Bisphenol A_Plate6' & casn=='80-05-7'),]
OECD_md<-OECD_md[!(OECD_md$date_chnm_plate =='20150610_Danazol_Plate40' & casn=='17230-88-5'),]

unique(OECD_md$chnm) #23 unique chemicals; flutamide and 2,4-dinitrophenol were not screened in concentration-response

# add empty rows for flutamide and 2,4-dinitrophenol

new.row<-data.frame(date_chnm_plate='NA',
              chnm = 'Flutamide',
              casn = '13311-84-7',
              maxD11P = NA,
              Scrit01 = NA,
              adj.maxmMd=NA,
              log.adj.maxmMd=NA,
              median.log.adj.maxmMd=NA)
OECD_md<-rbind.fill(OECD_md, new.row)

new.row2<-data.frame(date_chnm_plate='NA',
                    chnm = '2,4-Dinitrophenol',
                    casn = '51-28-5',
                    maxD11P = NA,
                    Scrit01 = NA,
                    adj.maxmMd=NA,
                    log.adj.maxmMd=NA,
                    median.log.adj.maxmMd=NA)
OECD_md<-rbind.fill(OECD_md, new.row2)

fwrite(OECD_md, 'OECD_maxmMd_values_10Aug2017.csv')

#-----------------------------------------------------------------------------------#
#create simple matrix for geometric tiling
#-----------------------------------------------------------------------------------#

colnames(OECD_hmap)
colnames(OECD_hmap)[2]<-'chnm'
OECD_hmap$mMD<-OECD_md$median.log.adj.maxmMd[match(OECD_hmap$CASRN, OECD_md$casn)]
OECD_hmap<-OECD_hmap[order(-mMD)]

hmap<-OECD_hmap[,!c(1:2,8:12,17:20)]

hmatrix<-data.matrix(hmap)
rownames(hmatrix) <- OECD_hmap[,chnm]

#-----------------------------------------------------------------------------------#
# Draw Figure 9
#-----------------------------------------------------------------------------------#

#setwd('L:/Lab/NCCT_ToxCast/Derik Haggard/OECD Analysis/Supp Tables or Figures 2/Figs')

pdf(file="Fig9_Geometric_tiling_comparison_18Aug2017.pdf",
    height = 8,
    width = 8)

annotation=data.frame(log10maxmMd = OECD_md$median.log.adj.maxmMd)
aheatmap(hmatrix, 
         color=c("#999999", "white", "#009E73", "#F0E442", "#0072B2"), 
         breaks=1,
         border_color = "grey",
         cellwidth = 15,
         cellheight = 15,
         scale ='none',
         cexRow=1,
         cexCol=0.75,
         Rowv=NA,
         Colv=NA,
         legend=FALSE,
         #distfun ="euclidean", 
         #hclustfun="complete",
         fontsize = 12,
         #txt=hmap_mat_text,
         annRow=annotation,
         annLegend=TRUE,
         annColors="YlOrRd"
         )

dev.off()



