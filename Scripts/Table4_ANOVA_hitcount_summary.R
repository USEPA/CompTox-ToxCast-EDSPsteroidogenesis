#-----------------------------------------------------------------------------------#
# Haggard et al. "High -throughput H295R steroidogenesis assay: utility as an
# alternative and a statistical approach to characterize effects on steroidogenesis"

# Summary of ANOVA results for the chemical library

# Table 4: Positive ANOVA results by steroid hormone analyte.

#Original 12 June 2017
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
library(tcpl)
library(stringr)
tcplConf()
tcplConfList()

# setwd() appropriately
# this text file is Supplemental File 5 and contains the hitcall (0,1) for each steroid hormone by chemical sample

anova.result<-fread('Supp5_Global_H295R_ANOVA_OECD_filtered_wide_output_2017-08-09.txt')

#-----------------------------------------------------------------------------------#
# Table 4 Number of chemicals that significantly affect each hormone
#-----------------------------------------------------------------------------------#

anova.result <- data.table(anova.result)
anova.result[, date_chnm_plate := str_replace(date_chnm_plate, "PharmaGSID_", "PharmaGSID")] # manipulation to catch all samples when using the date_chnm_plate string
sample.info<-strsplit(anova.result$date_chnm_plate, split="_", fixed=TRUE)
anova.result$chnm<-paste(sapply(sample.info, "[", 2))

chemdb<-tcplLoadChem() # load chemical information
anova.result$casn<-chemdb$casn[match(anova.result$chnm, chemdb$chnm)]

unique(anova.result[,chnm]) #654 unique chemicals by chnm

#generate Table 4: # and % of chemicals with positive ANOVA result by steroid hormone analyte

anova.sums<-numcolwise(sum)(anova.result)
anova.sums.long<-melt(anova.sums, na.rm=TRUE, value.name='value')
colnames(anova.sums.long)
anova.sums.long<-as.data.table(anova.sums.long)
anova.sums.long[,percent.library:= round((value/654)*100, 1)]

print(anova.sums.long)
