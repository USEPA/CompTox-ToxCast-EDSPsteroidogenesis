#-----------------------------------------------------------------------------------#
# Haggard et al. "Using a high-throughput steroidogenesis assay to predict effects on 
# ...estrogen and androgen production or overall perturbation of the steroid biosynthesis pathway"

# Figure 10 Boxplot and scatterplot of ANOVA hit calls and maxmMd values

# This figure is to show the distribution of the adjusted maxmMd values for all tested samples

#updated 08/11/2017

rm(list = ls())

library(reshape2)
library(grid)
library(gridExtra)
library(ggplot2)
library(GGally)
library(stringr)
library(data.table)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(ggiraphExtra)
library(corrplot)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(raster)

setwd("../Mahalanobis Distance")

#-------------------------------------------------------------------------------------------------#
#-----Load in all of the data from the mahalanobis_script.R output RData file
#-------------------------------------------------------------------------------------------------#

load("../RData/AllResps2017-08-09.RData")
load("../RData/Global_ANOVA_2017-08-09.RData")

#-------------------------------------------------------------------------------------------------#
#-----Generate plots looking at ANOVA hit calls and adjusted maxmMD
#-------------------------------------------------------------------------------------------------#

adj_mDist <- Mahalanobis_dat[,.(date_chnm_plate, maxD11P, adj_mDist = maxD11P - Scrit01)]

anova_sum <- OECD_ANOVA_filter_wide[,.(date_chnm_plate, casn, sum_anova = rowSums(OECD_ANOVA_filter_wide[,c(2:12), with = FALSE]))]

combined_dat <- merge(adj_mDist, anova_sum, by = "date_chnm_plate")

combined_dat <- combined_dat[order(sum_anova, decreasing = TRUE),]

#annotate OECD reference chemcials
casNOs <- c("1912-24-9", "17804-35-2", "94-26-8", "66575-29-9", "2212-67-1", "26027-38-3", "67747-09-5", "125-84-8", "4672-49-5",
            "112809-51-5", "13647-35-3", "80-05-7", "51-28-5", "17230-88-5", "117-81-7", "60-51-5", "60168-88-9", "98319-26-7",
            "13311-84-7", "446-72-0", "65277-42-1", "84371-65-3", "51-03-6", "1610-18-0", "52-01-7", "1330-78-5")


combined_dat[, reference := ifelse(casn %in% casNOs, 1, 0)]

ref_dat <- combined_dat[reference == 1]

ref_dat[, date_chnm_plate := factor(date_chnm_plate, levels = c("20130321_Prochloraz_Plate4", "20140402_Ketoconazole_Plate12",              
                                                                "20140409_Forskolin_Plate17", "20140903_Danazol_Plate21",                   
                                                                "20150610_Danazol_Plate40", "20170411_Fenarimol_Plate15",                 
                                                                "20130320_Atrazine_Plate14", "20130320_Prometon_Plate10",                
                                                                "20170411_Di(2-ethylhexyl) phthalate_Plate6", "20140326_Spironolactone_Plate10",            
                                                                "20140402_Genistein_Plate19", "20140326_Bisphenol A_Plate20",            
                                                                "20150610_Tricresyl phosphate_Plate16", "20170411_Letrozole_Plate4",              
                                                                "20140326_Bisphenol A_Plate6", "20140326_Piperonyl butoxide_Plate2",         
                                                                "20140409_Finasteride_Plate24", "20170411_Aminoglutethimide_Plate9",        
                                                                "20170411_Ethylene dimethanesulfonate_Plate3", "20130320_Molinate_Plate14",         
                                                                "20130321_Bisphenol A_Plate7", "20150610_Benomyl_Plate5",                 
                                                                "20150610_Nonoxynol_Plate17", "20140402_Mifepristone_Plate5",               
                                                                "20150610_Butylparaben_Plate12", "20150610_Dimethoate_Plate7"))]

#-------------------------------------------------------------------------------------------------#
#-----plot code
#-------------------------------------------------------------------------------------------------#

box_plot <- ggplot(data = combined_dat, aes(x = factor(sum_anova), y = maxD11P)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(data = combined_dat[adj_mDist >= 0 & reference == 0], width = 0.2, size = 2) +
  geom_jitter(data = combined_dat[adj_mDist < 0 & reference == 0], width = 0.2, shape = 1, size = 2) +
  geom_point(data = ref_dat, aes(x = factor(sum_anova), y = maxD11P, shape = factor(date_chnm_plate), color = factor(date_chnm_plate)), size = 4) +
  xlab("Steroid Hit Count") +
  ylab("maxmMd") +
  scale_y_log10() +
  scale_shape_manual(name = "OECD Reference Chemicals", labels = c("Prochloraz", "Ketoconazole", "Forskolin", "Danazol Rep 1", "Danazol Rep 2", "Fenarimol",
                                                                   "Atrazine", "Prometon", "Di(2-ethylhexyl phthalate", "Spironolactone", "Genistein", "Bisphenol A Rep 3",
                                                                   "Tricresyl phosphate", "Letrozole", "Bisphenol A Rep 2", "Piperonyl butoxide", "Finasteride",
                                                                   "Aminoglutethimide", "EDS", "Molinate", "Bisphenol A Rep 1", "Benomyl", "Nonoxynol 9", "Mifepristone", "Butylparaben", "Dimethoate"),
                     values = c(15:18, 8, 15:18, 8, 15:18, 8, 15:17, 6, 15:16, 17, 18, 8, 15, 8)) +
  scale_color_manual(name = "OECD Reference Chemicals", labels = c("Prochloraz", "Ketoconazole", "Forskolin", "Danazol Rep 1", "Danazol Rep 2", "Fenarimol",
                                                                   "Atrazine", "Prometon", "Di(2-ethylhexyl phthalate", "Spironolactone", "Genistein", "Bisphenol A Rep 3",
                                                                   "Tricresyl phosphate", "Letrozole", "Bisphenol A Rep 2", "Piperonyl butoxide", "Finasteride",
                                                                   "Aminoglutethimide", "EDS", "Molinate", "Bisphenol A Rep 1", "Benomyl", "Nonoxynol 9", "Mifepristone", "Butylparaben", "Dimethoate"),
                     values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#fb9a99", 
                                "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#fb9a99",
                                "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#fb9a99",
                                "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#fb9a99",
                                "#4daf4a", "#377eb8")) +
  annotation_logticks(sides = "l") +
  theme_few() +
  theme(axis.title = element_text(size = 16))

pdf(paste0("anovaHit_adjMD_boxplot_log10_", Sys.Date(), ".pdf"), width = 14, height = 8)
box_plot
dev.off()

