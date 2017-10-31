#-----------------------------------------------------------------------------------#
# Haggard et al. "High -throughput H295R steroidogenesis assay: utility as an
# alternative and a statistical approach to characterize effects on steroidogenesis"

# Supplemental File 9 Distribution of adjusted maxmMd for all 766 H295R samples

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
#-----Generate histogram to show distribution of chemicals above the critical limit
#-------------------------------------------------------------------------------------------------#

mDist_diff <- copy(Mahalanobis_dat)

mDist_diff[, dist_diff := maxD11P - Scrit01]

pdf(paste0("mMD_histogram_", Sys.Date(), ".pdf"), width = 8, height = 5)
histo_plot <- ggplot(mDist_diff, aes(x = dist_diff)) +
  geom_histogram(binwidth = 1, fill = "white", color = "black") +
  geom_vline(xintercept = median(mDist_diff[,dist_diff]), linetype = "dashed", color = "red") +
  ylab("Frequency") +
  xlab("Adjusted maxmMd (unitless)") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 135)) +
  scale_x_continuous(expand = c(0,0), limits = c(-2, 55)) +
  theme_few()

print(histo_plot)
dev.off()

tiff(paste0("mMD_histogram_", Sys.Date(), ".tiff"), width = 8, height = 5, units = "in", res = 300)
histo_plot <- ggplot(mDist_diff, aes(x = dist_diff)) +
  geom_histogram(binwidth = 1, fill = "white", color = "black") +
  geom_vline(xintercept = median(mDist_diff[,dist_diff]), linetype = "dashed", color = "red") +
  ylab("Frequency") +
  xlab("Adjusted maxmMd (unitless)") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 135)) +
  scale_x_continuous(expand = c(0,0), limits = c(-2, 55)) +
  theme_few()

print(histo_plot)
dev.off()