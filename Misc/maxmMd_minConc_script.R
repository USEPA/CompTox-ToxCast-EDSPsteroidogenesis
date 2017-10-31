#-----------------------------------------------------------------------------------#
# Haggard et al. "High -throughput H295R steroidogenesis assay: utility as an
# alternative and a statistical approach to characterize effects on steroidogenesis"

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
library(cowplot)


setwd("./Mahalanobis Distance")

#-------------------------------------------------------------------------------------------------#
#-----Load output from mahalanobis_distance_calculation_and_Supp9.R 
#-------------------------------------------------------------------------------------------------#

load("../RData/AllResps_outliersRemoved2017-08-09.RData")

#-------------------------------------------------------------------------------------------------#
#-----Generate table that contains the maxmMd and adjusted maxmMd for each conc
#-------------------------------------------------------------------------------------------------#

for(x in unique(Mahalanobis_dat[, date_chnm_plate])){
  Dists[CA == x, c("Scrit01", "maxmMd",  "adj_maxmMd") := .(Mahalanobis_dat[date_chnm_plate == x, Scrit01], Mahalanobis_dat[date_chnm_plate == x, maxD11P], Dists[CA == x, D11P] - Mahalanobis_dat[date_chnm_plate == x, Scrit01])]
}
Dists[, conc_index := as.character(conc_index)]

#-------------------------------------------------------------------------------------------------#
#-----Find minimum concentration where mMd is above the critical value for each test chemical
#-------------------------------------------------------------------------------------------------#

filtered_Dists <- Dists[adj_maxmMd > 0, .(Conc = min(Conc)), by = .(CA, maxmMd)]

for(x in unique(filtered_Dists[, CA])){
  filtered_Dists[CA == x, conc_index := Dists[CA == x & Conc == filtered_Dists[CA == x, Conc], conc_index]]
}

filtered_Dists[, conc_index := factor(conc_index, levels = c("CONC6", "CONC5", "CONC4", "CONC3", "CONC2", "CONC1"))]

#-------------------------------------------------------------------------------------------------#
#-----plot minimum concentration where mMd above the critical value by maxmMd
#-------------------------------------------------------------------------------------------------#

p1 <- ggplot(data = filtered_Dists, aes(x = conc_index, y = maxmMd)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_few() +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  xlab("Concentration Index") +
  ylab("maxmMd")
  
p2 <- ggplot(data = filtered_Dists, aes(x = maxmMd, y = Conc)) +
  geom_point() +
  geom_smooth() +
  theme_few() +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(sides = "lb") +
  xlab("maxmMd") +
  ylab("Concentration")

tiff("maxmMd_by_minConc_boxplot.tiff", width = 12, height = 8, units = "in", res = 300)
p1
dev.off()

pdf("maxmMd_by_minConc_boxplot.pdf", width = 12, height = 8)
p1
dev.off()

tiff("maxmMd_by_minConc.tiff", width = 12, height = 8, units = "in", res = 300)
p2
dev.off()

pdf("maxmMd_by_minConc.pdf", width = 12, height = 8)
p2
dev.off()

#-------------------------------------------------------------------------------------------------#
#-----Load output from the BMD analysis (Fits) of mMd data and apply logic for pathway effect and trend
#-------------------------------------------------------------------------------------------------#

load("../RData/mMdFits.RData")

Fits <- as.data.table(Fits)

Fits[MaxmMd >= Scrit & cor_pvalue <= 0.05, type := 1] #maxmMd over critical value and trend
Fits[MaxmMd >= Scrit & cor_pvalue >= 0.05, type := 2] #maxmMd over critical value but no trend
Fits[MaxmMd <= Scrit & cor_pvalue <= 0.05, type := 3] #maxmMd below critical value and trend
Fits[MaxmMd <= Scrit & cor_pvalue >= 0.05, type := 4] #maxmMd below critical value but no trend

#correct for small and infinite BMDs
Fits[BMD >= 150, BMD := 1e03]
Fits[BMD <= 1e-03, BMD := 1e-04]

#make table and table Grob
bmd_table <- data.table("No Trend" = c(48, 308), "Trend" = c(3, 407))
row.names(bmd_table) <- c("maxmMd < Critical Value", "maxmMd > Critical Value")

cols <- matrix(c("#92c5de", "#f4a582", "#313695", "#a50026"), nrow(bmd_table), ncol(bmd_table))
tt <- ttheme_minimal(
  core=list(fg_params = list(col = cols),
            bg_params = list(col=NA)),
  rowhead=list(bg_params = list(col=NA), fg_params = list(fontface = "bold")),
  colhead=list(bg_params = list(col=NA)))

t <- tableGrob(bmd_table, theme = tt)


p3 <- ggplot(data = Fits, aes(x = MaxmMd, y = BMD)) +
  geom_smooth(color = "black", alpha = 0.25) + #change acccordingly
  geom_point(data = Fits[type != 3], aes(color = factor(type), shape = factor(type)), size = 3, alpha = 0.5) +
  geom_point(data = Fits[type == 3], aes(color = factor(type), shape = factor(type)), size = 3, alpha = 0.5) +
  theme_few() +
  annotation_custom(grob = t, xmin = log10(15), xmax = log10(55), ymin = 4.25, ymax = 5.25) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(sides = "l", short = unit(0, "mm"), mid = unit(0, "mm")) +
  annotation_logticks(sides = "b") +
  scale_color_manual(values = c("#a50026", "#f4a582", "#313695", "#92c5de")) +
  scale_shape_manual(values = c(16, 17, 18, 15)) +
  xlab("maxmMd") +
  ylab("BMD (uM)") +
  guides(color = FALSE, shape = FALSE)
 

p3

tiff("maxmMd_by_BMD.tiff", width = 12, height = 8, units = "in", res = 300)
p3
dev.off()

pdf("maxmMd_by_BMD.pdf", width = 12, height = 8)
p3
dev.off()






