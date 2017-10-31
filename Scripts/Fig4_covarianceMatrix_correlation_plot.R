#-----------------------------------------------------------------------------------#
# Haggard et al. "High -throughput H295R steroidogenesis assay: utility as an
# alternative and a statistical approach to characterize effects on steroidogenesis"

# Figure 4 Correlation of covariance matrix

# This figure is to show the correlation of the covariance of the residaul matrix across block
#to better justify the use of Mahalanobis distance

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
#-----Corrplot of the 11x11 covariance matrix
#-------------------------------------------------------------------------------------------------#

load("../RData/CovTP0_2017-08-09.RData")

colnames(CovTP0) <- c("CORT", "CORTIC", "11DCORT", "ANDR", "DOC", "TESTO", "OHPROG", "OHPREG", "PROG",
                      "E1", "E2")

rownames(CovTP0) <- c("CORT", "CORTIC", "11DCORT", "ANDR", "DOC", "TESTO", "OHPROG", "OHPREG", "PROG",
                      "E1", "E2")

cor_cov <- cov2cor(CovTP0)

pdf(paste0("H295R_covariance_corrplot_", Sys.Date(), ".pdf"), width = 16, height = 12)
corrplot(corr = cor_cov, method = "color", tl.col = "black", tl.pos = "d", diag = FALSE,
         cl.align.text = "l", cl.cex = 1, tl.cex = 1, number.cex = 0.9, addgrid.col = "grey", cl.pos = "b",
         order = "hclust", addCoef.col = "black")
# corrplot(corr = cor_cov, order = "hclust", type = "lower", method = "color", tl.col = "black", tl.pos = "d",
#          cl.cex = 1, cl.align.text = "l", tl.cex = 1, number.cex = 0.9, addgrid.col = "grey")
# corrplot(corr = cor_cov, order = "hclust", type = "upper", method = "number",  diag = F, tl.pos = "n",
#          number.cex = 0.9, addgrid.col = "grey", cl.pos = "n", col = "black", add = TRUE)
dev.off()
