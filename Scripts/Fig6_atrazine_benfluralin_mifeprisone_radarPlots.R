#-----------------------------------------------------------------------------------#
# Haggard et al. "Using a high-throughput steroidogenesis assay to predict effects on 
# ...estrogen and androgen production or overall perturbation of the steroid biosynthesis pathway"

# Figure 6. Atrazine, benfluralin, and mifepristone radar plots

# This figure is to show the steroidogenesis profile and Mahalanobis distances for three example chemicals

#Original 08/14/2017

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
library(ggiraphExtra)
library(cowplot)

setwd("../Test Chemical Figures")

#-------------------------------------------------------------------------------------------------#
#-----Load in all of the data from the mahalanobis_script.R output RData file and the ANOVA tables
#-------------------------------------------------------------------------------------------------#

load("../RData/AllResps2017-08-09.RData")
load("../RData/Global_ANOVA_2017-08-09.RData")

#-------------------------------------------------------------------------------------------------#
#-----Add critical limits to Dists from Mahalanobis_dat
#-------------------------------------------------------------------------------------------------#

colnames(Dists)[4] <- "date_chnm_plate"

Dists[, date_chnm_plate := as.character(date_chnm_plate)]

for(x in unique(Dists$date_chnm_plate)){
  Dists[date_chnm_plate == x, c("Scrit05", "Scrit01") := .(Mahalanobis_dat[date_chnm_plate == x, Scrit05], Mahalanobis_dat[date_chnm_plate == x, Scrit01])]
}

#-------------------------------------------------------------------------------------------------#
#-----combine dat_mean and Dists
#-------------------------------------------------------------------------------------------------#

dat <- dat[, c(2, 15:19, 9:10, 6:8, 14)]
dat_wide <- dcast.data.table(dat, ... ~ steroid, value.var = "uM")
dat_wide[, c("D11P", "Scrit05", "Scrit01") := 0]

dat_wide_DMSO <- dat_wide[chnm == "DMSO",]
dat_wide_chems <- dat_wide[chnm != "DMSO",]

for(x in unique(dat_wide_chems[, date_chnm_plate])){
  for(y in unique(dat_wide_chems[date_chnm_plate == x, conc])){
    dat_wide_chems[date_chnm_plate == x & conc == y, c("D11P", "Scrit05", "Scrit01") := .(Dists[date_chnm_plate == x & Conc == y, D11P],
                                                                                          Dists[date_chnm_plate == x & Conc == y, Scrit05],
                                                                                          Dists[date_chnm_plate == x & Conc == y, Scrit01])]
  }
}

dat_combined_wide <- rbind(dat_wide_chems, dat_wide_DMSO)

#-------------------------------------------------------------------------------------------------#
#-----Convert to long format for plotting
#-------------------------------------------------------------------------------------------------#

dat_combined <- melt.data.table(dat_combined_wide, id.vars = colnames(dat_combined_wide)[c(1:10, 23:24)], value.name = "Value", variable.name = "Measurement")

#-------------------------------------------------------------------------------------------------#
#-------drop the NAs which show up due to NR values in the raw data
#-------------------------------------------------------------------------------------------------#

dat_combined <- dat_combined[!is.na(Value),]

#-------------------------------------------------------------------------------------------------#
#-----Fix the levels of the Measurement factors
#-------------------------------------------------------------------------------------------------#

dat_combined[Measurement == "D11P", Measurement := "Mahalanobis Distance",]

#-------------------------------------------------------------------------------------------------#
#-----Fix Measurement labels
#-------------------------------------------------------------------------------------------------#

dat_combined[Measurement == "OH-Progesterone", Measurement := "OHPROG"]
dat_combined[Measurement == "OH-Pregnenolone", Measurement := "OHPREG"]
dat_combined[Measurement == "Progesterone", Measurement := "PROG"]
dat_combined[Measurement == "Pregnenolone", Measurement := "PREG"]
dat_combined[Measurement == "Androstenedione", Measurement := "ANDR"]
dat_combined[Measurement == "11-deoxycortisol", Measurement := "11DCORT"]
dat_combined[Measurement == "Testosterone", Measurement := "TESTO"]
dat_combined[Measurement == "Corticosterone", Measurement := "CORTIC"]
dat_combined[Measurement == "Estrone", Measurement := "E1"]
dat_combined[Measurement == "Estradiol", Measurement := "E2"]
dat_combined[Measurement == "Cortisol", Measurement := "CORT"]

dat_combined[, Measurement := factor(Measurement, levels = c("OHPREG", "PROG", "OHPROG",
                                                             "DOC", "CORTIC", "11DCORT", "CORT",
                                                             "ANDR", "TESTO", "E1", "E2", "Mahalanobis Distance"))]

ANOVA_summary_out_long[steroid == "OH-Progesterone", steroid := "OHPROG"]
ANOVA_summary_out_long[steroid == "OH-Pregnenolone", steroid := "OHPREG"]
ANOVA_summary_out_long[steroid == "Progesterone", steroid := "PROG"]
ANOVA_summary_out_long[steroid == "Pregnenolone", steroid := "PREG"]
ANOVA_summary_out_long[steroid == "Androstenedione", steroid := "ANDR"]
ANOVA_summary_out_long[steroid == "11-deoxycortisol", steroid := "11DCORT"]
ANOVA_summary_out_long[steroid == "Testosterone", steroid := "TESTO"]
ANOVA_summary_out_long[steroid == "Corticosterone", steroid := "CORTIC"]
ANOVA_summary_out_long[steroid == "Estrone", steroid := "E1"]
ANOVA_summary_out_long[steroid == "Estradiol", steroid := "E2"]
ANOVA_summary_out_long[steroid == "Cortisol", steroid := "CORT"]

#-------------------------------------------------------------------------------------------------#
#-----Clean up ANOVA concentrations
#-------------------------------------------------------------------------------------------------#

ANOVA_summary_out_long[conc >= 99.9, conc := 100]
ANOVA_summary_out_long[str_detect(date_chnm_plate,"Triclosan") & conc == 0.004115226, conc := 0.004] 

ANOVA_summary_out_long[, conc := signif(conc, digits = 2)]

#-------------------------------------------------------------------------------------------------#
#-----Subset data for the chemicals of interest and add P-values
#-------------------------------------------------------------------------------------------------#

dat_combined[, Pvalue := 1]

chem_dat <- dat_combined[chnm == "Benfluralin" | chnm == "Atrazine" |
                           chnm == "Mifepristone" | chnm == "DMSO", ]

for(x in unique(chem_dat[chnm != "DMSO", date_chnm_plate])){
  for(z in unique(chem_dat[date_chnm_plate == x & Measurement != "Mahalanobis Distance", Measurement])){
    for(y in unique(chem_dat[date_chnm_plate == x & Measurement == z & chnm != "DMSO", conc])){
      chem_dat[date_chnm_plate == x & conc == y & Measurement == z, Pvalue := ANOVA_summary_out_long[date_chnm_plate == x & conc == y & steroid == z, PValue]]
    }
  }
}

chem_dat[, sig := factor(ifelse(Pvalue <= 0.05, 1, 0))]

#-------------------------------------------------------------------------------------------------#
#-----Generate radar plots for all test chems
#-------------------------------------------------------------------------------------------------#

dmso_block_average <- dat_combined[chnm == "DMSO", .(plate_dmso = mean(Value)), by = .(Measurement, date_chnm_plate)] #blocking by plate

dmso_block_average <- dmso_block_average[Measurement != "Mahalanobis Distance",]
dmso_block_average[, date_plate := str_replace(date_chnm_plate, "_DMSO_", "_")]

dat_radar <- copy(chem_dat)

dat_radar[, date_plate := paste(date, plate, sep = "_")]

for(x in unique(as.character(dat_radar[Measurement != "Mahalanobis Distance", Measurement]))){
  for(block in unique(dat_radar[chnm != "DMSO", date_plate])){
    for(y in unique(dat_radar[date_plate == block, date_chnm_plate]))
      dat_radar[Measurement == x & date_plate == block & date_chnm_plate == y, fold_change := Value/dmso_block_average[Measurement == x & date_plate == block, plate_dmso]]
  }
}

dat_radar[, log_fold_change := log2(fold_change)]


#individual chem plots
for(x in unique(chem_dat[chnm != "DMSO", date_chnm_plate])){
  pdf(paste0(x, "_radarPlot_mMD_", Sys.Date(), ".pdf"), width = 12, height = 8)
  p <- vector("list", 2)
  
  radar_dat <- dat_radar[date_chnm_plate == x & Measurement != "Mahalanobis Distance"]
  con <- paste(unique(radar_dat[, date]), "DMSO", unique(radar_dat[, plate]), sep = "_")
  radar_dat <- dat_radar[date_chnm_plate == x & Measurement != "Mahalanobis Distance" | date_chnm_plate == con & Measurement != "Mahalanobis Distance",]
  radar_dat[, conc := factor(conc)]
  
  lims <- c(min(radar_dat[,log_fold_change]) - 0.5, max(radar_dat[, log_fold_change]) + 0.5)
  
  radar_dat <- dcast.data.table(radar_dat[Measurement != "Mahalanobis Distance", c(1:6, 9:10, 13, 19)], ... ~ Measurement,
                                value.var = c("log_fold_change"), fun.aggregate = mean) #convert to wide format
  
  radar <- ggRadar(data = radar_dat[, -c(1:7)], aes(color = conc), rescale = FALSE, legend.position = "bottom", alpha = 0) +
    geom_hline(yintercept = c(log2(1.5), log2(1/1.5)), linetype = "dashed", color = "red") +
    coord_polar() +
    scale_color_brewer(palette = "YlGnBu") +
    ylim(lims) +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold")) +
    theme(axis.text = element_text(size = 10)) +
    theme(legend.text = element_text(size = 13)) +
    theme(legend.title = element_text(size = 13)) +
    labs(color = "Conc.", shape = "Conc.", fill = "Conc.") +
    theme(plot.margin=unit(c(0,0,0,0),"mm"))
  
  p[[1]] <- radar
  
  mMD_dat <- dat_combined[date_chnm_plate == x,]
  con <- paste(unique(mMD_dat[, date]), "DMSO", unique(mMD_dat[, plate]), sep = "_")
  mMD_dat <- dat_combined[date_chnm_plate == x & Measurement == "Mahalanobis Distance" | date_chnm_plate == con & Measurement == "Mahalanobis Distance",]
  
  mMD <- ggplot(data = mMD_dat[chnm != "DMSO",], aes(x = conc, y = Value)) +
    geom_point() +
    geom_hline(yintercept = unique(mMD_dat[chnm != "DMSO", Scrit01]), linetype = "dashed", color = "red") +
    #ggtitle(label = "Mahalanobis Distance") +
    scale_x_log10("Concentration (uM)") +
    ylab("Distance") +
    theme(plot.margin=unit(c(0,0,0,0),"mm")) +
    theme_few() +
    theme(axis.title = element_text(size = 13))
  
  p[[2]] <- mMD
  
  lay <- rbind(c(1, 1, NA), c(1, 1, 2), c(1, 1, 2), c(1, 1, NA))
  
  grid.arrange(grobs = p, layout_matrix = lay,
               top = textGrob(str_split_fixed(x, "_", 3)[,2], gp=gpar(fontsize=16)))
  dev.off()
}

