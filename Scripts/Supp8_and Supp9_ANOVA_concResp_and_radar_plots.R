#-----------------------------------------------------------------------------------#
# Haggard et al. "High -throughput H295R steroidogenesis assay: utility as an
# alternative and a statistical approach to characterize effects on steroidogenesis"

# Supp 8 ANOVA-annotated concentration response and Mahalanobis distance plots

# Supp 9 Radar and mMd plots

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

dat <- dat[, c(2,3, 15:19, 9:10, 6:8, 14)]
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
#-----Housekeeping: look at the number of chemicals with missing concs due to viability or NR
#-------------------------------------------------------------------------------------------------#

conc_counter <- dat_combined_wide[chnm != "DMSO", .(date_chnm_plate = date_chnm_plate, N = length(unique(conc))), by = "date_chnm_plate"]
conc_counter[N == 0, .N]#0
conc_counter[N == 1, .N]#0
conc_counter[N == 2, .N]#4 
conc_counter[N == 3, .N]#5
conc_counter[N == 4, .N]#6
conc_counter[N == 5, .N]#36
conc_counter[N == 6, .N]#715

#-------------------------------------------------------------------------------------------------#
#-----Convert to long format for plotting
#-------------------------------------------------------------------------------------------------#

dat_combined <- melt.data.table(dat_combined_wide, id.vars = colnames(dat_combined_wide)[c(1:11, 24:25)], value.name = "Value", variable.name = "Measurement")

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
dat_combined[Measurement == "Estrone", Measurement := "ESTRONE"]
dat_combined[Measurement == "Estradiol", Measurement := "ESTRADIOL"]
dat_combined[Measurement == "Cortisol", Measurement := "CORTISOL"]

dat_combined[, Measurement := factor(Measurement, levels = c("OHPREG", "PROG", "OHPROG",
                                                             "DOC", "CORTIC", "11DCORT", "CORTISOL",
                                                             "ANDR", "TESTO", "ESTRONE", "ESTRADIOL", "Mahalanobis Distance"))]

ANOVA_summary_out_long[steroid == "OH-Progesterone", steroid := "OHPROG"]
ANOVA_summary_out_long[steroid == "OH-Pregnenolone", steroid := "OHPREG"]
ANOVA_summary_out_long[steroid == "Progesterone", steroid := "PROG"]
ANOVA_summary_out_long[steroid == "Pregnenolone", steroid := "PREG"]
ANOVA_summary_out_long[steroid == "Androstenedione", steroid := "ANDR"]
ANOVA_summary_out_long[steroid == "11-deoxycortisol", steroid := "11DCORT"]
ANOVA_summary_out_long[steroid == "Testosterone", steroid := "TESTO"]
ANOVA_summary_out_long[steroid == "Corticosterone", steroid := "CORTIC"]
ANOVA_summary_out_long[steroid == "Estrone", steroid := "ESTRONE"]
ANOVA_summary_out_long[steroid == "Estradiol", steroid := "ESTRADIOL"]
ANOVA_summary_out_long[steroid == "Cortisol", steroid := "CORTISOL"]

#-------------------------------------------------------------------------------------------------#
#-----Clean up ANOVA table
#-------------------------------------------------------------------------------------------------#

ANOVA_summary_out_long[, date_chnm_plate := str_replace(date_chnm_plate, "PharmaGSID_", "PharmaGSID")]
ANOVA_summary_out_long[conc >= 99.9, conc := 100]
ANOVA_summary_out_long[str_detect(date_chnm_plate,"Triclosan") & conc == 0.004115226, conc := 0.004] 

ANOVA_summary_out_long[, conc := signif(conc, digits = 2)]

#-------------------------------------------------------------------------------------------------#
#-----Add p-value information to dat0_combined to annotate significant effects in plots
#-------------------------------------------------------------------------------------------------#

ANOVA_summary_out_long[is.na(PValue), PValue := 2]

dat_combined[, Pvalue := 1]


for(x in unique(dat_combined[chnm != "DMSO", date_chnm_plate])){
  for(z in unique(dat_combined[date_chnm_plate == x & Measurement != "Mahalanobis Distance", Measurement])){
    for(y in unique(dat_combined[date_chnm_plate == x & Measurement == z & chnm != "DMSO", conc])){
      dat_combined[date_chnm_plate == x & conc == y & Measurement == z, Pvalue := ANOVA_summary_out_long[date_chnm_plate == x & conc == y & steroid == z, PValue]]
    }
  }
}

dat_combined <- dat_combined[Pvalue != 2,]

dat_combined[, sig := factor(ifelse(Pvalue <= 0.05, 1, 0))]
dat_combined[,conc_index := factor(conc_index, levels = c("CONC6", "CONC5", "CONC4", "CONC3", "CONC2", "CONC1"))]

#-------------------------------------------------------------------------------------------------#
#-----Generate Supp 8 file (ANOVA plots with mMD plot as well)
#-------------------------------------------------------------------------------------------------#

pdf(paste0("Supp7_OECD_global_concResp_mMD_", Sys.Date(), ".pdf"), width = 12, height = 8)
for(x in unique(dat_combined[chnm != "DMSO", date_chnm_plate])){
  p <- vector("list", 12)
  plotting_dat <- dat_combined[date_chnm_plate == x]
  con <- paste(unique(plotting_dat[, date]), "DMSO", unique(plotting_dat[, plate]), sep = "_")
  plotting_dat <- dat_combined[date_chnm_plate == x | date_chnm_plate == con,]
  i <- 1
  for(y in unique(levels(plotting_dat[, Measurement]))){
    plot_dat <- plotting_dat[Measurement == y]
    plot_dat[, Fconc := factor(conc)]
    if(unique(plot_dat[, Measurement]) != "Mahalanobis Distance"){
      ft <- lm(Value ~ Fconc, data = plot_dat[Measurement != "Mahalanobis Distance"])
      pdta <- unique(plot_dat[,c("conc", "Fconc"), with = FALSE])
      pdta$mn <- predict(ft, newdata=pdta)
    }
    if(unique(plot_dat[, Measurement]) != "Mahalanobis Distance"){
      conc_resp_plot <- ggplot(data = plot_dat, aes(x = conc, y = (Value))) +
        geom_point(aes(color = sig)) +
        geom_point(data = pdta, aes(x = conc, y = (mn)), shape = 3) +
        geom_hline(yintercept = (1.5*mean(plot_dat[conc == 0, Value])), linetype = "dashed") +
        geom_hline(yintercept = (1/1.5*mean(plot_dat[conc == 0, Value])), linetype = "dashed") +
        ggtitle(label = y) +
        scale_x_log10("") +
        scale_color_manual(values = c("black", "red"), guide = FALSE) +
        ylab("") +
        theme_few()
      p[[i]] <- conc_resp_plot
    } else {
      conc_resp_plot <- ggplot(data = plot_dat[chnm != "DMSO",], aes(x = conc, y = Value)) +
        geom_point() +
        geom_hline(yintercept = unique(plot_dat[chnm != "DMSO", Scrit01]), linetype = "dashed", color = "red") +
        #geom_hline(yintercept = unique(plot_dat[chnm != "DMSO", Scrit05]), linetype = "dashed", color = "blue") +
        ggtitle(label = y) +
        scale_x_log10("") +
        ylab("Distance") +
        theme_few()
      p[[i]] <- conc_resp_plot
    }
    i <- i +1
  }
  grid.arrange(grobs = p, #layout_matrix=layout,
               top = textGrob(x, gp=gpar(fontsize=16)),
               bottom=textGrob("Concentration (uM)", gp=gpar(fontsize=16)),
               left=textGrob("Measured Analyte (uM)", gp=gpar(fontsize=16), rot=90))
  #break
}
dev.off()

#-------------------------------------------------------------------------------------------------#
#-----Generate fold-change of H295R response data relative to plate-wise DMSO samples
#-------------------------------------------------------------------------------------------------#

dmso_block_average <- dat_combined[chnm == "DMSO", .(plate_dmso = mean(Value)), by = .(Measurement, date_chnm_plate)]

dmso_block_average <- dmso_block_average[Measurement != "Mahalanobis Distance",]
dmso_block_average[, date_plate := str_replace(date_chnm_plate, "_DMSO_", "_")]

dat_radar <- copy(dat_combined)

dat_radar[, date_plate := paste(date, plate, sep = "_")]

for(x in unique(as.character(dat_radar[Measurement != "Mahalanobis Distance", Measurement]))){
  for(block in unique(dat_radar[chnm != "DMSO", date_plate])){
    for(y in unique(dat_radar[date_plate == block, date_chnm_plate]))
      dat_radar[Measurement == x & date_plate == block & date_chnm_plate == y, fold_change := Value/dmso_block_average[Measurement == x & date_plate == block, plate_dmso]]
  }
}

dat_radar[, log_fold_change := log2(fold_change)]
#dat_radar[, log_fold_change := log10(fold_change)]

dat_radar_wide <- dcast.data.table(dat_radar[Measurement != "Mahalanobis Distance", c(1:7, 10:11, 14, 20)], ... ~ Measurement,
                                   value.var = c("log_fold_change"), fun.aggregate = mean) #convert to wide format

#-------------------------------------------------------------------------------------------------#
#-----Generate Supp 9 file (Radar plots and maxmMd plot)
#-------------------------------------------------------------------------------------------------#
pdf(paste0("Supp8_OECD_GLOBAL_radar_mDist_plots_", Sys.Date(), ".pdf"), width = 16, height = 10)
for(x in unique(dat_combined[chnm != "DMSO", date_chnm_plate])){
  p <- vector("list", 2)
  
  radar_dat <- dat_radar[date_chnm_plate == x & Measurement != "Mahalanobis Distance"]
  con <- paste(unique(radar_dat[, date]), "DMSO", unique(radar_dat[, plate]), sep = "_")
  radar_dat <- dat_radar[date_chnm_plate == x & Measurement != "Mahalanobis Distance" | date_chnm_plate == con & Measurement != "Mahalanobis Distance",]
  radar_dat[, conc := factor(conc)]
  
  lims <- c(min(radar_dat[,log_fold_change]) - 0.5, max(radar_dat[, log_fold_change]) + 0.5)
  
  radar_dat <- dcast.data.table(radar_dat[Measurement != "Mahalanobis Distance", c(1:7, 10:11, 14, 20)], ... ~ Measurement,
                                value.var = c("log_fold_change"), fun.aggregate = mean) #convert to wide format
  
  radar <- ggRadar(data = radar_dat[, -c(1:7)], aes(col = conc), rescale = FALSE, legend = "right", alpha = 0) +
    geom_hline(yintercept = c(log2(1.5), log2(1/1.5)), linetype = "dashed", color = "red") +
    coord_polar() +
    scale_color_brewer(palette = "YlGnBu") +
    ylim(lims) +
    theme_minimal()
  
  p[[1]] <- radar
  
  mMD_dat <- dat_combined[date_chnm_plate == x,]
  con <- paste(unique(mMD_dat[, date]), "DMSO", unique(mMD_dat[, plate]), sep = "_")
  mMD_dat <- dat_combined[date_chnm_plate == x & Measurement == "Mahalanobis Distance" | date_chnm_plate == con & Measurement == "Mahalanobis Distance",]
  
  mMD <- ggplot(data = mMD_dat[chnm != "DMSO",], aes(x = conc, y = Value)) +
    geom_point() +
    geom_hline(yintercept = unique(mMD_dat[chnm != "DMSO", Scrit01]), linetype = "dashed", color = "red") +
    ggtitle(label = "Mahalanobis Distance") +
    scale_x_log10("Concentration (uM)") +
    ylab("Distance") +
    theme_few()
  
  p[[2]] <- mMD
  
  lay <- rbind(c(1, 1, NA), c(1, 1, 2), c(1, 1, 2), c(1, 1, NA))
  
  grid.arrange(grobs = p, layout_matrix = lay,
               top = textGrob(x, gp=gpar(fontsize=16)))
}
dev.off()
