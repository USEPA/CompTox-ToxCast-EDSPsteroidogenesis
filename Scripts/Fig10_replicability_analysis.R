#-----------------------------------------------------------------------------------#
# Haggard et al. "High -throughput H295R steroidogenesis assay: utility as an
# alternative and a statistical approach to characterize effects on steroidogenesis"


#Figure 10.
#Replicability analysis

# This script with perform the replicability analysis for the internally replicated chemicals
# across blocks, and generate the appropriate Supp 12 figures and files

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

#-------------------------------------------------------------------------------------------------#
#-----Look at recall sensitivity and maxmMd of replicated chemicals
#-------------------------------------------------------------------------------------------------#

rep_chems <- unique(Mahalanobis_dat[duplicated(casn), casn])

#-------------------------------------------------------------------------------------------------#
#-----filter OECD_ANOVA_filter_wide by the replicated chemicals
#-------------------------------------------------------------------------------------------------#

rep_OECD_anova <- OECD_ANOVA_filter_wide[casn %in% rep_chems, -1]

rep_OECD_anova_recall <- vector("list", 107)
names(rep_OECD_anova_recall) <- unique(rep_OECD_anova[,chnm])

for(x in unique(rep_OECD_anova[,chnm])){
  rep_OECD_anova_recall[[x]] <- data.table(chnm = x, casn = unique(rep_OECD_anova[chnm == x, casn]), rep_OECD_anova[chnm == x,lapply(.SD, function(n) sum(n)/rep_OECD_anova[chnm == x, .N]), .SDcols = colnames(rep_OECD_anova)[c(-12, -13)]]) 
}

rep_OECD_anova_recall <- do.call(rbind, rep_OECD_anova_recall)

#any replicates that matched with no significant effects will have a zero value, change this to 1 to reflect a correct match
#for tripicate chemicals, a 2/3 agreement is considered a positive match
rep_OECD_anova_recall[rep_OECD_anova_recall == 0] <- 1
rep_OECD_anova_recall[rep_OECD_anova_recall < 0.6] <- 0

rep_OECD_anova_recall[, ave_chem_recall := (rowSums(.SD)/11)*100, .SDcols = colnames(rep_OECD_anova_recall)[c(-1, -2)]]

fwrite(rep_OECD_anova_recall, file = paste0("Supp12_OECD_anova_recall_", Sys.Date(), ".txt"), sep = "\t")

#-------------------------------------------------------------------------------------------------#
#-----plot maxmMD for each replicated chemical
#-------------------------------------------------------------------------------------------------#

rep_maxmMd <- Mahalanobis_dat[casn %in% rep_chems,]
rep_maxmMd[, adj_maxmMd := maxD11P - Scrit01]

rep_maxmMd <- rep_maxmMd[order(maxD11P)]
rep_maxmMd[adj_maxmMd > 0, hitc := 1]
rep_maxmMd[adj_maxmMd < 0, hitc := 0]

for(x in unique(rep_maxmMd[, chnm])){
  rep_maxmMd[chnm == x, CV := cv(rep_maxmMd[chnm == x, maxD11P])]
  rep_maxmMd[chnm == x, log_SD := sd(log(rep_maxmMd[chnm == x, maxD11P]))]
  rep_maxmMd[chnm == x, log_mean := mean((rep_maxmMd[chnm == x, maxD11P]))]
  rep_maxmMd[chnm == x, geo_mean := exp(mean(log(rep_maxmMd[chnm == x, maxD11P])))]
}

#-------------------------------------------------------------------------------------------------#
#-----fit data to a linear model with no intercept to determine the standard deviation
#-----of the residuals
#-------------------------------------------------------------------------------------------------#

maxmMd_lm <- summary(lm(log(maxD11P) ~ 0 + chnm, data = rep_maxmMd)) #0.33 is the residual st. error
#maxmMd_lm_corrected <- summary(lm(log(maxD11P) ~ 0 + chnm, data = rep_maxmMd[chnm != "DINP branched" & chnm != "1,2,4-Butanetriol"])) #0.2773 is the residaul standard error

#-------------------------------------------------------------------------------------------------#
#-----Generate Figure 10
#-------------------------------------------------------------------------------------------------#

maxmMd_plot <- ggplot(data = rep_maxmMd, aes(x = maxD11P, y = reorder(chnm, maxD11P), shape = factor(hitc))) +
  geom_point() +
  annotate("text", x = 1.1, y = 103, label = "Residual Standard Deviation = 0.33", size = 7, hjust = "left") +
  scale_x_log10(name = "maxmMd") +
  #scale_x_continuous(trans = "log") +
  annotation_logticks(sides = "b") +
  ylab("") +
  scale_shape_manual(values = c(1, 16)) +
  theme_few() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 20))


#Figure 
pdf(paste0("replicated_maxmMd_plot_log_", Sys.Date(), ".pdf"), width = 12, height = 15)
maxmMd_plot
dev.off()

tiff(paste0("replicated_maxmMd_plot_log_", Sys.Date(), ".tiff"), width = 12, height = 15, units = "in", res = 300)
maxmMd_plot
dev.off()

#-------------------------------------------------------------------------------------------------#
#-----Generate additional replicability figures
#-------------------------------------------------------------------------------------------------#
#collapse replicate samples to a single value for the standard deviation for the histogram
rep_maxmMd_diff <- rep_maxmMd[, .(maxmMd = max(maxD11P), minmMd = min(maxD11P), log_SD = min(log_SD)), by = .(chnm)]
#rep_maxmMd_diff <- rep_maxmMd[, .(maxmMd = max(log(maxD11P)), minmMd = min(log(maxD11P)), maxmMd_sd = sd(log(maxD11P)), maxmMd_diff = max((maxD11P)) - min((maxD11P))), by = .(casn, chnm)]


maxmMd_diff_histo <- ggplot(data = rep_maxmMd_diff, aes(x = log_SD)) +
  geom_histogram(binwidth = 0.03) +
  geom_vline(xintercept = 0.3308, linetype = "dashed") +
  #scale_x_log10() +
  theme_few() +
  xlab("Standard Deviation of ln(maxmMd)") +
  scale_x_continuous(limits = c(0, 1.5), breaks = c(0, 0.5, 1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 15), expand = c(0,0)) +
  ylab("Count") +
  theme(axis.title = element_text(size = 18))


# pdf(paste0("replicated_maxmMd_diff_distribution_", Sys.Date(), ".pdf"), width = 12, height = 8)
# maxmMd_diff_histo
# dev.off()

maxmMd_meanSD_plot <- ggplot(data = rep_maxmMd, aes(x = geo_mean, y = log_SD)) +
  geom_point() +
  geom_smooth() +
  geom_hline(yintercept = 0.3308, linetype = "dashed") +
  theme_few() +
  ylab("Standard Deviation of ln(maxmMd)") +
  xlab("maxmMd") +
  scale_x_continuous(trans = "log10", expand = c(0,0)) +
  #scale_y_continuous(trans = "log") +
  annotation_logticks(sides = "b") +
  theme(axis.title = element_text(size = 18))

# pdf(paste0("replicated_maxmMd_meanSD_", Sys.Date(), ".pdf"), width = 12, height = 8)
# maxmMd_meanSD_plot
# dev.off()

#supplemental plot
pdf(paste0("maxmMd_replicability_plots_supplement", Sys.Date(), ".pdf"), width = 12, height = 6)
plot_grid(maxmMd_diff_histo, maxmMd_meanSD_plot, labels = c("A.", "B."))
dev.off()

# pdf(paste0("maxmMd_replicability_plots_", Sys.Date(), ".pdf"), width = 20, height = 15)
# ggdraw() +
#   draw_plot(maxmMd_plot, 0,0, 0.6, 1) +
#   draw_plot(maxmMd_diff_histo, 0.61, 0.55, 0.38, 0.35) +
#   draw_plot(maxmMd_meanSD_plot, 0.61, 0.1, 0.38, 0.35) +
#   draw_plot_label(label = c("A", "B", "C"), x = c(0, 0.61, 0.61), y = c(1, 0.95, 0.5), size = 25)
# dev.off()