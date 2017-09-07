#-----------------------------------------------------------------------------------#
# Haggard et al. "Using a high-throughput steroidogenesis assay to predict effects on 
# ...estrogen and androgen production or overall perturbation of the steroid biosynthesis pathway"

# Figure 3 prochloraz anova profile

# This figure is to show the steroidogenesis profile for prochloraz

#Original 08/11/2017

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
library(ggiraphExtra)
library(cowplot)

setwd("../Test Chemical Figures")

#-------------------------------------------------------------------------------------------------#
#-----Load in all of the data from the mahalanobis_script.R output RData file and the ANOVA tables
#-------------------------------------------------------------------------------------------------#

load("../RData/AllResps2017-08-09.RData")
load("../RData/Global_ANOVA_2017-08-09.RData")

#-----Add critical limits to Dists from Mahalanobis_dat
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

chem_dat <- dat_combined[chnm == "DMSO" | chnm == "Prochloraz", ]

for(x in unique(chem_dat[chnm != "DMSO", date_chnm_plate])){
  for(z in unique(chem_dat[date_chnm_plate == x & Measurement != "Mahalanobis Distance", Measurement])){
    for(y in unique(chem_dat[date_chnm_plate == x & Measurement == z & chnm != "DMSO", conc])){
      chem_dat[date_chnm_plate == x & conc == y & Measurement == z, Pvalue := ANOVA_summary_out_long[date_chnm_plate == x & conc == y & steroid == z, PValue]]
    }
  }
}

chem_dat[, sig := factor(ifelse(Pvalue <= 0.05, 1, 0))]

#-------------------------------------------------------------------------------------------------#
#-----Generate ANOVA plot for prochloraz
#-------------------------------------------------------------------------------------------------#

pdf(paste0("Prochloraz_concResp_", Sys.Date(), ".pdf"), width = 12, height = 8)
p <- vector("list", 11)
plotting_dat <- chem_dat[chnm == "Prochloraz"]
con <- paste(unique(plotting_dat[, date]), "DMSO", unique(plotting_dat[, plate]), sep = "_")
plotting_dat <- chem_dat[chnm == "Prochloraz" | date_chnm_plate == con,]
i <- 1
for(y in unique(levels(plotting_dat[Measurement != "Mahalanobis Distance", Measurement])[-12])){
  plot_dat <- plotting_dat[Measurement == y]
  plot_dat[, Fconc := factor(conc)]
  ft <- lm(Value ~ Fconc, data = plot_dat[Measurement != "Mahalanobis Distance"])
  pdta <- unique(plot_dat[,c("conc", "Fconc"), with = FALSE])
  pdta$mn <- predict(ft, newdata=pdta)
  
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
  
  i <- i +1
}
grid.arrange(grobs = p, #layout_matrix=layout,
             bottom=textGrob("Concentration (uM)", gp=gpar(fontsize=16)),
             left=textGrob("Measured Analyte (uM)", gp=gpar(fontsize=16), rot=90))

dev.off()