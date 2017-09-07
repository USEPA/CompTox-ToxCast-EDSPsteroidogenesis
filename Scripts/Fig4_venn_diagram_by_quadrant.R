#-----------------------------------------------------------------------------------#
# Haggard et al. "Using a high-throughput steroidogenesis assay to predict effects on 
# ...estrogen and androgen production or overall perturbation of the steroid biosynthesis pathway"

# Figure 4 Venn Diagram for all of the ANOVA results

#Original 07/10/2017
#Updated  07/10/2017

#-----------------------------------------------------------------------------------#
#loading libraries and setting directories
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
library(VennDiagram)


setwd('../ANOVA/Global ANOVA')

load('../../RData/Global_ANOVA_2017-08-09.RData')

#-----------------------------------------------------------------------------------#
#determine activity by quadrant
#-----------------------------------------------------------------------------------#

OECD_ANOVA_filter_wide[,ProgestagenQ :=sum(`OH-Pregnenolone`,
                                           Progesterone,
                                           `OH-Progesterone`), by=list(casn)]

OECD_ANOVA_filter_wide[,GlucocorticoidQ :=sum(DOC,
                                              Corticosterone,
                                              `11-deoxycortisol`,
                                              Cortisol), by=list(casn)]

OECD_ANOVA_filter_wide[,AndrogenQ :=sum(Androstenedione,
                                        Testosterone), by=list(casn)]

OECD_ANOVA_filter_wide[,EstrogenQ :=sum(Estrone,
                                        Estradiol), by=list(casn)]

#-----------------------------------------------------------------------------------#
#binarize for counting
#-----------------------------------------------------------------------------------#

OECD_ANOVA_filter_wide[ProgestagenQ > 0,ProgQbinary := 1]
OECD_ANOVA_filter_wide[GlucocorticoidQ > 0,GlucQbinary := 1]
OECD_ANOVA_filter_wide[AndrogenQ > 0,AndrQbinary := 1]
OECD_ANOVA_filter_wide[EstrogenQ > 0,EstrQbinary := 1]

OECD_ANOVA_filter_wide[ProgestagenQ == 0,ProgQbinary := 0]
OECD_ANOVA_filter_wide[GlucocorticoidQ == 0,GlucQbinary := 0]
OECD_ANOVA_filter_wide[AndrogenQ == 0,AndrQbinary := 0]
OECD_ANOVA_filter_wide[EstrogenQ == 0,EstrQbinary := 0]

#-----------------------------------------------------------------------------------#
#remove duplicated chemicals since we are looking at unique chemicals not samples
#-----------------------------------------------------------------------------------#

OECD_ANOVA_filter_wide <- OECD_ANOVA_filter_wide[!duplicated(OECD_ANOVA_filter_wide[,casn])]

#look at OECD ref chems separately
OECD_ANOVA_ref <-OECD_ANOVA_filter_wide[casn %in% casns]
colnames(OECD_ANOVA_ref)
is.data.table(OECD_ANOVA_ref)
OECD_ANOVA_ref[,Steroid_Hit_Count := sum(`OH-Pregnenolone`,
                                         Progesterone,
                                         `OH-Progesterone`, 
                                         DOC,
                                         Corticosterone,
                                         `11-deoxycortisol`,
                                         Cortisol,
                                         Androstenedione,
                                         Testosterone,
                                         Estrone,
                                         Estradiol), by=list(chnm)]
fwrite(OECD_ANOVA_ref, file='OECD_Ref_Chem_Hit&Quadrant_Count_11JUly2017.csv')
#-----------------------------------------------------------------------------------#
#sum quadrant hits for each chemical for reference
#-----------------------------------------------------------------------------------#

venn.count<-data.table(count(count(OECD_ANOVA_filter_wide, c("ProgQbinary",
                                      "GlucQbinary",
                                      "AndrQbinary",
                                      "EstrQbinary"))))
sum(venn.count$freq)
#653 casn samples

#-----------------------------------------------------------------------------------#
#try using venn.diagram() which automatically computes areas
#-----------------------------------------------------------------------------------#

area1<-OECD_ANOVA_filter_wide[ProgQbinary==1, chnm]

area2<-OECD_ANOVA_filter_wide[GlucQbinary==1, chnm]

area3<-OECD_ANOVA_filter_wide[AndrQbinary==1, chnm]

area4<-OECD_ANOVA_filter_wide[EstrQbinary==1, chnm]

venn.plot <- venn.diagram(x = list(Progestagen = area1, Glucocorticoid = area2, Androgen = area3, Estrogen = area4),
                          filename = NULL,
                          #fill = c("green", "yellow", "blue", "red"),
                          lty = "dashed",
                          cex = 2,
                          cat.cex = 2
                          #cat.col = c("green", "yellow", "blue", "red")
                          )

#-----------------------------------------------------------------------------------#
#Save plots
#-----------------------------------------------------------------------------------#

tiff(filename = "Quad_Venn_diagram_nocolor.tiff", compression = "lzw", width = 1000, height = 850)
grid.draw(venn.plot)
dev.off()

pdf("Quad_Venn_diagram_nocolor.pdf", width = 14, height = 10)
grid.draw(venn.plot)
dev.off()

