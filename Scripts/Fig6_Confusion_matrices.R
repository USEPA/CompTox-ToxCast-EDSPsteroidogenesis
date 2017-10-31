#-----------------------------------------------------------------------------------#
# Haggard et al. "High -throughput H295R steroidogenesis assay: utility as an
# alternative and a statistical approach to characterize effects on steroidogenesis"

# Figure 6 confusion matrix

# This figure is to compare the results for E2 and T for the OECD reference chemicals

#Original 09 June 2017
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
library(grid)
library(gridExtra)
library(tcpl)
tcplConf()
tcplConfList()

#NOTE: set working directory appropriately

#-----------------------------------------------------------------------------------#
# Figure 8 confusion matrix
#-----------------------------------------------------------------------------------#

# First set: includes all except equivocals from OECD
#-----------------------------------------------------------------------------------#
#Testosterone_up information
OECD<-factor(c(0,0,1,1))
HT<-factor(c(0,1,0,1))
Ytu=c(18,2,0,1)
Testosterone_up<-data.frame(OECD, HT, Ytu)

tu<-ggplot(data=Testosterone_up, mapping=aes(x=OECD, y=HT))+
  geom_tile(aes(fill=Ytu), color='white')+
  geom_text(aes(label=sprintf('%1.0f', Ytu)),fontface='bold', lineheight=2, vjust=1)+
  scale_fill_gradient(low='#e5f5e0', high='#31a354')+
  theme_bw()+
  theme(legend.position='none',
        legend.title=element_blank())+
  theme(axis.title.x = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.y = element_text(size=12),
        plot.title=element_text(lineheight=0.6, face='bold', hjust=0.5))+
  ggtitle('Testosterone up')

#Testosterone_dn information
Ytd=c(12,1,5,6)
Testosterone_dn<-data.frame(OECD, HT, Ytd)

td<-ggplot(data=Testosterone_dn, mapping=aes(x=OECD, y=HT))+
  geom_tile(aes(fill=Ytd), color='white')+
  geom_text(aes(label=sprintf('%1.0f', Ytd)),fontface='bold', lineheight=2, vjust=1)+
  scale_fill_gradient(low='#e5f5e0', high='#31a354')+
  theme_bw()+
  theme(legend.position='none',
        legend.title=element_blank())+
  theme(axis.title.x = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.y = element_text(size=12),
        plot.title=element_text(lineheight=0.6, face='bold', hjust=0.5))+
  ggtitle('Testosterone dn')

#Estradiol_up information
Yeu=c(11,2,2,6)
Estradiol_up<-data.frame(OECD, HT, Yeu)

eu<-ggplot(data=Estradiol_up, mapping=aes(x=OECD, y=HT))+
  geom_tile(aes(fill=Yeu), color='white')+
  geom_text(aes(label=sprintf('%1.0f', Yeu)),fontface='bold', lineheight=2, vjust=1)+
  scale_fill_gradient(low='#F0E442', high='#FF6600')+
  theme_bw()+
  theme(legend.position='none',
        legend.title=element_blank())+
  theme(axis.title.x = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.y = element_text(size=12),
        plot.title=element_text(lineheight=0.6, face='bold', hjust=0.5))+
  ggtitle('Estradiol up')

#Estradiol_dn information
Yed=c(17,1,1,4)
Estradiol_dn<-data.frame(OECD, HT, Yed)

ed<-ggplot(data=Estradiol_dn, mapping=aes(x=OECD, y=HT))+
  geom_tile(aes(fill=Yed), color='white')+
  geom_text(aes(label=sprintf('%1.0f', Yed)),fontface='bold', lineheight=2, vjust=1)+
  scale_fill_gradient(low='#F0E442', high='#FF6600')+
  theme_bw()+
  theme(legend.position='none',
        legend.title=element_blank())+
  theme(axis.title.x = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.y = element_text(size=12),
        plot.title=element_text(lineheight=0.6, face='bold', hjust=0.5))+
  ggtitle('Estradiol dn')


#Confusion matrix set that includes revised matrices (drop nonoxynol-9 from all and letrozole from T_dn)
#-----------------------------------------------------------------------------------#

Ytu_rev=c(17,2,0,1)
Testosterone_up_rev<-data.frame(OECD, HT, Ytu_rev)

tu_rev<-ggplot(data=Testosterone_up_rev, mapping=aes(x=OECD, y=HT))+
  geom_tile(aes(fill=Ytu_rev), color='white')+
  geom_text(aes(label=sprintf('%1.0f', Ytu_rev)),fontface='bold', lineheight=2, vjust=1)+
  scale_fill_gradient(low='#e5f5e0', high='#31a354')+
  theme_bw()+
  theme(legend.position='none',
        legend.title=element_blank())+
  theme(axis.title.x = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.y = element_text(size=12),
        plot.title=element_text(lineheight=0.6, face='bold', hjust=0.5))+
  ggtitle('Testosterone up revised')

#Testosterone_dn information
Ytd_rev=c(12,1,3,6)
Testosterone_dn_rev<-data.frame(OECD, HT, Ytd_rev)

td_rev<-ggplot(data=Testosterone_dn_rev, mapping=aes(x=OECD, y=HT))+
  geom_tile(aes(fill=Ytd_rev), color='white')+
  geom_text(aes(label=sprintf('%1.0f', Ytd_rev)),fontface='bold', lineheight=2, vjust=1)+
  scale_fill_gradient(low='#e5f5e0', high='#31a354')+
  theme_bw()+
  theme(legend.position='none',
        legend.title=element_blank())+
  theme(axis.title.x = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.y = element_text(size=12),
        plot.title=element_text(lineheight=0.6, face='bold', hjust=0.5))+
  ggtitle('Testosterone dn revised')

#Estradiol_up information
Yeu_rev=c(10,2,2,6)
Estradiol_up_rev<-data.frame(OECD, HT, Yeu_rev)

eu_rev<-ggplot(data=Estradiol_up_rev, mapping=aes(x=OECD, y=HT))+
  geom_tile(aes(fill=Yeu_rev), color='white')+
  geom_text(aes(label=sprintf('%1.0f', Yeu_rev)),fontface='bold', lineheight=2, vjust=1)+
  scale_fill_gradient(low='#F0E442', high='#FF6600')+
  theme_bw()+
  theme(legend.position='none',
        legend.title=element_blank())+
  theme(axis.title.x = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.y = element_text(size=12),
        plot.title=element_text(lineheight=0.6, face='bold', hjust=0.5))+
  ggtitle('Estradiol up revised')

#Estradiol_dn information
Yed_rev=c(17,0,1,4)
Estradiol_dn_rev<-data.frame(OECD, HT, Yed_rev)

ed_rev<-ggplot(data=Estradiol_dn_rev, mapping=aes(x=OECD, y=HT))+
  geom_tile(aes(fill=Yed_rev), color='white')+
  geom_text(aes(label=sprintf('%1.0f', Yed_rev)),fontface='bold', lineheight=2, vjust=1)+
  scale_fill_gradient(low='#F0E442', high='#FF6600')+
  theme_bw()+
  theme(legend.position='none',
        legend.title=element_blank())+
  theme(axis.title.x = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.y = element_text(size=12),
        plot.title=element_text(lineheight=0.6, face='bold', hjust=0.5))+
  ggtitle('Estradiol dn revised')


summary<-data.frame(Effect=c('Testosterone up',
                             'Testosterone dn',
                             'Estradiol up',
                             'Estradiol dn'),
                    Sensitivity=c(1,0.55, 0.75, 0.8),
                    Specificity=c(0.9, 0.92, 0.85, 0.94),
                    Accuracy=c(0.90, 0.75, 0.81, 0.91),
                    Sensitivity_Revised = c(1,0.67,0.75, 0.8),
                    Specificity_Revised = c(0.89, 0.92, 0.83, 1),
                    Accuracy_Revised = c(0.90, 0.82, 0.80, 0.95))

names(summary) <- c("Effect", "Sensitivity", "Specificity", "Accuracy", "Revised Sensitivity", "Revised Specificity", "Revised Accuracy")
# Set theme to allow for plotmath expressions
ttheme_minimal <- ttheme_minimal(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(summary, rows=NULL, theme=ttheme_minimal)
# Plot chart and table into one object

pdf(file="Draft_Confusion_Matrix_Figure_20June2017.pdf",
    height = 10,
    width = 9)

grid.arrange(tu,td,eu,ed, tbl, tu_rev,td_rev,eu_rev,ed_rev,
             ncol=2,
             as.table=TRUE,
             heights=c(2,2,2,2,2), 
             layout_matrix=rbind(c(1,2),
                                 c(3,4),
                                 c(5,5),
                                 c(6,7),
                                 c(8,9)))

dev.off()
