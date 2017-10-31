#-----------------------------------------------------------------------------------#
# Haggard et al. "High -throughput H295R steroidogenesis assay: utility as an
# alternative and a statistical approach to characterize effects on steroidogenesis"

# Code to run 4-parameter Hill models on mean Mahalanobis Distance values to calculate the BMD (conce where
# mMd goes above the critical limit).

# Figure 9.

rm(list = ls())

library(parallel)
library(ggplot2)
library(RColorBrewer)
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
#-----Functions used
#-------------------------------------------------------------------------------------------------#

logistic4par  <- function(lx, T, cc, d) {
    (1 + (cc - 1) / (1 + exp(-d * lx + T)))
}

BMD <- function(Z, T, cc, d) {
    lD <- (log((cc - 1) / (Z - 1) - 1) - T) / (- d)
    if (Z > cc | d < 0.01) Inf else if (is.finite(lD)) exp(lD) else lD
}

dofit <- function(Dsnm) {
    of <- function(parms) {
        T <- parms[1]
        cc <- exp(parms[2])
        d <- exp(parms[3])
        sum((DS[[Dsnm]]$y - logistic4par(DS[[Dsnm]]$lx, T, cc, d))^2)
    }
    Spcor <- cor(DS[[Dsnm]]$y, DS[[Dsnm]]$lx, method="spearman")
    Spcortest <- cor.test(DS[[Dsnm]]$y, DS[[Dsnm]]$lx, method="spearman")$p.value
    out <- optim(c(0, 1, 0), of, control=list(maxit=2000))
    Parms <- c(T = out$par[1], cc = exp(out$par[2]), d = exp(out$par[3]))
    list(Parms=Parms, BMD=unname(BMD(DS[[Dsnm]]$Scrit01, Parms[1], Parms[2], Parms[3])),
         MaxmMd = DS[[Dsnm]]$MaxmMd, Scrit = DS[[Dsnm]]$Scrit01,
         cor=Spcor, cor_pvalue=Spcortest, 
         convergence = out$convergence, Name = Dsnm)
}

doplots <- function(zz) {
    pdta1 <- data.frame(Concentration = exp(DS[[zz[["Name"]]]]$lx),
                        mMd = DS[[zz[["Name"]]]]$y)
    xmin <- if(is.finite(zz$BMD) & zz$BMD > 0) min(c(DS[[zz[["Name"]]]]$lx, log(zz$BMD))) else min(DS[[zz[["Name"]]]]$lx)
    pdta2 <- data.frame(x = exp(z <- seq(xmin,
                                         max(DS[[zz[["Name"]]]]$lx), length=300)),
                        y = logistic4par(z, zz$Parms["T"], zz$Parms["cc"], zz$Parms["d"]))
    p <- ggplot() +
        geom_point(data=pdta1, aes(x = Concentration, y=mMd)) +
        geom_line(data=pdta2, aes(x=x, y=y)) +
        annotate(geom="segment", x=min(pdta2$x), xend=zz$BMD, y=zz$Scrit, yend=zz$Scrit,
                 lty=2, color=brewer.pal(5, "YlGnBu")[3]) +
        annotate(geom="segment", x=zz$BMD, xend=zz$BMD, y=zz$Scrit, yend=0,
                 color=brewer.pal(5, "YlGnBu")[3]) +
        scale_x_log10() +
        ggtitle(zz[["Name"]])
    p
}



#-------------------------------------------------------------------------------------------------#
#-----Load output from mahalanobis_distance_calculation_and_Supp9.R 
#-------------------------------------------------------------------------------------------------#

load("./RData/AllResps_outliersRemoved2017-08-09.rdata")

## split out data 

ChemNames <- unique(Mahalanobis_dat$date_chnm_plate)

DS <- lapply(ChemNames, function(nm) {
    ds <- Dists[Dists$CA == nm,]
    list(lx = log(ds$Conc),
         y = ds$D11P,
         Scrit01 = Mahalanobis_dat[Mahalanobis_dat$date_chnm_plate == nm,"Scrit01"],
         MaxmMd = Mahalanobis_dat[Mahalanobis_dat$date_chnm_plate == nm,"maxD11P"]
         )
    })
names(DS) <- ChemNames

## Delete the objects we don't need (to reduce the amount of stuff copied over to the child
## processes
rm(CritLim, dat, dat_mean, Dists, Mahalanobis_dat, Residuals, Resps)
gc()

## -----------------------------------------------------------------
## Do the fits and construct plots

out <- mclapply(ChemNames, function(nm) {
    zz <- try(dofit(nm))
    p <- if (!is(zz, "try-error")) {
        doplots(zz)
    } else NA
        
    list(Name = nm,
         Fit = zz,
         Plot = p)
}, mc.cores=35)

Fits <- do.call(rbind,
                lapply(out,
                       function(x) {
                           cbind(data.frame(Name=x$Name,
                                            as.data.frame(matrix(c(x$Fit$Parms,BMD=x$Fit$BMD,
                                                                   MaxmMd=x$Fit$MaxmMd,
                                                                   Scrit=x$Fit$Scrit,
                                                                   cor=x$Fit$cor,
                                                                   cor_pvalue=x$Fit$cor_pvalue,
                                                                   convergence=x$Fit$convergence),
                                                                 nrow=1,
                                                                 dimnames=list(NULL, c(names(x$Fit$Parms),
                                                                                       "BMD","MaxmMd","Scrit","cor","cor_pvalue",
                                                                                       "convergence"))))))}))

save(Fits, file="./RData/mMdFits.RData")

sum(sapply(out, function(x) x$Fit$convergence != 0))

x <- sapply(out, function(x) x$Fit$MaxmMd)
y <- sapply(out, function(x) x$Fit$BMD)
y <- ifelse((yy <- ifelse(y <= 100, y, 1000)) < 1e-2, 1e-3,yy)
pdta <- data.frame(MaxmMd = x, BMD=y)

p <- ggplot(data=pdta, aes(x=MaxmMd, y=BMD)) +
    geom_point() + geom_smooth() +
    scale_y_log10("BMD") +
    scale_x_log10("MaxmMd")


pdf("MaxmMd_plots.pdf")
print(p)

indx <- which(y == 1e-3)
for (i in 1:length(indx)) print(out[[indx[i]]]$Plot)

indx <- which(y == 1000)
for (i in 1:length(indx)) print(out[[indx[i]]]$Plot)

indx <- which(y > 1e-3 & y < 1000)
for (i in 1:length(indx)) print(out[[indx[i]]]$Plot)

dev.off()

#-------------------------------------------------------------------------------------------------#
#-----Apply logic for pathway effect and trend for the Fits data frame
#-------------------------------------------------------------------------------------------------#

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

#-------------------------------------------------------------------------------------------------#
#-----Figure 9
#-------------------------------------------------------------------------------------------------#

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
  guides(color = FALSE, shape = FALSE) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))


p3

tiff("Fig9_maxmMd_by_BMD.tiff", width = 12, height = 8, units = "in", res = 300)
p3
dev.off()

pdf("Fig9_maxmMd_by_BMD.pdf", width = 12, height = 8)
p3
dev.off()

