#-----------------------------------------------------------------------------------#
# Haggard et al. "Using a high-throughput steroidogenesis assay to predict effects on 
# ...estrogen and androgen production or overall perturbation of the steroid biosynthesis pathway"


#This script will compute the mahalanobis distance and generate all of the necessary plots
#based on the manova residual matrix to determine the variance and covariance matrices to use
#in the mahalanobis distance calculation for all of the H295R data

#This script also generates the Supp. 9 table

#-----------------NOTES-----------------#
##The loaded manova_output.RData file contains three R data objects:
##(1) dat: the manova-corrected H295R master data file (see the manova_script.R to see how the original H295R master data file was altered)##
##(2) Pvalues: a data.table having all of the calculated p-values for the manova by block
##(3) Models: a list of the model output from the manova, which consists of a matrix of the residuals, which will be used in the covariance matrix estiamtes
#---------------------------------------#

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

setwd("../Mahalanobis Distance")

#-------------------------------------------------------------------------------------------------#
#-----Functions used in the script
#-------------------------------------------------------------------------------------------------#

#dated_file_name: takes a name for a file, and appends
#a date to it: input myfile.txt -> output myfile_2017-02-17.txt
dated_file_name <- function(fn) {
  stamp <- paste0(format(Sys.time(), "_%Y-%m-%d"),"-",
                  round(as.double(difftime(Sys.time(), strptime("00:00:00", format="%H:%M:%S"), units = "secs"))*100)
  )
  prts <- strsplit(fn, "\\.")[[1]]
  extens <- prts[length(prts)]
  front <- if (length(prts) > 2) paste(prts[1:(length(prts)-1)],collapse=".") else prts[1]
  paste0(front,stamp,".",extens)
}

#Implements algorithm in:
#Nakamura 2005,
#author = {Nakamura, Tomohiro and Imada, Tsunehisa},
#file = {:C$\backslash$:/Users/wsetzer/Downloads/18{\_}21.pdf:pdf},
#journal = {Journal of the Japanese Society of Computational Statistics},
#keywords = {Analysis of variance and covariance,Paired and multiple comparisons,distribution of statistics},
#number = {1},
#pages = {21--32},
#title = {{Multiple Comparison Procedure of Dunnett's Type for Multivariate Normal Means}},
#volume = {18},
#year = {2005}

PT2maxmv  <-  function(c, K, p) {
  fn <- function(v, c, p, K) {
    exp((K - 1) * pchisq(c, df=p, ncp=v, log.p=TRUE)) *
      dchisq(v, p)
  }
  tmp <- integrate(fn, 0, Inf, c=c, p=p, K=K, rel.tol=1e-12)
  tmp$value
}

critval <- function(alpha, K, p) {
  if (K < 2 || p < 1) {
    list(root=NA)
  } else {
    ##    Find c in PT2maxmv such that PT2maxmv == 1 - alpha
    obfun <- function(c) {
      PT2maxmv(c, K=K, p=p) - (1 - alpha)
    }
    ## When c is small, obfun < 0, increases with c.
    ## Start with c=1
    cc <- 1
    Z <- obfun(cc)
    Niter <- 0
    while ((Z <- obfun(cc)) < 0 && ((Niter <- Niter + 1) < 100)) cc <- cc * 2
    if (Z > 0) {
      out <- uniroot(obfun, interval=c(cc/2,cc))
      out
    } else {
      list(root=NA)
    }
  }
}

#Empirical cdfs of the distribution of the various distance measures between the lowest
#pair of concentrations. Here is a little function to do that
getDist <- function(dmeasure){
  D <- sapply(unique(Dists$CA), function(ca) {
    tdta <- Dists[CA == ca,]
    tdta[which.min(Conc), D11P]
  })
  D
}

#-------------------------------------------------------------------------------------------------#
#-----Load output from manova_script.R which has three R objects: dat, Pvalues, and Models
#-------------------------------------------------------------------------------------------------#

load("../RData/manova_output2017-08-09.RData")

#-------------------------------------------------------------------------------------------------#
#-----Clean up dat
#-------------------------------------------------------------------------------------------------#

dat[spid == "DMSO", chnm := "DMSO"]

dat[spid == "DMSO", date_chnm_plate := str_replace(date_chnm_plate, "NA", "DMSO")]

dat <- dat[steroid != "DHEA" & steroid != "Pregnenolone"]

for(x in unique(dat[, date_chnm_plate])){
  for(y in unique(dat[,steroid])){
    for(z in unique(dat[date_chnm_plate == x & steroid == y, conc])){
      dat[date_chnm_plate == x & steroid == y & conc == z, N := dat[date_chnm_plate == x & steroid == y & conc == z, .N]]
    }
  }
}

#-------------------------------------------------------------------------#
#-----OPTIONAL: Filter the samples with residual outliers-----------------#
#-------------------------------------------------------------------------#
# dat[, conc_date_chnm_plate := paste(conc, date_chnm_plate, sep = "_")]
# 
# dat <- dat[!c(conc_date_chnm_plate %in% Residuals[maxSD == 1, conc_date_chnm_plate]),] #filter based on Residuals with the maxSD filter == 1
# dat[, conc_date_chnm_plate := NULL]

#-------------------------------------------------------------------------------------------------#
#-----Average the uM values between replicates across test chemicals
#-------------------------------------------------------------------------------------------------#

dat_mean <- dat[, .(mean_log_uM = mean(log_uM)), by = .(date_chnm_plate, steroid, chnm, conc, conc_index, wllt, date, N)]

CA <- unique(dat_mean$date_chnm_plate)
CA <- CA[!str_detect(CA, "DMSO")]

hormones_sub <- colnames(Pvalues)[5:17]
hormones_sub <- hormones_sub[c(-8, -11)]

#-------------------------------------------------------------------------------------------------#
#-----Generate covariance matrix based on manova residuals stored in Models
#-------------------------------------------------------------------------------------------------#

Blocks <- sort(names(Models))
Covs <- vector("list", length(Blocks))
names(Covs) <- Blocks
for (blk in Blocks) {
  Covs[[blk]] <- estVar(Models[[blk]]$fit)
}

#-------------------------------------------------------------------------------------------------#
#-----Generate pooled covariance matrix
#-------------------------------------------------------------------------------------------------#

nms <- unique(unlist(lapply(Covs, colnames)))
CovT <- array(NA, dim=c(length(nms), length(nms), length(Blocks)), dimnames=list(nms,nms,Blocks))
for (blk in Blocks) CovT[rownames(Covs[[blk]]),colnames(Covs[[blk]]),blk] <- Covs[[blk]]
CovTP <- apply(CovT,c(1,2), mean, na.rm=TRUE)
#CI11P <- solve(CovTP)
#VI11P <- 1/diag(CovTP)
#VIP <- 1/diag(CovTP)

CovTP0 <- CovTP
CovTP0[is.na(CovTP0)]  <- 0
save(CovTP0, file = paste0("../RData/CovTP0_", Sys.Date(), ".RData"))
#-------------------------------------------------------------------------------------------------#
#-----Let's try to calculate the Mahalanobis distance 
#-------------------------------------------------------------------------------------------------#

#Oxadiazon.Plate.4_1 and Thiabendazole.Plate.4_100 are missing Estradiol and Estrone values
out <- vector("list", length=length(CA))
names(out) <- CA
SSnms <- paste("N",nms,sep="_")
for (ca in CA) {
  tdta <- dat_mean[date_chnm_plate == ca,]
  con <- str_split_fixed(ca, "_", 3)[,3]
  con <- paste0(unique(tdta[,date]), "_", "DMSO_", con)
  blk <- unique(tdta[,date])
  tdta <- dat_mean[date_chnm_plate == ca | date == blk & date_chnm_plate == con]
  #tdta <- tdta[tdta$acsn %in% colnames(Covs[[blk]]),]
  tdta[chnm == "DMSO", conc_index := "CONC_DMSO"]
  tdta$id <- paste0(tdta$date_chnm_plate,"_X_",tdta$conc)
  tdta <- dcast.data.table(tdta, ... ~ steroid, value.var = c("mean_log_uM", "N"))
  indx <- order(tdta$conc)
  tdta <- tdta[indx,]
  colnames(tdta) <- gsub("^mean_log_uM_","",colnames(tdta))
  ## Hormones are in columns 7:ncol(tdta)
  ## Drop hormones which are NA in the bottom row
  refrow <- 1
  maxconc <- nrow(tdta)
  y0 <- data.matrix(tdta[1,c(hormones_sub), with = FALSE])
  y1 <- data.matrix(tdta[-1,c(hormones_sub), with = FALSE])
  y1 <- sweep(y1,2,y0,"-")
  y1 <- y1[, apply(y1, 2, function(x) !any(is.na(x)))]
  NH <- length(which(hormones_sub %in% colnames(y1)))
  CI11P <- solve(CovTP[which(rownames(CovTP) %in% colnames(y1)),
                       which(colnames(CovTP) %in% colnames(y1))])
  SSmeasured <- colnames(tdta)[colnames(tdta) %in% SSnms]
  maxsqrtN <- sqrt(apply(tdta[, SSmeasured, with = FALSE],1,max, na.rm=TRUE))
  out[[ca]] <- data.frame(Conc = tdta$conc[-1],
                          conc_index = tdta$conc_index[-1],
                          #D11 = apply(y1,1, function(x) sqrt((x %*% CI11[[blk]] %*% x)/NH)),
                          D11P = apply(y1,1, function(x) sqrt((x %*%  CI11P %*% x)/NH)),
                          #sD11 = apply(y1,1, function(x) sqrt(sum(x^2 * VI11[[blk]])/NH)),
                          #sD11P = apply(y1,1, function(x) sqrt(sum(x^2 *  VI11P[which(names(VI11P) %in% colnames(y1))])/NH)),
                          CA = rep(ca, nrow(y1)),
                          row.names=paste(rep(ca, nrow(y1)),tdta$conc[-1],sep="_X_"),
                          NH = NH,
                          N = max(maxsqrtN)^2)
}

#-------------------------------------------------------------------------------------------------#
#-----Save output as Dists
#-------------------------------------------------------------------------------------------------#

Dists <- data.table(do.call("rbind", out))
Dists[, chnm := str_split_fixed(CA, "_", 3)[,2]]
Dists$nChem <- make.names(Dists$chnm)
rownames(Dists) <- as.character(1:nrow(Dists))
Pvalues$BM.fdr <- p.adjust(Pvalues$Manova_block, method="fdr")

DMeasures <- c("D11P")


#-------------------------------------------------------------------------------------------------#
#-----Calculate the critical limit for the Mahalanobis distance for the whole data set
#-------------------------------------------------------------------------------------------------#

## Pairwise plot all the log10(distances) against each other
## On the diagonals, plot the empirical CDF, identifying the upper
## 90th percentile.
CritLim <- structure(numeric(length(DMeasures)), names=DMeasures)
for (i in seq_along(DMeasures)) {
  pdta <- data.frame(x = getDist(DMeasures[i]))
  pdta2 <- data.frame(x = CritLim[i] <- quantile(pdta$x, p=0.90, na.rm=TRUE))
  pdta3 <- data.frame(x = quantile(pdta$x, p=0.91, na.rm=TRUE),
                      y = 0.0,
                      text = as.character(signif(pdta2$x, digits=2)))
  ps <- ggplot() + stat_ecdf(data=pdta, aes(x=x), geom="step") +
    geom_vline(data=pdta2, aes(xintercept=x)) +
    geom_text(data=pdta3, aes(x=x,y=y,label=text), hjust=-0.1, vjust=0) +
    scale_x_log10() +
    theme_bw()
  p[i,i] <- ps
}

## Convert CritLim into a dataframe
CritLim <- data.frame(dist_type=names(CritLim),
                      dist=CritLim)

#-------------------------------------------------------------------------------------------------#
#-----Identify, for each test chemical, summary Mahalanobis distance statistics as well as the 
#-----number of concs above the critical limit
#-------------------------------------------------------------------------------------------------#

## How many D's > CritLim for each CA and D-type? What's the conc associated
## with the first D > 2 for each D-type?
CA2 <- unique(Dists$CA)
Resps <- as.data.frame(t(sapply(CA2, function(ca) {
  tdta <- Dists[CA == ca,]
  CA <- as.character(unique(Dists[CA == ca, CA]))
  Nconc <- sum(!is.na(tdta[,D11P]))
  NH <- max(tdta$NH)
  crit05 <- critval(0.05, (Nconc + 1), NH)$root
  crit01 <- critval(0.01, (Nconc + 1), NH)$root
  NgtCrit <- apply(sweep(tdta[,DMeasures, with = FALSE],2,CritLim$dist,"-"),2,function(x) sum(x > 0, na.rm=TRUE))
  names(NgtCrit) <- paste0("NgtCrit_",names(NgtCrit))
  
  minconc <- apply(sweep(tdta[,DMeasures, with = FALSE],2,CritLim$dist,"-"),2,
                   function(x) if (any(x[!is.na(x)] > 0)) min(tdta$Conc[x > 0], na.rm=TRUE) else Inf)
  names(minconc) <- paste0("minconc_",names(minconc))
  c(NgtCrit, minconc, Nconc = Nconc,
    Maxconc = max(tdta$Conc, na.rm=TRUE),
    meanD11P = mean(tdta$D11P, na.Rm=TRUE),
    maxD11P = max(tdta$D11P, na.rm=TRUE),
    ## maxT2stat = max(tdta$T2stat),
    meanD11P02 = mean(tdta$D11P^2, na.rm=TRUE),
    N = max(tdta$N),
    NH = NH,
    crit05 = crit05,
    crit01 = crit01,
    Scrit05 = scrit05 <- median(sqrt(crit05 / (tdta$N * tdta$NH))),
    Scrit01 = scrit01 <- median(sqrt(crit01 / (tdta$N * tdta$NH))),
    ## NgtCritT2_05 = sum(tdta$T2stat > crit05),
    ## NgtCritT2_01 = sum(tdta$T2stat > crit01),
    NgtScrit_05 = sum(tdta$D11P > scrit05),
    NgtScrit_01 = sum(tdta$D11P > scrit01),
    date_chnm_plate = CA)
})))

#sapply makes everything a factor, so convert all numeric columns to a numeric vector
Resps[, 1:15] <- sapply(Resps[, 1:15], function(x) as.numeric(levels(x))[x])
#sapply(Resps, class)

for (i in seq_along(DMeasures)) {
  Resps[,nm <- paste0("minD_",DMeasures[i])] <- numeric(nrow(Resps))
  tmp <- getDist(DMeasures[i])
  Resps[,nm] <- tmp[Resps$date_chnm_plate]
}


## For each D-type, how many Ngtcrit when D > CritLim?
for (i in seq_along(DMeasures)) {
  writeLines(paste("\n------------ ",DMeasures[i]))
  tmp1 <- table(Resps[(Resps[,paste0("minD_",DMeasures[i])] > CritLim[DMeasures[i], "dist"]),paste0("NgtCrit_",DMeasures[i])])
  #dimnames(tmp1)[[1]] <- as.character(as.numeric(dimnames(tmp1)[[1]]) - 1)
  N1 <- sum(tmp1)
  tmp2 <- table(Resps[(Resps[,paste0("minD_",DMeasures[i])] < CritLim[DMeasures[i], "dist"]),paste0("NgtCrit_",DMeasures[i])])
  N2 <- sum(tmp2)
  tbl <- matrix(0, nrow=2, ncol=max(dim(tmp1), dim(tmp2)),
                dimnames=list(c("D < crit", "D > crit"), sort(unique(unlist(c(dimnames(tmp1),dimnames(tmp2)))))))
  tbl[1,dimnames(tmp2)[[1]]] <- tmp2
  tbl[2,dimnames(tmp1)[[1]]] <- tmp1
  tbl[1,] <- tbl[1,] / N2
  tbl[2,] <- tbl[2,] / N1
  tbl <- cbind(tbl, N=c(N2,N1))
  print(signif(tbl, digits=3))
}

#-------------------------------------------------------------------------------------------------#
#-----combine Pvalues and Resps data frames
#-------------------------------------------------------------------------------------------------#

Pvalues_corrected <- Pvalues[chnm != "DMSO" & chnm != "Colchicine",]
Resps <- as.data.table(Resps)
Resps[, date_chnm_plate := as.character(date_chnm_plate)]

Pvalues_corrected <- Pvalues_corrected[order(date_chnm_plate),]
Resps <- Resps[order(date_chnm_plate),]

#Mahalanobis_dat <- cbind(Pvalues_corrected, Resps)
Mahalanobis_dat <- merge(Pvalues_corrected, Resps, by = "date_chnm_plate")

#add Cas numbers
casns <- unique(dat[spid != "DMSO",.(casn, date_chnm_plate)])
casns <- casns[order(date_chnm_plate)]

Mahalanobis_dat[, casn := casns[date_chnm_plate %in% Mahalanobis_dat[,date_chnm_plate], casn]]

Mahalanobis_dat <- Mahalanobis_dat[,c(1:2, 35, 3:34)]

#-------------------------------------------------------------------------------------------------#
#-----Save output
#-------------------------------------------------------------------------------------------------#

save(Mahalanobis_dat, CritLim, Resps, Dists, dat_mean, dat, Residuals, file=paste0("../RData/AllResps", Sys.Date(), ".RData"))
#save(Mahalanobis_dat, CritLim, Resps, Dists, dat_mean, dat, Residuals, file=paste0("../RData/AllResps_outliersRemoved", Sys.Date(), ".RData"))

#-------------------------------------------------------------------------------------------------#
#-----Make Supp. 9 table with Mahalanobis distances
#-------------------------------------------------------------------------------------------------#

Mahalanobis_dat_supp <- copy(Mahalanobis_dat[,c(1:3, 25, 32)])
Mahalanobis_dat_supp[, adj_maxmMD := maxD11P - Scrit01]

fwrite(Mahalanobis_dat_supp, file = paste0("Supp9_Global_OECD_Mahalanobis_distances_", Sys.Date(), ".txt"), sep = "\t")
