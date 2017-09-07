#-----------------------------------------------------------------------------------#
# Haggard et al. "Using a high-throughput steroidogenesis assay to predict effects on 
# ...estrogen and androgen production or overall perturbation of the steroid biosynthesis pathway"


#This code will generate the residual matrices from the H295R
#MC data for use in the Mahalanobis distance calculations


#-----------------NOTES-----------------#
##The loaded H295R master data file has been corrected for LLOQ/sqrt(2) for ND/NQ values and all NR values were removed##
##Similarly, Forskolin and Prochloraz controls were removed##
##rvals have also already been converted to uM##
#---------------------------------------#

rm(list = ls())

library(data.table)
library(multcomp)
library(stringr)
library(reshape2)

setwd("../Mahalanobis Distance")

#-------------------------------------------------------------------------------------------------#
#-----Functions used in the script
#-------------------------------------------------------------------------------------------------#

Ngood <- function(X) {
  inx <- intersect(colnames(X), hormones_sub)
  x2 <- X[apply(X,1,function(x) all(!is.na(x[inx]))),]
  clF <- unique(x2$classFull)
  clF <- clF[-grep("^DMSO",clF)]
  length(clF)
}

delcols <- function(x) {
  CCols <- intersect(colnames(x), hormones_sub)
  lastNgood <- Ngood(x)
  lasti <- NA
  lastj <- NA
  for (i in CCols)
    for (j in CCols) {
      tmp <- Ngood(x[,-match(unique(c(i,j)),colnames(x)), with = FALSE])
      if (tmp > lastNgood) {
        lastNgood <- tmp
        lasti <- i
        lastj <- j
      }
    }
  list(ngood=lastNgood, dropcols=c(lasti, lastj))
}

#-------------------------------------------------------------------------------------------------#
#-----Load H295R master data file
#-------------------------------------------------------------------------------------------------#
load("../RData/H295R_master_table_2017-08-08.RData")

dat <- copy(dat0_combined)

#-------------------------------------------------------------------------------------------------#
#-----Give DMSO control a chnm
#-------------------------------------------------------------------------------------------------#

dat[spid == "DMSO", chnm := "DMSO"]

#-------------------------------------------------------------------------------------------------#
#-----fix Triadimenol to have two chnms reflecting different conc ranges
#-------------------------------------------------------------------------------------------------#

dat[spid == "TP0000884B12" | spid == "TP0000884A04" | spid == "TP0000883A04", chnm := "Triadimenol2"]

#-------------------------------------------------------------------------------------------------#
#-----Extract date component of each apid and reformat to a consistent format
#-------------------------------------------------------------------------------------------------#

dat[, date := str_split_fixed(apid, "\\.", n = 2)[,1]]
dat[, plate := str_split_fixed(apid, "\\.", n = 2)[,2]]
dat[, plate := str_replace(plate, "\\.", "")]

reform <- c("02Apr14" = "20140402",
            "02Apr2014" = "20140402",
            "03Sep2014" = "20140903",
            "09Apr2014" = "20140409",
            "10Jun2015" = "20150610",
            "20130320" = "20130320",
            "20130321" = "20130321",
            "26Mar2014" = "20140326",
            "04112017" = "20170411")

dat[, date := reform[date]]

#-------------------------------------------------------------------------------------------------#
#-----make chnm_date variable since the date will be the blocking variable for the MANOVA
#-------------------------------------------------------------------------------------------------#

dat[, chnm := str_replace(chnm, "PharmaGSID_", "PharmaGSID")]

dat[, chnm_date := paste0(chnm, "_", date)]

dat[, date_chnm_plate := paste(date, chnm, plate, sep = "_")]

dat[, chnm_plate := paste0(chnm, "_", plate)]

#-------------------------------------------------------------------------------------------------#
#-----log uM
#-------------------------------------------------------------------------------------------------#

dat[, log_uM := log(uM)]

#-------------------------------------------------------------------------------------------------#
#-----fix rounding errors for chemcial concentrations
#-------------------------------------------------------------------------------------------------#

dat[conc >= 99.9, conc := 100]

#fix odd Triclosan conc
dat[chnm == "Triclosan" & conc == 0.004115226, conc := 0.004] 

dat[, conc := signif(conc, digits = 2)]

#dat[, fixed_conc := signif(conc, digits = 2)]
#conc_count <- dat[, .(sum_conc = length(unique(conc)), sum_fixed_conc = length(unique(fixed_conc))), by = chnm]

#-------------------------------------------------------------------------------------------------#
#-----Select which hormones we want to use
#-------------------------------------------------------------------------------------------------#

hormones_all <- unique(dat[,steroid])

#remove Pregnenolone and DHEA
hormones_sub <- hormones_all[c(-8, -11)]

#-------------------------------------------------------------------------------------------------#
#-----Make empty Pvalues data table
#-------------------------------------------------------------------------------------------------#

Pvalues <- data.table(chnm=dat[,str_split_fixed(unique(date_chnm_plate), "_", n = 3)[,2]], #use date_chnm_plate as blocking variable
                      date_chnm_plate=unique(dat[,date_chnm_plate]),
                      Manova_block=numeric(length(unique(dat[,date_chnm_plate]))),
                      IndivAdjusted=numeric(length(unique(dat[,date_chnm_plate]))),
                      H1=numeric(length(unique(dat[,date_chnm_plate]))),
                      H2=numeric(length(unique(dat[,date_chnm_plate]))),
                      H3=numeric(length(unique(dat[,date_chnm_plate]))),
                      H4=numeric(length(unique(dat[,date_chnm_plate]))),
                      H5=numeric(length(unique(dat[,date_chnm_plate]))),
                      H6=numeric(length(unique(dat[,date_chnm_plate]))),
                      H7=numeric(length(unique(dat[,date_chnm_plate]))),
                      H8=numeric(length(unique(dat[,date_chnm_plate]))),
                      H9=numeric(length(unique(dat[,date_chnm_plate]))),
                      H10=numeric(length(unique(dat[,date_chnm_plate]))),
                      H11=numeric(length(unique(dat[,date_chnm_plate]))),
                      H12=numeric(length(unique(dat[,date_chnm_plate]))),
                      H13=numeric(length(unique(dat[,date_chnm_plate]))),
                      stringsAsFactors=FALSE)

Pvalues <- Pvalues[order(date_chnm_plate),]
colnames(Pvalues)[5:17] <- hormones_all
rownames(Pvalues) <- unique(dat[order(date_chnm_plate),date_chnm_plate])

setkey(Pvalues, date_chnm_plate)

#-------------------------------------------------------------------------------------------------#
#-----FOR loop for MANOVA by block (date)
#-------------------------------------------------------------------------------------------------#

Blocks <- unique(dat$date)
Models <- vector("list", length=length(Blocks))
names(Models) <- Blocks

Residuals <- vector("list", length=length(Blocks))
names(Residuals) <- Blocks

for (block in unique(dat[, date])) {
  dt <- copy(dat[date == block, c("conc", "steroid", "chnm","log_uM", "date", "coli", "rowi", "wllt","date_chnm_plate", "plate"), with = FALSE])
  dt$classFull <- factor(ifelse(dt$wllt == "t", paste0(make.names(dt[,date_chnm_plate]),"_",dt[,conc]), paste0("DMSO_",dt[,wllt])))
  chnm_apids <- unique(dt[,date_chnm_plate])
  chnm_apids <- chnm_apids[-grep("DMSO", chnm_apids)]
  for (i in seq_along(chnm_apids)) {
    chnmapid <- chnm_apids[i]
    dt[,paste0("Reduced_",i)] <-
      factor(ifelse(dt$date_chnm_plate != chnmapid & dt$wllt == "t", paste0(make.names(dt$date_chnm_plate), "_", dt$conc),
                    paste0("DMSO_n")))
  }
  #dt[, ID := paste0(coli, "_", rowi, "_", chnm)]
  dt[, ID := paste0(coli, "_", rowi, "_", chnm, "_", plate)]
  #dtw <- reshape(dt, v.names="log_uM", timevar="steroid", idvar="ID", direction="wide", drop=c("coli","rowi"))
  #dtw <- dcast.data.table(dt[,c(-6, -7, -10)], ID + ... ~ steroid, value.var = "log_uM", fun.aggregate = mean)
  dtw <- dcast.data.table(dt[,c(-6, -7)], ID + ... ~ steroid, value.var = "log_uM")
  colnames(dtw) <- gsub("^log_uM\\.","",colnames(dtw))
  dtw <- dtw[,-match(c("DHEA","Pregnenolone"), colnames(dtw)), with = FALSE]
  ## For the Manova, drop columns that have an excessive number of NA's.
  ## Criterion is to maximize the number of chemxconcxapid levels.
  out <- delcols(dtw)
  if (any(!is.na(out$dropcols))) {
    dtw2 <- dtw[,-match(out$dropcols[!is.na(out$dropcols)], colnames(dtw)), with = FALSE]
    LHS <- paste0("cbind(",paste(paste0("`", intersect(hormones_sub, colnames(dtw2)),"`"), collapse=", "),")")
  } else {
    dtw2 <- dtw
    LHS <- paste0("cbind(",paste(paste0("`", hormones_sub, "`"), collapse=", "),")")
  }
  ## Initial test to find bad residuals
  zz <- list()
  Formula <- as.formula(paste0(LHS," ~ 0 + classFull"))
  zz[[1]] <- try(lm(Formula, data=dtw2), silent=TRUE)
  
  #filter dtw to remove singleton samples and outlier paired samples with standard deviation of the residuals > 1
  inx <- intersect(colnames(dtw2), hormones_sub)
  x2 <- dtw2[apply(dtw2,1,function(x) all(!is.na(x[inx]))),] #NAs are ommitted and allows for proper cbind of residuals
  residual_list <- cbind(data.table(zz[[1]]$residuals), ID = x2[,ID], date_chnm_plate = x2[,date_chnm_plate], conc = x2[,conc])
  residual_list[, conc_date_chnm_plate := paste(conc, date_chnm_plate, sep = "_")]
  
  #count up the number of replicates for each conc_date_chnm_plate
  for(x in unique(residual_list[, conc_date_chnm_plate])){
    residual_list[conc_date_chnm_plate == x, reps := residual_list[conc_date_chnm_plate == x, .N]]
  }
  
  residual_list_filtered <- residual_list[reps > 1,] #remove singleton samples
  
  #calculate standard deviation of the residuals
  residual_diff_data <- residual_list_filtered[, lapply(.SD, function(x) sd(x)),
                                               by = .(conc_date_chnm_plate, date_chnm_plate),
                                               .SDcols = colnames(residual_list_filtered)[-c((length(residual_list_filtered)-4):length(residual_list_filtered))]]
  
  
  residual_diff_data[, maxSD := ifelse(apply(.SD, 1, function(x) max(x)) > 1, 1, 0), .SDcols = colnames(residual_diff_data)[3:length(residual_diff_data)]]
  
  Residuals[[block]] <- residual_diff_data
  
  #make new filtered dtw data table
  dtw3 <- copy(dtw2)
  dtw3[, conc_date_chnm_plate := paste(conc, date_chnm_plate, sep = "_")]
  dtw3 <- dtw3[conc_date_chnm_plate %in% residual_diff_data[maxSD == 0, conc_date_chnm_plate],] #filter out bad data
  dtw3[, conc_date_chnm_plate := NULL]
  
  ## Real Tests
  zz <- list()
  Formula <- as.formula(paste0(LHS," ~ 0 + classFull"))
  zz[[1]] <- try(lm(Formula, data=dtw3), silent=TRUE)
  
  Models[[block]] <- list(fit = zz[[1]],
                        Hormones = if (any(!is.na(out$dropcols))) setdiff(names(hormones_all), c(out$dropcols, "Pregnenolone","DHEA")) else setdiff(names(hormones_all), c("Pregnenolone","DHEA")))
  
  if (!is(zz[[1]], "try-error")) {
    for (i in 1:(length(chnm_apids))) {
      Formula <- as.formula(paste0(LHS," ~ 0 + Reduced_",(i)))
      zz[[i+1]] <- lm(Formula, data=dtw3)
    }
    ## Now, save the P-values for the tests of individual chemicals
    for (i in 2:length(zz)) {
      out <- try(anova(zz[[i]], zz[[1]], test="Pillai"), silent=TRUE)
      Pvalues[date_chnm_plate == chnm_apids[i-1], Manova_block :=  ifelse(!is(out, "try-error"), out$`Pr(>F)`[2], NA)]
    }
  } else {
    Pvalues[chnm_apids, Manova_block := NA]
  }
  ## Now, do the individual hormone tests for each chemical
  for (ch in chnm_apids) {
    rwnm <- ch
    for (hm in hormones_sub) {
      dts <- dt[steroid == hm & date_chnm_plate == ch | steroid == hm & wllt == "n",]
      dts <- dts[plate == unique(dts[steroid == hm & date_chnm_plate == ch, plate]),]
      dts$classFull <- factor(ifelse(dts$wllt == "t", paste0("X_",dts$conc), paste0("DMSO_",dts$wllt)))
      dts$classReduced <- factor(ifelse(dts$wllt == "t", paste0("DMSO_n"), paste0("DMSO_",dts$wllt)))
      zz1 <- try(lm(log_uM ~ classFull, data=dts))
      if (!is(zz1, "try-error")) {
        zz0 <- lm(log_uM ~ 1, data=dts)
        out <- try(anova(zz1,zz0))
        Pvalues[date_chnm_plate == rwnm, which(colnames(Pvalues) == hm)] <- if(!is(out, "try-error")) out$`Pr(>F)`[2] else NA
      }
    }
    Pvalues[date_chnm_plate == rwnm, IndivAdjusted := min(p.adjust(Pvalues[date_chnm_plate == rwnm, hormones_sub, with = FALSE], method="holm"),na.rm=TRUE)]
  }
}

Residuals <- do.call(rbind, c(l = Residuals, fill = TRUE))

#-------------------------------------------------------------------------------------------------#
#-----save output for Mahalanobis distance calculation
#-------------------------------------------------------------------------------------------------#

save(dat, Pvalues, Models, Residuals, file = paste0("../RData/manova_output", Sys.Date(), ".RData"))

