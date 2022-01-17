# PTBi ANALYSIS
# Using TMLEs compare estimation for cluster-level & indv-level effects
# in the real data set

# Author: Laura B. Balzer
rm(list=ls())

set.seed(1)

library('ltmle')
source('MainFunctions.R')
source('Stage2_Functions.R')
source('Adapt_Functions.R')

dt <- read.csv('ptbi.csv')
# add dummy for intercept-only regression + weights
dt$alpha <- dt$U <- 1
O <- NULL
# add cluster size
for(j in 1:20){
  scotts.tots <- dt[dt$id %in% j,]
  scotts.tots$n <- sum(scotts.tots$U)
  O <- rbind(O, scotts.tots)
}

O <- O[,c('id','pair','U','alpha','n', 'W2','A','Y')]

# UNADJUSTED
# break matches - point ests & se should match ltmle
unadjB <- suppressWarnings( suppressMessages( do.estimation(O=O, break.match=T, return.ltmle=T) ))
# keep matches - 
unadjP <- suppressWarnings( suppressMessages( do.estimation(O=O, break.match=F) ))

# ADAPTIVE PRE-SPECIFIC FOR {U, W2 }
# subjset dataset on those with measured covariates (C-section status ==W2)
O2 <- O[!is.na(O$W2),]
dim(O); dim(O2)
# break matches
tmleB <- suppressWarnings( suppressMessages(
  do.estimation(O=O2, do.data.adapt=T, do.cv.variance=T, cand.adj=c('U','W2'), break.match=T) ))
# keep matches - 
tmleP <- suppressWarnings( suppressMessages(
  do.estimation(O=O2, do.data.adapt=T, do.cv.variance=T, cand.adj=c('U','W2'), break.match=F) ))



save(unadjB, unadjP, tmleB, tmleP, file='PTBi_results.Rdata')

