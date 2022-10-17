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

goal <- 'aRR'
#goal <- 'RD'

# UNADJUSTED
# break matches - point ests & se should match ltmle
unadjB <- suppressWarnings( suppressMessages( do.estimation(goal=goal, O=O, break.match=T) ))
# keep matches - 
unadjP <- suppressWarnings( suppressMessages( do.estimation(goal=goal,O=O, break.match=F) ))

round( unadjB[,c('Txt.est', 'Con.est', 'est', 'CI.lo', 'CI.hi','pval')], 2)
round( unadjP[,c('Txt.est', 'Con.est', 'est', 'CI.lo', 'CI.hi','pval')], 2)

# ADAPTIVE PRE-SPECIFIC FOR {U, W2 }
# subset the dataset on those with measured covariates (C-section status ==W2)
O2 <- O[!is.na(O$W2),]
dim(O); dim(O2)
# unadjB2 <- suppressWarnings( suppressMessages( do.estimation(goal=goal,O=O2, break.match=T) ))
# unadjP2 <- suppressWarnings( suppressMessages( do.estimation(goal=goal,O=O2, break.match=F) ))

# break matches
tmleB <- suppressWarnings( suppressMessages(
  do.estimation(goal=goal, O=O2, do.data.adapt=T, do.cv.variance=T, cand.adj=c('U','W2'), break.match=T) ))
# keep matches - 
tmleP <- suppressWarnings( suppressMessages(
  do.estimation(goal=goal, O=O2, do.data.adapt=T, do.cv.variance=T, cand.adj=c('U','W2'), break.match=F) ))

file.name<- paste('PTBi_goal_', goal,
                  '_v', 
                  format(Sys.time(), "%d%b%Y"), 
                  '.Rdata', sep='')

save(unadjB, unadjP, tmleB, tmleP, file=file.name)

this.order <- c('dataC.effectC',  'dataC.effectI','dataI.effectC', 'dataI.effectI')


round( tmleB[,c('Txt.est', 'Con.est', 'est', 'CI.lo', 'CI.hi','pval')], 2)
round( tmleP[,c('Txt.est', 'Con.est', 'est', 'CI.lo', 'CI.hi','pval')], 2)


round( (unadjB$se^2)/(tmleB$se^2), 1)
