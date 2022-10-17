# SIMULATION 1 - SUPPLEMENTARY ESTIMATORS//SETTINGS
# Using TMLEs compare estimation for cluster-level & indv-level effects
#   with an effect and under the null
#   breaking and keeping the matches used for analysis 

# Author: Laura B. Balzer

rm(list=ls())

library('nbpMatching')
library('ltmle')
source('gen_data.R')

source('MainFunctions.R')
source('Stage2_Functions.R')
source('Adapt_Functions.R')
# source('Estimators.R')

set.seed(1)
# simulation indicator
sim <- 1 
nReps <- 500
# Indicator if there is an effect (vs. the null)
effect<- F
# Indicator break the match pairs used for randomization
break.match <- T

# Number of clusters in the target population
nPop<- 2500

# Number of clusters in each study
J <- 10
# Average number of individuals per cluster
N.mean <- 150
N.sd <- 80
# OCT 2022 - update to examine impact of limited cluster sizes
# N.mean <- 20
# N.sd <- 10

# print output
verbose<- F

# consider more candidates if breaking the match (supported by more inpt units)
if(break.match){
  cand.adj <- c('U', 'W1','W2')
} else{
  cand.adj <- c('U', 'W2')
}
# OCT 2022 - add in dummy variables W3 and W4
# cand.adj <- c('U','W1','W2','W3','W4')

# SUPP FILES ONLY FOCUS ON TMLES WITH AND WITHOUT PAIR-MATCHING 
# ANALOGOUS TO SIMULATION 2
file.name<- paste('SUPP_Sim', sim, '_effect', effect, '_break', break.match,
                  '_J',J, '_nMean', N.mean, '_nSD', N.sd, '_nReps', nReps,"_nPop", nPop, '_',
                  paste(cand.adj, sep='', collapse=''), 
                  '_v',  format(Sys.time(), "%d%b%Y"),
                  '.Rdata', sep='')
print(file.name)


pop <- get.full.data(sim=sim, J=nPop, N.mean=N.mean, N.sd=N.sd, effect=effect, verbose=F)

truth <- get.truth(pop)

print(truth)

UNADJ <- ADJ.AP <- ADJ.AP.CV <- QADJ <- GADJ <- NULL 

for(r in 1:nReps){
  
  X.all <- get.full.data(sim=sim, J=J, N.mean=N.mean, N.sd=N.sd, effect=effect, verbose=verbose)
  
  O <- subset(X.all, select=-c(Y1,Y0))
  
  # UNADJUSTED EFFECT ESTIMATOR
  unadj <- suppressWarnings( suppressMessages( 
                do.estimation(O=O, break.match=break.match) ))
  UNADJ <- rbind(UNADJ, extract.me(unadj, truth))

  # TMLES WITH ADAPTIVE PRE-SPEC
  adj.AP <- suppressWarnings( suppressMessages( 
                do.estimation(O=O, do.data.adapt=T, break.match=break.match,
                              do.cv.variance = T, cand.adj=cand.adj)  ))
  
  # non.CV.inference
  std.ic <- adj.AP[,c('est', 'CI.lo', 'CI.hi', 'se','pval')]
  ADJ.AP <- rbind(ADJ.AP, extract.me(std.ic, truth))
  
  # CV.inference (for comparison; results not shown)
  cv.ic <- adj.AP[,c('est.1', 'CI.lo.1', 'CI.hi.1', 'se.1','pval.1')]
  colnames(cv.ic) <- colnames(std.ic)
  ADJ.AP.CV <- rbind(ADJ.AP.CV, extract.me(cv.ic, truth))
  
  # selected adjustment sets
  QADJ <- rbind(QADJ, t(adj.AP$QAdj))
  GADJ <- rbind(GADJ, t(adj.AP$gAdj))
  
  print(r)
  
}

colnames(QADJ) <- colnames(GADJ) <- rownames(adj.AP)

save(truth, UNADJ, ADJ.AP, ADJ.AP.CV, QADJ, GADJ, 
     generate.cluster.sim1, 
     file = file.name)

