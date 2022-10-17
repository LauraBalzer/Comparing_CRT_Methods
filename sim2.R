# SIMULATION 2
# Using TMLEs compare estimation for cluster-level & indv-level effects
# when there is highly informative cluster-size
#   with an effect and under the null
#   breaking the matches used for analysis 

# Author: Laura B. Balzer

rm(list=ls())

library('nbpMatching')
library('ltmle')
source('gen_data.R')
source('MainFunctions.R')
source('Stage2_Functions.R')
source('Adapt_Functions.R')

set.seed(1)

sim <- 2

# scale of inference
goal <- 'aRR'
#goal <- 'RD'

# Indicator if there is an effect (vs. the null)
effect<- T
# Indicator break the match pairs used for randomization
break.match <- T

# number of iterations
nReps <- 500
# Number of clusters in the target population
nPop<- 1000
# Number of clusters in each study
J <- 20
# Average (sd) number of individuals per cluster
N.mean <- 400
N.sd <- 250

if(break.match){ 
  # consider more covariates if breaking the match
  cand.adj <- c('U', 'W1','W2')
} else{
  # consider fewer covariates if keeping the match in the analysis (fewer indpt units)
  cand.adj <- c('U', 'W2')
}
verbose <- F

file.name<- paste('Sim', sim, '_effect', effect, '_goal', goal, '_break', break.match,
                  '_J',J, '_nMean', N.mean, '_nSD', N.sd, '_nReps', nReps,"_nPop", nPop, '_',
                  paste(cand.adj, sep='', collapse=''), 
                  '_v',  format(Sys.time(), "%d%b%Y"),
                  '.Rdata', sep='')
print(file.name)

pop <- get.full.data(sim=sim, J=nPop, N.mean=N.mean, N.sd=N.sd, effect=effect, verbose=F)

truth <- get.truth(pop)

print(truth)

UNADJ <- ADJW1 <-  ADJ.AP <- ADJ.AP.CV <- QADJ <- GADJ <- NULL 

for(r in 1:nReps){
  X.all <- get.full.data(sim=sim, J=J, N.mean=N.mean, N.sd=N.sd, effect=effect, verbose=verbose)

  O <- subset(X.all, select=-c(Y1,Y0))
    
  # UNADJUSTED 
  unadj <- suppressWarnings( suppressMessages( 
                do.estimation(goal=goal, O=O, break.match=break.match) ))
  UNADJ <- rbind(UNADJ, extract.me(unadj, truth))
  
  # FIXED ADJUSTMENT FOR W1 in the outcome regressin
  adjW1 <- suppressWarnings( suppressMessages( 
                do.estimation(goal=goal, O=O, break.match=break.match, QAdj='W1', gAdj='U') ))
  ADJW1 <- rbind(ADJW1, extract.me(adjW1, truth))

  # TMLE WITH ADAPTIVE PRESPEC
  adj.AP <- suppressWarnings( suppressMessages( 
                do.estimation(goal=goal, O=O, do.data.adapt=T, break.match=break.match,
                              do.cv.variance = T, cand.adj=cand.adj)  ))
  
  # non.CV.inference
  std.ic <- adj.AP[,c('est', 'CI.lo', 'CI.hi', 'se','pval')]
  ADJ.AP <- rbind(ADJ.AP, extract.me(std.ic, truth))
  # CV.inference
  cv.ic <- adj.AP[,c('est.1', 'CI.lo.1', 'CI.hi.1', 'se.1','pval.1')]
  colnames(cv.ic) <- colnames(std.ic)
  ADJ.AP.CV <- rbind(ADJ.AP.CV, extract.me(cv.ic, truth))
  
  QADJ <- rbind(QADJ, t(adj.AP$QAdj))
  GADJ <- rbind(GADJ, t(adj.AP$gAdj))
  
  print(r)
  
}


file.name
colnames(QADJ) <- colnames(GADJ) <- rownames(adj.AP)

save(truth, UNADJ, ADJW1, ADJ.AP, ADJ.AP.CV, QADJ, GADJ, 
     generate.cluster.sim2, 
     file = file.name)

################



