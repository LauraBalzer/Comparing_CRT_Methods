# SIMULATION 1 -
# Compare CRT estimators for their natural target of inference 
#   with an effect and under the null
#   breakingthe matches used for analysis 

# Author: Laura B. Balzer
rm(list=ls())

library('nbpMatching')
library('ltmle')
source('gen_data.R')

source('MainFunctions.R')
source('Stage2_Functions.R')
source('Adapt_Functions.R')
source('Estimators.R')


set.seed(1)

nReps <- 500
# Indicator if there is an effect (vs. the null)
effect<- F
# Number of clusters in the target population
nPop<- 2500

# Number of clusters in each study
J <- 20
# Average number of individuals per cluster
N.mean <- 150
N.sd <- 80
sim <- 1 


# individual-level adjustment variables for GEE, Aug-GEE, and CARE
ind.cov <- c('W1','W2', 'E1','E2')

file.name<- paste('Sim', sim, '_effect', effect, 
                  '_J',J, '_nMean', N.mean, '_nSD', N.sd, '_nReps', nReps,"_nPop", nPop,
                  '_v', 
                  format(Sys.time(), "%d%b%Y"), 
                  '.Rdata', sep='')
print(file.name)

#load(file.name)
verbose=F


pop <- get.full.data(sim=sim, J=nPop, N.mean=N.mean, N.sd=N.sd, effect=effect, verbose=F)
pop.C <- aggregate(pop, by=list(pop$id), mean)

truth <- data.frame( cbind(
                     int.ind=mean(pop$Y1),
                     con.ind=mean(pop$Y0),
                     RR.ind= mean(pop$Y1)/mean(pop$Y0),
                     int.clust= mean(pop.C$Y1),
                     con.clust=mean(pop.C$Y0),
                     RR.clust=mean(pop.C$Y1)/mean(pop.C$Y0),
                     # also need geometric incidence ratio 
                     gRR.clust = exp( mean(log(pop.C$Y1)) - mean(log(pop.C$Y0)) ),
                     getMeasuresVariability(X1=pop.C$Y1, X0=pop.C$Y0, pairs=pop.C$pair)
                     ))
print( truth  )


UNADJ <- TMLE.CC <- TMLE.CI <- TTEST <- CARE <- GEE <- A.GEE <- ADJ <-   NULL 

for(r in 1:nReps){
  X.all <- get.full.data(sim=sim, J=J, N.mean=N.mean, N.sd=N.sd, effect=effect, verbose=verbose)
  
  O <- subset(X.all, select=-c(Y1,Y0))
  
  temp <- suppressWarnings( suppressMessages( do.estimation.sim1( O=O, ind.cov=ind.cov, break.match=T, verbose=F) ))
  # temp
  
  # cluster-level effect
  UNADJ <- rbind(UNADJ, extract.me.goal(temp$YAY['unadj',], psi=truth$RR.clust))
  TMLE.CC <- rbind(TMLE.CC, extract.me.goal(temp$YAY['dataC.effectC',], psi=truth$RR.clust))
  TMLE.CI <- rbind(TMLE.CI, extract.me.goal(temp$YAY['dataI.effectC',], psi=truth$RR.clust))
  
  # geometric scale
  TTEST <- rbind(TTEST, extract.me.goal(temp$YAY['ttest',], psi=truth$gRR.clust))
  CARE <- rbind(CARE, extract.me.goal(temp$YAY['care',], psi=truth$gRR.clust))
  
  # ind-effect
  GEE <- rbind(GEE, extract.me.goal(temp$YAY['gee',], psi=truth$RR.ind))
  A.GEE <- rbind(A.GEE, extract.me.goal(temp$YAY['aug.gee',], psi=truth$RR.ind))
  
  # adjustment set fot TMLE
  ADJ <- rbind(ADJ, temp$ADJ)
  print(r)
  
}


save(truth, ind.cov, UNADJ, TMLE.CC, TMLE.CI, TTEST, CARE , GEE , A.GEE ,ADJ,
     file = file.name)


