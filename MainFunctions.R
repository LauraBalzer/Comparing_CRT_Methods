## Author: Laura Balzer 

# Wrapper functions to do the estimation (including aggregation and weights) 
# and inference 
#--------------------

# Standard outputs: 


# do.estimation.sim1: estimation and inference with 
# - unadjusted effect estimator: cluster-level TMLE adjusting for nothing; targeting cluster-level effect
# - dataC.effectC: cluster-level TMLE for the cluster-level effect
# - dataI.effectC: Hierarchical TMLE for the cluster-level effect 
# - ttest: standard t-test on log cluster-level outcomes (estimates geometric risk ratio)
# - CARE: covariate adjusted residuals (estimates geometric risk ratio)
# - GEE: generalized estimating equations for the indv-level risk ratio
# - Augmented-GEE: for the indv-level risk ratio

do.estimation.sim1 <- function(O, ind.cov, break.match=T, verbose=F){
  
  #--------------------
  # CLUSTER-LEVEL EFFECT
  
  # aggregate with the empirical mean & then add weights for effect: cluster
  Occ <- get.Yc(O=O, target='cluster')
  
  
  unadj <- Stage2(data.input=Occ, target='clust', QAdj=NULL, gAdj=NULL, 
                  do.data.adapt =F,  do.cv.variance=F, cand.adj=NULL, 
                  break.match=T, verbose=verbose)
  
  # OCT 2022 - consider expanded set
  cand.adj <- c('U', 'W1', 'W2', 'W3', 'W4')
  
  dataC.effectC <- Stage2(data.input=Occ, target='clust', QAdj=NULL, gAdj=NULL, 
                          do.data.adapt =T,  do.cv.variance=F, cand.adj=cand.adj, 
                          break.match=T, verbose=verbose)
  
  # using individual-level data
  O$alpha <- 1/O$n
  dataI.effectC <- Stage2(data.input=O, target='clust', QAdj=NULL, gAdj=NULL, 
                          do.data.adapt =T,  do.cv.variance=F, cand.adj=cand.adj, 
                          break.match=T, verbose=verbose)
  
  #--------------------
  # GEOMETRIC RISK RATIO 
  ttest <- output.ttest(psi=NA, 
                        g=t.test(x=log(Occ[Occ$A==1, 'Y']),
                                 y=log(Occ[Occ$A==0, 'Y']), var.equal=T), 
                        gRR=T, paired=F)
  
  
  
  care <- do.care(gRR=T, train=O, ind.cov=ind.cov, psi=NA, paired=F)
  
  #--------------------
  # INDV-EFFECT
  # Oct2022 update: do Fay/Graubard small sample correction 
  #gee <- do.gee(train=O, ind.cov=ind.cov, psi=NA, paired=F, link='poisson')
  gee <- do.gee.fay(train=O, ind.cov=ind.cov, psi=NA, paired=F, link='poisson')
  aug.gee <- do.gee.aug(train=O, ind.cov=ind.cov, psi=NA, link='poisson')
  
  #--------------------
  YAY <- data.frame(rbind(unadj=unadj[,colnames(care)],
                          dataC.effectC=dataC.effectC[,colnames(care)], 
                          dataI.effectC=dataI.effectC[,colnames(care)],
                          ttest=ttest,
                          care=care,
                          gee=gee, aug.gee=aug.gee
  ))
  
  # RETURN ADJUSTMENT SET FOR TMLEs 
  ADJ <- data.frame(cbind( QC = dataC.effectC$QAdj, gC = dataC.effectC$gAdj,
                           QI = dataI.effectC$QAdj, gI = dataI.effectC$gAdj
  ))
  
  list(YAY=YAY, ADJ=ADJ)
}



# do.estimation: Wrapper function implemente TMLEs (including aggregation and weights) 
# for PTBi analysis and for Simulation 2

# Standard outputs: 
#   dataC.effectC: cluster-level TMLE estimate cluster-level effect
#   dataI.effectC: Hierarchical TMLE (with indv-level data) to estimate cluster-level effect                           
#   dataC.effectI: cluster-level TMLE to estimate ind-level effect
#   dataI.effectI: Hierarchical TMLE (with ind-level data) to estimate ind-level effect

do.estimation <- function(goal='aRR', O, QAdj=NULL, gAdj=NULL, 
                          do.data.adapt =F, do.cv.variance=F, cand.adj=NULL, 
                          break.match=T, verbose=F){

  
  #--------------------
  # CLUSTER-LEVEL EFFECT
  
  # using cluster-level data
  # aggregate with the empirical mean & then add weights for effect: cluster
  Occ <- get.Yc(O=O, target='cluster')
  dataC.effectC <- Stage2(goal=goal, data.input=Occ, target='clust', 
                          QAdj=QAdj, gAdj=gAdj, do.data.adapt =do.data.adapt, 
                          do.cv.variance=do.cv.variance, cand.adj=cand.adj, 
                          break.match=break.match, verbose=verbose)
  
  # using individual-level data
  O$alpha <- 1/O$n
  dataI.effectC <- Stage2(goal=goal, data.input=O, target='clust', 
                          QAdj=QAdj, gAdj=gAdj, do.data.adapt =do.data.adapt, 
                          do.cv.variance=do.cv.variance, cand.adj=cand.adj, 
                          break.match=break.match, verbose=verbose)
  
  
  #---------------------------
  # INDIVIDUAL-LEVEL EFFECT 
  
  # using cluster-level data 
  # aggregate with the empirical mean & then add weights for effect: indv
  Oci <- get.Yc(O=O, target='indv') #weights= J/nTot*Nj
  dataC.effectI <- Stage2(goal=goal, data.input=Oci, target='indv', 
                          QAdj=QAdj, gAdj=gAdj, do.data.adapt =do.data.adapt, 
                          do.cv.variance=do.cv.variance, cand.adj=cand.adj, 
                          break.match=break.match, verbose=verbose)
  
  # using individual-level data 
  O$alpha <- 1
  dataI.effectI <- Stage2(goal=goal, data.input=O,  target='indv', 
                          QAdj=QAdj, gAdj=gAdj, do.data.adapt =do.data.adapt, 
                          do.cv.variance=do.cv.variance, cand.adj=cand.adj, 
                          break.match=break.match, verbose=verbose)
  
  
  YAY <- data.frame(rbind(dataC.effectC=dataC.effectC, dataI.effectC=dataI.effectC, 
                          dataC.effectI=dataC.effectI, dataI.effectI=dataI.effectI))
  YAY
}


# HELPER FUNCTIONS


# get.Yc: aggregates data to the cluster-level
# add weights appropriate for the target of inference
get.Yc <- function(O, target='cluster'){
  
  # aggregate to the cluster by taking an empirical mean in each cluster 
  # equivalent to weighted sum with 1/N_j
  Oc <- aggregate(O, by=list(O$id), mean)[,-1]
  J <- nrow(Oc)
  nTot <- sum(Oc$n)
  
  if(target!='cluster'){
    # if target is 'indv', then need weights: J/nTot*1/alpha_ij
    # where alpha_ij is the weights used in the weighted sum (i.e )
    Oc$alpha <- J/nTot*Oc$n
  }
  # target='cluster', then no action
  
  # weights should sum up to the number of clusters
  # print( sum(Oc$alpha) ) 
  Oc
}


# extract.me.goal: return point estimate, standard error estimate, 
# CI covers truth, and  reject null
#   with respect to the target of inference
extract.me.goal <- function(est, psi){
  psi.hat <- est[,'est']
  data.frame(cbind( 
    est=psi.hat,  se = est[,'se'],
    cover= (psi >= est[,'CI.lo'] & psi <= est[,'CI.hi']),
    reject = as.numeric( est[,'pval'] < 0.05)
  ))
  
}


# extract.me: wrapper function to summarize results of simulation study 2
extract.me <- function(est.matrix, truth){
  
  data.frame(cbind(
    dCeC = extract.me.goal(est=est.matrix['dataC.effectC',], psi=truth$RR.clust),
    dIeC = extract.me.goal(est=est.matrix['dataI.effectC',], psi=truth$RR.clust),
    dCeI = extract.me.goal(est=est.matrix['dataC.effectI',], psi=truth$RR.ind),
    dIeI = extract.me.goal(est=est.matrix['dataI.effectI',], psi=truth$RR.ind)
  ))
  
  
}



# extract.ltmle - make ltmle output match output from Stage2
extract.ltmle <- function(est, psi){
  extract.extract <- function(est.g){
    Z <- data.frame( est.g$estimate, est.g$CI, est.g$std.dev)
    colnames(Z) <- c('est', 'CI.lo', 'CI.hi', 'se')
    Z
  }
  data.frame(cbind(Txt= extract.extract(est$treatment),
                   Con=extract.extract(est$control),
                   extract.extract(est$RR), pval=est$RR$pvalue, QAdj='U', gAdj='U'
  ))
}




