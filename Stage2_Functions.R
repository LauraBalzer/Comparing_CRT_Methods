##############
# Author: Laura B. Balzer
#
# Stage2_Functions.R 
# R code to implement all Stage2 analyses to compare intervention effects between arms
# Using TMLE with and without Adaptive Prespecification 

# Modified from https://github.com/LauraBalzer/SEARCH_Analysis_Adults

# Stage2: Main function for estimation and inference 
# input: 
#	  goal;  aRR= arithmetic risk ratio; RD=risk difference; OR= odds ratio (not recommended)
#	  target; "clust" for cluster-level effect; otherwise indv-level effect
#	  data.input: observed data,
#	  QAdj: prespecified adjustment variables for estimation of the conditional mean outcome with main terms (e.g. Y~A+W1+W2), 
#	  gAdj: prespecified adjustment variables for estimation of the propensity score with main terms (e.g. A~W1+W2),
#	  do.data.adapt: indicator to do adaptive prespecification,
#   do.cv.variance: indicator to get cross-validated variance only applies if do.data.adapt=T,
#	  cand.adj: set of candidate adjustment variables if do.data.adapt=T,
#	  break.match: indicator to break the match pairs during the analysis,
#	  verbose: indicator to print updates 
#   psi: true value of effect (if known)
#
# output: point estimate and inference

Stage2 <- function(goal='aRR', target='clust', data.input, QAdj=NULL, gAdj=NULL, 
                   do.data.adapt =F, do.cv.variance=F, cand.adj=NULL, 
                   break.match=F, verbose=F, psi=NULL){	
  
  
  if(do.data.adapt){
    # implement adaptive pre-specification
    select <- do.adaptive.prespec(goal=goal, target=target, break.match = break.match,
                                 Ldata= data.input, clust.adj=cand.adj, QAdj=QAdj, gAdj=gAdj)
    
    QAdj = select$QAdj
    gAdj = select$gAdj	
  }
  
  # Run full TMLE algorithm with prespecified adjustment set
  est <- do.TMLE(goal=goal, target=target, train=data.input, QAdj=QAdj, 
                 gAdj=gAdj, verbose=verbose)
  
  # Get point estimates and inference
  R1 <- est$R1
  R0 <- est$R0
  
  if( goal=='aRR' ){
    psi.hat <- log(R1/R0)
  } else if (goal=='RD'){
    psi.hat <- R1- R0
  } else if (goal=='OR'){
    psi.hat <- log( R1/(1-R1)*(1-R0)/R0)
  }
  
  n.clust <- length(unique(data.input$id)) 
  n.pair <- length(unique(data.input$pair))
  
  # INFERENCE FOR THE TREATMENT-SPEC MEANS 
  # note: only implemented with standard (non-CV) inference
  Txt <- get.inference(psi.hat=R1, se=sqrt(est$var.R1), df=(n.clust-2))[,c('est','CI.lo','CI.hi','se')]
  Con <- get.inference(psi.hat=R0, se=sqrt(est$var.R0), df=(n.clust-2))[,c('est','CI.lo','CI.hi','se')]
  
  # NOW FOR THE INTERVENTION EFFECT
  if(break.match){
    # if breaking the match, set df to (#clusters -2)
    df <-n.clust - 2
    var.hat <- est$var.break
  } else{
    # if preserving the match, set df to (#pairs-1)
    df <- n.pair-1 
    var.hat <- est$var.pair
  }
  
  inference <- get.inference(goal=goal, psi=psi, psi.hat=psi.hat, se=sqrt(var.hat), df=df )

  if(do.cv.variance){
    # also return cross-validated inferencec for the intervention effect
    inference.CV <- get.inference(goal=goal, psi=psi, psi.hat=psi.hat, se=sqrt(select$var.CV), df=df )
    
    RETURN<-  data.frame(Txt=Txt, Con=Con, inference, inference.CV, 
                         QAdj=est$QAdj, gAdj=est$gAdj)  
  } else{
    RETURN<-  data.frame(Txt=Txt, Con=Con, inference, 
                         QAdj=est$QAdj, gAdj=est$gAdj)
  }
  
  RETURN
}

#-----------------------------------------------------#-----------------------------------------------------
# do.TMLE: master function to run TMLE 
#
# input: 
#		goal - aRR: arithmetic risk ratio; O/W risk difference,
#	  target; "clust" for cluster-level effect; otherwise indv-level effect
#   training data (train),
#		prespecified adjustment variables for the conditional mean outcome (QAdj), 
#		prespecified adjustment variables for the propensity score (gAdj),
#		initial estimator of the conditional mean outcome (Q.out),
#		estimator of the propensity score (p.out),
#		indicator to print updates (verbose)
#
# output: list 
#		training data augmented with estimates,
#		prespecified adjustment variables for the conditional mean outcome (QAdj), 
#		initial estimator of the conditional mean outcome (Q.out),
#		prespecified adjustment variables for the propensity score (gAdj),
#		estimator of the propensity score (p.out),
#		estimated fluctuation coef (eps),
#		estimated treatment-specific means(R1, R0),
#   estimated variance for each treatment-specifc mean (var.R1, var.R0)
#		estimated variance preserving the match (var.pair),
# 	estimated variance breaking the match (var.break)
#-----------------------------------------------------#-----------------------------------------------------

do.TMLE <- function(goal, target, train, QAdj, gAdj=NULL, Q.out=NULL, p.out=NULL,
                    verbose=F){	
  
  J <- length(unique(train$id))
  
  #=====================================================
  # Step1 - initial estimation of E(Y|A,W)= Qbar(A,W)
  #=====================================================
  
  # run glm on the adjustment set
  Q <- do.Init.Qbar(train=train, QAdj=QAdj, glm.out=Q.out , verbose=verbose)
  train <- Q$train
  
  #==========================================================
  # Step2: Calculate the clever covariate
  #==========================================================	
  
  G <- get.clever.cov(train=train, gAdj=gAdj, p.out=p.out,  verbose=verbose)
  train <- G$train
  
  #==========================================================
  # Step3: Targeting
  #==========================================================
  
  eps <- get.epsilon(train=train, goal=goal, verbose=verbose)
  
  train <- do.targeting(train=train, eps=eps, goal=goal)
  
  #==========================================================
  # Step4: Parameter estimation
  #==========================================================
  
  # IF DATA ARE AT THE CLUSTER-LEVEL: the weights are 1 for cluster-effect
  #   & n_j*J/nTot for individual-target effect 
  # note weights are normalized to sum to J	
  if(nrow(train)> J & target=='clust')	{
    # IF DATA ARE AT THE INDV-LEVEL, but the goal is cluster-level effect
    # Aggregate now 
    train <- aggregate(train, by=list(train$id), mean)[,-1]
    train$alpha <- 1 
  }

  R1 <- mean( train$alpha*train$Qbar1W.star )
  R0 <- mean( train$alpha*train$Qbar0W.star ) 
  
  variance.out <- get.IC.variance(goal=goal, Vdata=train, R1=R1, R0=R0)
  
  RETURN<- list(train=train, 	
                QAdj=Q$QAdj, Q.out=Q$glm.out,
                gAdj=G$gAdj, p.out=G$p.out, 
                eps=eps, 	R1=R1, R0=R0, 
                var.R1=variance.out$var.R1, 
                var.R0=variance.out$var.R0,
                var.pair=variance.out$var.pair, 
                var.break=variance.out$var.break)	
  RETURN
}

#-----------------------------------------------------#-----------------------------------------------------
# do.Init.Qbar - function to do initial estimation of E[Y|A,W] = Qbar(A,W)
# 	input: data set, adjustment variable(s), outcome regression fit, verbose
# 	output:	adjustment variable(s),	outcome regression fit 
#		  data set augmented w/ initial predictions: Qbar(A,W), Qbar(1,W) and Qbar(0,W)
#-----------------------------------------------------#-----------------------------------------------------
do.Init.Qbar<- function(train, QAdj, glm.out=NULL, verbose=F){
  
  if( is.null(QAdj) ){
    QAdj<- 'U'
  }
  
  train.temp <- train[, c(QAdj, 'A', 'Y')]
  
  X1<- X0<- train.temp
  X1$A<-1; X0$A<- 0	
  
  if( is.null(glm.out) ){
    # fit using the training data
    glm.out<- suppressWarnings( glm( Y~. , family='binomial', data=train.temp, weights=train$alpha ) )	
    if(verbose) print(glm.out)
  }	
  
  # get initial predictions
  QbarAW <- predict(glm.out, newdata=train.temp, type='response')
  Qbar1W <- predict(glm.out, newdata=X1, type='response')
  Qbar0W <- predict(glm.out, newdata=X0, type='response')
  
  list(QAdj=QAdj, glm.out=glm.out, train=data.frame(train, QbarAW , Qbar1W , Qbar0W ))
}

#-----------------------------------------------------#-----------------------------------------------------
# get.clever.cov - function calculate the clever covariate
# 	input: 
#		  data set, adjustment variable(s), pscore regression, verbose
# 	output: 
#		  adjustment variable(s), pscore regression, 
#		  data set augmented with pscore & clever covariate (H.AW, H.1W, H.0W)
#-----------------------------------------------------#-----------------------------------------------------
get.clever.cov<- function(train, gAdj, p.out=NULL, verbose=F){
  
  if( is.null(gAdj) ){
    gAdj <- 'U'
  }
  
  train.temp <- train[, c(gAdj, 'A')]  
  
  if( is.null(p.out) ){
    
    # fit pscore on training set 	
    p.out<-   suppressWarnings( glm( A~. , family='binomial', data= train.temp, weights=train$alpha) )
    if(verbose){ print(p.out)}
  }
  
  # now use p.out to get estimated pscores
  pscore <- predict(p.out, newdata= train.temp,  type="response")
  # Note if gAdj=U & train$alpha!=1, then this will differ from 0.5
  
  # bound g - should not apply for a randomized trial
  pscore [pscore < 0.025] <- 0.025
  pscore [pscore > 0.975] <- 0.975
  
  A.train <- train$A
  # Clever covariate is two-dimensional; 
  H.1W <- A.train/pscore 
  H.0W <- (1-A.train)/(1-pscore )
  # via delta method
  H.AW <- H.1W - H.0W
  
  
  list(gAdj=gAdj, p.out=p.out,  train=data.frame(train, pscore, H.1W , H.0W , H.AW) ) 
}	


#-----------------------------------------------------#-----------------------------------------------------
# get.epsilon - function calculate the fluctuation coefficient
# 	input: 
#		  data set, goal with 'aRR'=arithmetic RR, verbose
# 	output: 
#	  	estimated fluctuation coefficient (eps)
#-----------------------------------------------------#-----------------------------------------------------
get.epsilon <- function(train, goal, verbose=F){
  
  A.train<- train$A
  Y.train<- train$Y
  
  # Skip fitting if outcome=0 or outcome=1 for all observations in either txt or control group
  Skip.update <-  mean(Y.train[A.train==1])==0 | mean(Y.train[A.train==0])==0 |  
    mean(Y.train[A.train==1])==1 | mean(Y.train[A.train==0])==1 
  
  if(goal=='RD'){
    
    # if going after RD, then use a 1-dim clever covariate
    if(!Skip.update){
      logitUpdate<- suppressWarnings( 
        glm(Y.train ~ -1 +offset(qlogis(train$QbarAW )) + train$H.AW, family="binomial",  weights=train$alpha))
      eps<-logitUpdate$coef
    } else{
      eps<- 0
    }
    names(eps) <- 'H.AW'
    
  }else{
    # targeting the risk or odds ratio requires a two-dimensional clever covariate
    
    if( !Skip.update  ){
      logitUpdate<- suppressWarnings(
        glm(Y.train ~ -1 +offset(qlogis(train$QbarAW )) + train$H.0W + train$H.1W, family="binomial", weights=train$alpha))
      eps<-logitUpdate$coef
    } else{
      eps <- c(0,0)
    }
    names(eps)<- c('H.0W', 'H.1W')	
  }
  if(verbose) print(eps)
  
  eps
}

#-----------------------------------------------------#-----------------------------------------------------
# do.targeting - function to update initial estimators of QbarAW
# 	input: 
#		  data set (train), fluctuation coefficient (eps), goal (aRR= arithmetic risk ratio; otherwise RD)
# 	output: data.frame w/ targeted predictions: Qbar*(A,W), Qbar*(1,W), Qbar*(0,W)
#-----------------------------------------------------#-----------------------------------------------------

do.targeting <- function(train, eps, goal){
  
  g1W<- train$pscore
  g0W<- (1 - g1W)
  
  if(goal=='RD'){
    
    # updated QbarAW estimates for training set. 
    QbarAW.star <- plogis( qlogis(train$QbarAW ) + eps*train$H.AW)	
    Qbar1W.star <- plogis( qlogis(train$Qbar1W ) + eps/g1W )
    Qbar0W.star <- plogis( qlogis(train$Qbar0W ) - eps/g0W )
    
  }else{
    # updated QbarAW estimates for training set. 
    QbarAW.star <- plogis( qlogis(train$QbarAW) + eps['H.0W']*train$H.0W + eps['H.1W']*train$H.1W)	
    Qbar0W.star <- plogis( qlogis(train$Qbar0W) + eps['H.0W']/g0W )
    Qbar1W.star <- plogis( qlogis(train$Qbar1W) + eps['H.1W']/g1W )
  }
  train <- data.frame(train, QbarAW.star, Qbar1W.star, Qbar0W.star)		
  train
}


#-----------------------------------------------------#-----------------------------------------------------
# get.IC.variance - function to do influence curve-based variance estimate 
# 	input: 
#		  goal (aRR= arithmetic risk ratio; otherwise RD)
#		  dataset (Vdata)
#		  risk estimates under txt and control: R1 & R0
# 	output: 
#		  estimated IC & variance - preserving/breaking the match
#     on log scale for if goal!='RD'
#-----------------------------------------------------#-----------------------------------------------------
get.IC.variance <- function(goal, Vdata, R1, R0){
  
  # CHANGE from SEARCH CODE to population-level effect 
  DY1 <- Vdata$alpha*(Vdata$H.1W*(Vdata$Y - Vdata$Qbar1W.star) + Vdata$Qbar1W.star - R1)
  DY0 <- Vdata$alpha*(Vdata$H.0W*(Vdata$Y - Vdata$Qbar0W.star) + Vdata$Qbar0W.star - R0)	
  
  if( length(DY1)>length(unique(Vdata$id) ) ) {
    # CHANGE from SEARCH CODE - if individual-level data, then aggregate to the cluster-level 
    # only should kick-in if the goal is the individual effect and the data are at the indv level
    DY1 <- c(ltmle:::HouseholdIC(as.matrix(DY1), id = Vdata$id))
    DY0 <- c(ltmle:::HouseholdIC(as.matrix(DY0), id = Vdata$id))
    
    # for the pair-matched IC also need to aggregate to the cluster-level
    Vdata <- aggregate(Vdata, by=list(Vdata$id), mean)[,-1]
  }
  
  if(goal=='RD'){
    # going after RD, easy IC
    DY <-  DY1 - DY0
  } else if (goal=='aRR'){ 
    
    # going after aRR, then get IC estimate on log scale
    #	i.e. Delta method for log(aRR) = log(R1) - log(R0)
    DY <- 1/R1*DY1 - 1/R0*DY0
  } else if(goal=='OR'){
    # Delta method for log(OR)
    DY <- 1/R1*DY1 + 1/(1-R1)*DY1 - 1/(1-R0)*DY0 - 1/R0*DY0
  }
  
  # print(mean(DY))
  
  # estimated variance for txt specific means or if break the match	
  J <- length( unique(Vdata$id) )
  var.R1 <- var(DY1) /J
  var.R0 <- var(DY0) / J
  var.break <- var(DY) /J
  
  
  # estimated variance if preserve the match
  pairs <- unique(Vdata$pair)
  n.pairs <- length(pairs)
  DY.paired <-  rep(NA, n.pairs)
  for(i in 1:n.pairs){		
    these<- Vdata$pair== pairs[i] 
    DY.paired[i]<- 0.5*sum(DY[ these] )			
  }
  
  var.pair <- var(DY.paired) / n.pairs
  
  
  list(var.R1=var.R1, var.R0=var.R0, DY=DY, var.break=var.break, 
       DY.paired=DY.paired, var.pair=var.pair)
}





#-----------------------------------------------------#-----------------------------------------------------
# get.inference: function to calculate (1-sig.level)% confidence intervals & test the null hypothesis
#	input: 
#		goal (aRR= arithmetic risk ratio; otherwise RD)
#   psi (true value)
#   psi.hat (estimate)
#   se (standard error)
#		df (degrees of freedom if using a Student's t-dist ) 
#		sig.level (significance level)
# output: 
#		  point estimate, 95%CI, std error estimate, pval
#-----------------------------------------------------#-----------------------------------------------------	

get.inference <- function(goal='RD', psi=NA, psi.hat, se, df=NA, sig.level=0.05){
  
  # test statistic (on the transformed scale if goal= aRR or OR )
  tstat <- psi.hat/se
  
  if(is.na(df)){
    # ass normal
    cutoff <- qnorm(sig.level/2, lower.tail=F)
    pval<- 2*pnorm(abs(tstat), lower.tail=F) 
  }else{
    # cutoff based on t-dist for testing and CI	
    cutoff <- qt(sig.level/2, df=df, lower.tail=F)
    pval<- 2*pt(abs(tstat), df=df, lower.tail=F)
  }
  
  # 95% confidence interval 
  CI.lo <- (psi.hat - cutoff*se)
  CI.hi <- (psi.hat + cutoff*se)
  
  # transform back 
  if(goal!='RD'){
    psi.hat<- exp(psi.hat)
    CI.lo <- exp(CI.lo)
    CI.hi <- exp(CI.hi)
  }  
  
  # confidence interval coverage
  cover<- ( CI.lo <= psi & psi <= CI.hi )
  # reject the null
  reject <- pval < sig.level 
  
  data.frame(est=psi.hat,  CI.lo, CI.hi, se=se,  pval )
  
}



