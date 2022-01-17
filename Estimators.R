##############
# Laura B. Balzer
#
# Run t-test, CARE, GEE, and Aug-GEE
# (Functions for running TMLEs are in Stage2_Functions.R & Adapt_Functions.R)

# T-TEST on cluster-level means 
output.ttest<- function(psi, g, gRR=F, paired=F){
  
  if(paired){
    psi.hat<- g$est
  } else{
    psi.hat <- g$est[1]- g$est[2]
  }
  pval <- g$p.value
  reject<- pval< 0.05
  CI.lo <- g$conf.int[1]
  CI.hi <- g$conf.int[2]
  if(gRR){
    psi.hat <- exp(psi.hat)
    CI.lo <-  exp(CI.lo)
    CI.hi <- exp(CI.hi)
  } 
  
  cover<- ( CI.lo <= psi & psi <= CI.hi )
  
  out<- data.frame(est=psi.hat,  CI.lo, CI.hi,  se=g$stderr, pval#,
                   #cover, reject
  )
  out
}

# COVARIATE ADJUSTED RESIDUALS ESTIMATOR (CARE)
do.care <- function(gRR=T, train, ind.cov, psi, paired){
  
  clust <- aggregate(train, by=list(train$id), mean )[,c('id', 'pair', 'A', 'Y')]
  
  # clust.s <- aggregate(train, by=list(train$id), sum )[,c('id', 'pair', 'A', 'Y')]
  
  if(paired){
    n.pairs <- sum(!duplicated(clust$pair)) 
    temp <- data.frame(  train[, c('pair',ind.cov, 'Y') ])
    temp$pair <- as.factor(temp$pair)
  }else{
    temp <- train[, c(ind.cov, 'Y') ]
  }
  # pooled indv regression of outcome on covariates but not txt 
  ind.reg <- glm(Y~  ., data=temp, family='binomial')
  ind.risk <- predict(ind.reg, type='response')
  pred <- aggregate(ind.risk, by=list(train$id), mean)$x
  #pred.s <- aggregate(ind.risk, by=list(train$id), sum)$x
  
  if(gRR){ # relative scale
    resid <- clust$Y/ pred
    # resid.s <- clust.s$Y/pred.s
  }else{ # absolute 
    resid <- clust$Y- pred
  }
  clust <- cbind(clust, pred, resid)
  txt <- clust[clust$A==1, ]
  con <- clust[clust$A==0,]
  
  if(paired){
    red1 <- red0 <- rep(NA, n.pairs)
    temp<- unique(clust$pair)
    for(i in 1:n.pairs){		
      red1[i]<- txt[ txt$pair==temp[i], 'resid'] 
      red0[i]<- con[ con$pair==temp[i], 'resid'] 
    }
  }else{
    red1 <- txt$resid
    red0 <- con$resid
  }
  if(gRR){
    red1 <- log(red1)
    red0 <- log(red0)
  }
  out <- output.ttest(psi=psi, g= t.test(red1, red0, paired=paired, var.equal=T),
                      gRR=gRR, paired=paired)
  out
  
}

# GEE: Generalized estimating equations 
library(geepack)
do.gee <- function(train, ind.cov, psi, paired=F, link){
  
  N <- sum(!duplicated(train$id))
  id <- train$id
  pair <- train$pair
  if(paired){
    train$pair <- as.factor(train$pair)
    train <- train[, c('pair', ind.cov, 'A', 'Y')]
  }else{
    train <- train[, c(ind.cov, 'A', 'Y')]
  }
  m <- geeglm(Y~ ., data=train, family=link, id=id)
  psi.hat<- coef(m)['A']
  se <- summary(m)$coefficients['A', 'Std.err']
  get.inference(goal='aRR', psi=psi, psi.hat=psi.hat, se=se, df=NA)
}


# AUGMENTED GEE
library('CRTgeeDR')
do.gee.aug <- function(train, ind.cov, psi, link){
  
  N <- sum(!duplicated(train$id))
  
  #  get a formula for augmentation terms
  # WARNING: this only works if length(ind.cov) is 2+
  
  temp <- paste(ind.cov[1:(length(ind.cov)-1)], '+ ', sep='', collapse='')
  temp <- paste('Y', '~', temp, ind.cov[length(ind.cov)], sep=' ', collapes='')
  out <- geeDREstimation(formula=Y~A, nameY='Y', id='id', nameTRT='A', 
                         data=train, family=link,
                         model.augmentation.ctrl = as.formula(temp),
                         model.augmentation.trt = as.formula(temp))
  # alternatively can directly specify
                       #  model.augmentation.ctrl = Y~ W1 + W2 + E1 + E2 ,
                        # model.augmentation.trt = Y ~ W1 + W2 + E1 + E2)
  
  psi.hat <- summary(out)$beta[2]
  se <- summary(out)$se.robust[2]    
  get.inference(goal='aRR', psi=psi, psi.hat=psi.hat, se=se, df=NA ) 
  
}


