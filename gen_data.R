## Author: Laura Balzer & Alejandra Benitez 


# get.full.data:  
# input:
#   sim: =1 for sim1 and otherwise sim2
#   J: number of clusters
#   N.mean: average cluster size
#   N.ds: stand dev of cluster sizes
#   effect: indicator of intervention effect (or null)
#   verbose: indicator to prinkt updates
# output: 
#   individual-level data with counterfactual outcomes

get.full.data <- function(sim=2, J=20, N.mean=400, N.sd=200, effect, verbose=F){
  
  n.indv<- round(rnorm(J, N.mean, N.sd))
  n.indv[n.indv < 30] <- 30
  
  full.data<- NULL
  
  if(sim==2){
    for(j in 1:J){
      data.comm.j <- generate.cluster.sim2(effect=effect, N=n.indv[j], j, verbose=verbose) 
      full.data <- rbind(full.data, data.comm.j)
    }
    
  }else{
    for(j in 1:J){
      data.comm.j <- generate.cluster.sim1(effect=effect, N=n.indv[j], j, verbose=verbose) 
      full.data <- rbind(full.data, data.comm.j)
    }
  }
  
  
  X <- aggregate(x = full.data, by = list(id = full.data$id), FUN = mean)
  
  ### PAIR-MATCHED TRIAL
  matchOn<- c('E2')
  dist<- distancematrix(gendistance(data.frame(X[, matchOn])))
  matches<- nonbimatch(dist)
  # matches contains ids for the pair match & the distance measure
  grpA<- as.numeric(matches$halves[,'Group1.Row'])
  grpB<- as.numeric(matches$halves[,'Group2.Row'])
  
  # re-organize the data into paired data & randomize the txt within pairs
  X1<- data.frame(X[grpA, ], Pair=1:(J/2), A= rbinom(J/2, 1, .5))
  X2<- data.frame(X[grpB, ], Pair=1:(J/2), A= ifelse(X1$A==1, 0, 1 ))
  
  X.all<- NULL
  for(i in 1:(J/2)){
    X.all<- rbind(X.all, X1[i,], X2[i,]) 
  }
  #match(X.all$id, full.data$id)
  
  full.data$pair = X.all$Pair[match(full.data$id, X.all$id)]
  full.data$A = X.all$A[match(full.data$id, X.all$id)]
  full.data$Y = ifelse(full.data$A == 1, full.data$Y1, full.data$Y0)
  full.data$alpha <- 1
  full.data
  
}





# generate.cluster.sim1: generate data for simulation 1
# input:
#   effect: indicator of intervention effect (or null)
#   N: cluster size
#   j: cluster number
# output: 
#   individual-level data with counterfactual outcomes
generate.cluster.sim1 <- function(effect=T, N, j, verbose){
  
  UE <- rep(runif(1, -0.2, 1.5), N)
  
  E1 <- rep(rnorm(1, 2, 1), N)
  E2 <- rep(rnorm(1, 0, 1), N)
  
  W1 <- rnorm(N, 2*UE, 0.35)
  W2 <- rnorm(N, 4*UE, 0.9)
  
  UY <- runif(N, 0, 1)
  
  # getY.sim1: generate outcome 
  getY.sim1 <- function(A, W1, W2, E1, E2, UY, verbose=F){

    p <- plogis(-0.75 -0.35*A + 0.8*W1 + 0.4*W2 - 0.3*E1 + 0*E2 - 0.2*A*W2)
    
    if(verbose) { print(summary(p)); hist(p) }
    
    t <- as.numeric(UY <  p)
    t
    
  }
  
  Y0 <- getY.sim1(A=0, W1=W1, W2=W2, E1=E2, E2=E2, UY=UY, verbose=verbose)
  
  if(effect){
    Y1 <- getY.sim1(A=1, W1=W1, W2=W2, E1=E2, E2=E2,  UY=UY, verbose=verbose)
  }else{
    #if under the null, then set the counterfactual outcome Y1 to Y0
    Y1 <- Y0
  }
  if(verbose) print(round(c(mean(Y1), mean(Y0)),2))
  
  
  # return data.frame with id as the cluster indicator and dummy variable U=1 (for unadjusted)
  data.frame( cbind(id=j, n=N, U=1, W1 = W1, W2 = W2, E1 = E1, E2 = E2, Y1, Y0) )
}



# generate.cluster.sim2: generate data for simulation 2 
# input:
#   effect: indicator of intervention effect (or null)
#   N: cluster size
#   j: cluster number
# output: 
#   individual-level data with counterfactual outcomes
generate.cluster.sim2 <- function(effect=T, N, j, verbose){
  
  UE <- runif(1, -1, 1)
  UE2 <- runif(1, -1, 1)
  UE3 <- runif(1, -1, 1) 
  
  E1 <- rep(rnorm(1, 0, 1), N)
  E2 <- rep(rnorm(1, 0, 1), N)
  W1 <- rnorm(N, UE, 0.5)
  W2 <- rnorm(N, UE2, 0.5)
  W3 <- rnorm(N, UE3, 0.5)
  
  UY <- runif(N, 0, 1)

  # getY: generate outcome 
  getY.sim2 <- function(A, W1, W2, W3, E1, E2, N, UY){
    
    N.scale <- N/150
    p <- plogis( .5 + W1/6 + W2/2 + W3/4 + E1/5 + E2/5 - N.scale/8- A*N.scale/5 )
    t <- as.numeric(UY <  p)
    t
    
  }
  
  Y0 <- getY.sim2(A=0, W1=W1, W2=W2, W3=W3, E1=E2, E2=E2, N=N, UY=UY)
  
  if(effect){
    Y1 <- getY.sim2(A=1, W1=W1, W2=W2, W3=W3, E1=E2, E2=E2, N=N, UY=UY)
  }else{
    #if under the null, then set the counterfactual outcome Y1 to Y0
    Y1 <- Y0
  }
  # print(round(c(mean(Y1), mean(Y0)),2))
  # return data.frame with id as the cluster indicator and dummy variable U=1 (for unadjusted)
  data.frame( cbind(id=j, n=N, U=1, W1 = W1, W2 = W2, W3=W3, E1 = E1, E2 = E2, Y1, Y0) )
}



get.truth <- function(X.all){
  X.all.C <- aggregate(X.all, by=list(X.all$id), mean)
  
  truth <- data.frame( cbind(
    int.ind=mean(X.all$Y1),
    con.ind=mean(X.all$Y0),
    RR.ind= mean(X.all$Y1)/mean(X.all$Y0),
    int.clust= mean(X.all.C$Y1),
    con.clust=mean(X.all.C$Y0),
    RR.clust=mean(X.all.C$Y1)/mean(X.all.C$Y0),
    getMeasuresVariability(X1=X.all.C$Y1, X0=X.all.C$Y0, pairs=X.all.C$pair)
  ))
  truth
}


getMeasuresVariability<- function(X1, X0, pairs){	
  
  k.con<- sd(X0)/mean(X0)
  k.txt<- sd(X1)/mean(X1)
  
  temp<- unique(pairs)
  n.pairs<- length(temp)
  varm<- pi <- rep(NA, n.pairs)
  for(i in 1:n.pairs){		
    varm[i]<- var(X0[ pairs==temp[i]])
    pi[i] <- mean(X0[ pairs==temp[i]])
  }
  km <- mean(sqrt(varm)/pi)
  data.frame(k.con, k.txt, km)
}

