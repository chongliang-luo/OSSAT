## Marks-Anglin, et al. 2024. Optimal surrogate-assisted sampling for cost-efficient validation of electronic health record outcomes
## simulation functions for the novel method OSSAT (Optimal Subsampling strategy with Surrogate-Assisted Two-step procedure)


require(logistf)
# rm(list=ls()) 

## MLE with firth correction for pilot sample estimation, directly use logistf...
getMLE.firth <- function(x, y, w, n) {
  x = data.frame(x)
  w = 1/w # glm model weight is 1/ssp 
  Loop  <- 10
  msg <- "NA"
  
  fit = logistf(y ~ ., firth = T, data=data.frame(y,x[,-1]), weights=w) 
  list(coefficients=fit$coef, message=msg, iter=fit$iter, converged=ifelse(fit$iter<=Loop,T,F), covH=fit$var)
}
 
## generate synthetic data for simulation, dist is X distribution
X.gen <- function(n, dim.X,  
                  dist=c("mzNormal", "nzNormal", "imbNormal", "unNormal", "mixNormal", "T3", "Exp"), lambda=2){
  
  sigma <- matrix(0.5,nr=dim.X-1,nc=dim.X-1)
  diag(sigma) <- rep(1,dim.X-1)
  
  if (dist=="mzNormal"){
    X <-cbind(1, mvrnorm(n, mu=rep(0,dim.X-1), Sigma= sigma))      
  }
  if (dist=="nzNormal"){
    mean.shift = ifelse(dim.X==4, -1, -0.8)
    X <-cbind(1, mvrnorm(n, mu=rep(mean.shift,dim.X-1), Sigma= sigma))     
  }
  if (dist=="imbNormal"){
    mean.shift = ifelse(dim.X==4, -2.5, -1.6)
    X <-cbind(1, mvrnorm(n, mu=rep(mean.shift,dim.X-1), Sigma= sigma))   
  }
  if (dist=="unNormal"){
    d <- seq(1,dim.X-1,1)
    U <- diag(1/d) # diag(sapply(paste("1/",d, sep=""), function(x) eval(parse(text=x))))
    X <- cbind(1, mvrnorm(n, mu=rep(0,dim.X-1), Sigma= U%*%sigma%*%U))      
  }
  if (dist=="mixNormal"){
    mu.shift <- sample(x = c(-1,1), size = n, replace=T, prob=c(0.5,0.5))
    X <- cbind(1, mvrnorm(n, mu=rep(0,dim.X-1),Sigma= sigma) ) 
    X[,-1] <- X[,-1] + mu.shift
    # hist(X[,2])
  }
  if (dist=="T3"){
    X <- cbind(1, rmvt(n, sigma=sigma, df=3)/10)  ## 
  }
  if (dist=="Exp"){
    X <- matrix(NA, nr=n, ncol=dim.X) 
    X[,1] <- 1
    for (p in 2:dim.X){
      X[,p] <- rexp(n, rate=lambda)
      # X[,p] <- scale(X[,p], scale=F)-1/(dim.X-1)
    } 
  }
  
  return(X)
}


## pilot, by default case-control sampling
weighted.model.set1 <- function(set1,n,y,X,stage1.weights){
  #fits weighted model specifically for sequential sampling approach using surrogate outcome
  y.final <- y[c(set1)]
  X.final <- X[c(set1),]
  # stage1.weights = 1/n
  
  # case-control sampling
  if(is.null(stage1.weights)){
    stage1.weights <- rep(NA,n)
    n0 = sum(y==0)
    n1 = sum(y==1)
    stage1.weights[y==0] = 1/n0/2
    stage1.weights[y==1] = 1/n1/2
  }
  model.weights <- stage1.weights[set1]
  
  weighted.model <-getMLE.firth(x=X.final, y=y.final, w=model.weights,n)
    
  hat.p <- 1/(1 + exp(-weighted.model$coefficients %*%t(X.final)))
  w1 <- c(hat.p*(1-hat.p))
  Mx1 <- t(c(1/stage1.weights[set1])*w1*X.final)%*%X.final/(n*length(set1))
  
  psi1 <- t(c(1/stage1.weights[set1])^2*c(y.final-hat.p)^2*X.final)%*%X.final/(n^2*(length(set1))^2) # Vc 
  V1 <- solve(Mx1)%*%psi1%*%solve(Mx1)
  weighted.model$cov <- V1

  return(weighted.model)
}

## this is pilot+stage2 combined SSP weighted est
weighted.model.seq4 <- function(set1,n,r2,replace=TRUE,ssp,y,X,stage1.weights){
  #fits weighted model specifically for sequential sampling approach using surrogate outcome
  #set1 = very first subsample of observations
  #set2 = remaining r2 observations that were sampled sequentially
  #set2 <- sample(seq(1,n,1)[(!seq(1,n,1) %in% set1) & (!is.na(ssp))], r2, replace=TRUE, prob=ssp[(!seq(1,n,1) %in% set1) & (!is.na(ssp))])
  
  model.weights <- rep(0,n)
  # pilot: case-control sampling
  if(is.null(stage1.weights)){
    stage1.weights <- rep(NA,n)
    n0 = sum(y==0)
    n1 = sum(y==1)
    stage1.weights[y==0] = 1/n0/2
    stage1.weights[y==1] = 1/n1/2
  }
  model.weights[set1] <- stage1.weights[set1]
  
  set2 <- sample(seq(1,n,1), r2, replace=TRUE, prob=ssp)
  # model.weights[set2] <- model.weights[set2] + ssp[set2] # add ssps for r2 individuals
  # set12 = unique(c(set1, set2)) # 
  set12 = c(set1, set2)
  n12 = length(set12)
  y.final <- y[c(set12)]
  X.final <- X[c(set12),]
  model.weights <- c(model.weights[set1], ssp[set2])

  weighted.model <-getMLE.firth(x=X.final, y=y.final, w=model.weights,n)

  hat.p <- 1/(1 + exp(-weighted.model$coefficients %*%t(X.final)))
  w1 <- c(hat.p*(1-hat.p))
  Mx1 <- t(c(1/model.weights)*w1*X.final)%*%X.final/(n*n12)
  psi1 <- t(c(1/model.weights)^2*c(y.final-hat.p)^2*X.final)%*%X.final/(n^2*n12^2)
  V1 <- solve(Mx1)%*%psi1%*%solve(Mx1)
  weighted.model$cov = V1
  weighted.model$n12 = length(unique(set12))
  weighted.model$id = unique(set12)
  return(weighted.model)
}
