
## Marks-Anglin, et al. 2024. Optimal surrogate-assisted sampling for cost-efficient validation of electronic health record outcomes
## simulation code for the novel method OSSAT (Optimal Subsampling strategy with Surrogate-Assisted Two-step procedure)

# X,Y,S
# beta=y~X
# gamma=s~X
# alpha=y~s+X

# yOBS = OSMAC(y)
# sAUG =  proposed (OSSAT)
# sSUB = OSMAC(s)
# yMISS = OSUMC
# Unif = y~X with pilot
# MLE = y~X with all

# rm(list=ls())
library(MASS)
# library(survey)
library(mvtnorm)
# library(psych)
# library(expm)

setwd('/Users/chongliangluo/Library/CloudStorage/Dropbox/R/Sampling/sampling/AA-SiM/rev/')
source("ossat_functions")

nsim = 500
n <- 10000 #full data size
r.size <- seq(800,1600,200)
r1.size <- c(200, 600) #stage 1 sample size
# r2.size <- seq(200,1000,200) #stage 2 sample size
sesp = matrix(c(0.9, 0.7,   # sens1
                0.99, 0.95, # spec1
                0.8, 0.6,   # sens2
                0.95, 0.90),# spec2
              4,2,byrow = T) 

dim.X <- 4 + 4 # number of covariates + intercept
X.dist <- c("mzNormal", "imbNormal", "unNormal", "mixNormal", "T3", "Exp")  

BETA <- array(NA,dim=c(nsim,7,length(r.size),length(X.dist),ncol(sesp),length(r1.size),dim.X)) 
VAR = VAR.MLE <- array(NA,dim=c(nsim,7,length(r.size),length(X.dist),ncol(sesp),length(r1.size),dim.X)) # var(b), var(b-bMLE)
SSP = array(NA,dim=c(n,4,length(X.dist),ncol(sesp),length(r1.size))) 

# ixd = 2
# iss = 2
set.seed(1234)

stage1.weights = rep(NA, n)
par(mfrow=c(2,2))
tt=c
for (c in 3:20){ # rep
  print(c)
  for(iss in 1:2){ # se/sp
    for(ir1 in 1:2){ # r1
    # ir1=1
      sens1 <- sesp[1,iss]
      spec1 <- sesp[2,iss]
      sens2 <- sesp[3,iss]
      spec2 <- sesp[4,iss] 
      r1 = r1.size[ir1] # r1 is the pilot sample size
      for(ixd in 1:length(X.dist)){
        # ixd=2
        beta <- c(0.5,rep(0.5, dim.X-1)) 
        if(X.dist[ixd]=='Exp') beta[1] = -2 # mean(y)=0.45
        X <- X.gen(n, dim.X, dist=X.dist[ixd])
        
        # generate true y and surrogate s
        y <- rbinom(n, 1, 1/(1 + exp(-beta%*%t(X)))) 
        # print(mean(y)) 
        # [1] 0.5652
        # [1] 0.0553
        # [1] 0.6001
        # [1] 0.5389
        # [1] 0.6224
        # [1] 0.4423
        
        # non-differential classification: sens/spec depends on X2
        x2b = X[,2]<= quantile(X[,2],0.3) 
        s <- y
        pr_s = rep(NA, n)
        pr_s[x2b] = sens1*(y[x2b]==1) + (1-spec1)*(y[x2b]==0)
        pr_s[!x2b] = sens2*(y[!x2b]==1) + (1-spec2)*(y[!x2b]==0) 
        s = rbinom(n,1,pr_s)
        # table(y,s) 
        
        ## full sample MLE (benchmark)
        beta.model.y <- glm(y ~ 0+X, family=binomial(link="logit"))
        b.MLE = beta.model.y$coef
        hat.p.y <- 1/(1 + exp(-b.MLE %*%t(X)))
        w.y <- c(hat.p.y * (1-hat.p.y))
        Mx.y <- t(w.y*X) %*% X / n
        
        ## case-control sampling, pilot index # stage1.weights <- rep(1/n,n)
        stage1.weights[s==1] <- 0.5*(1/length(s[s==1]))
        stage1.weights[s==0] <- 0.5*(1/length(s[s==0]))
        set1 <- sample(seq(1,n,1), size=r1,replace=TRUE, prob=stage1.weights) 
        # table(y[set1], s[set1]) 
        while(min(table(y[set1], s[set1]))<3){ # avoid numerical error
          set1 <- sample(seq(1,n,1), size=r1,replace=TRUE, prob=stage1.weights) 
          cat('.')
        }
        
        ## pilot:  y ~ X to obtain approx p^(beta) 
        beta.model <- weighted.model.set1(set1,n,y,X,stage1.weights)
        # fit = logistf(y~., data=data.frame(y,X[,-1])[set1,], weights = 1/stage1.weights[set1],family='binomial')
        hat.p <- 1/(1 + exp(-beta.model$coef %*%t(X)))
        w <- c(hat.p * (1-hat.p))
        Mx <- t(w*X) %*% X / n 
        
        ## full s ~ X to obtain p^*(gamma)
        gamma.model <- glm(s ~ 0+X, family=binomial(link="logit"))  
        hat.p.star <- 1/(1 + exp(-gamma.model$coef %*%t(X)))
        w.hat <- c(hat.p.star*(1-hat.p.star))
        Qx <- t(w.hat*X) %*% X / n   
        
        ## use Bayes rule to calculate p(y=1|s,X), instead of above pilot y~S+X
        # p(y=1|s,X) = p(s|y=1,X)p(y=1|X)/{p(s|y=1,X)p(y=1|X)+p(s|y=0,X)p(y=0|X)}
        hat.p.SX <- rep(NA, n)
        # estimated prevalence (p) ppv npv se sp 
        p = mean(s)
        ppv = mean(y[set1][s[set1]==1]==1)
        npv = mean(y[set1][s[set1]==0]==0)
        se = (1-p-p*npv/(1-npv)) / (p*((1-ppv)/ppv-npv/(1-npv)))
        sp = (p-(1-p)*ppv/(1-ppv)) / ((1-p)*((1-npv)/npv-ppv/(1-ppv)))
        se = max(min(se, 0.99), 0.5) # truncate within 0.5-1 in extreme cases...
        sp = max(min(sp, 0.99), 0.5) 
        hat.p.SX[s==1] = se*hat.p[s==1] / (se*hat.p[s==1] + (1-sp)*(1-hat.p[s==1]))
        hat.p.SX[s==0] = (1-se)*hat.p[s==0] / ((1-se)*hat.p[s==0] + sp*(1-hat.p[s==0])) 
        
        ## sub-sampling prob
        ssp.OSMAC <- abs(c(y-hat.p))*apply(solve(Mx,t(X)), 2, function(a) sqrt(sum(a^2)) )
        ssp.OSMAC <- ssp.OSMAC/sum(ssp.OSMAC)
        ssp.OSMAC.S <- abs(s - hat.p.star) * apply(solve(Qx,t(X)), 2, function(a) sqrt(sum(a^2)) )
        ssp.OSMAC.S <- c(ssp.OSMAC.S / sum(ssp.OSMAC.S) )
        numer01 <- c(hat.p.SX - 2*hat.p.SX*hat.p + hat.p^2) # OSSAT SSP (Proposition 1)
        ssp.OSSAT <- sqrt(numer01)*apply(solve(Mx,t(X)), 2, function(a) sqrt(sum(a^2)))
        ssp.OSSAT <- ssp.OSSAT / sum(ssp.OSSAT) 
        ssp.OSUMC <- sqrt(w) * apply(solve(Mx,t(X)),2, function(a) sqrt(sum(a^2)))
        ssp.OSUMC <- ssp.OSUMC / sum(ssp.OSUMC) 
        
        if(c==1){ # some ssp plots
          SSP[,1:4,ixd,iss,ir1] = cbind(ssp.OSMAC, ssp.OSMAC.S, ssp.OSSAT, ssp.OSUMC)
          print(mean(y))
          print(table(y[set1], s[set1]))
          plot(log10(ssp.OSMAC), log10(ssp.OSMAC.S), main=cor(ssp.OSMAC, ssp.OSMAC.S,method = 'spearman'))
          plot(log10(ssp.OSMAC), log10(ssp.OSSAT), main=cor(ssp.OSMAC, ssp.OSSAT,method = 'spearman'))
        }
        
        cat('<')
        for (k in 1:length(r.size)){ # r is total subsample size
          r <- r.size[k]
          r2 <- r-r1 # r2 is the second step size
          
          #1# benchmark (MLE)
          BETA[c,1,k,ixd,iss,ir1,] <- b.MLE
          
          #2# OSMAC(y) y~X (yOBS)
          model1 <- weighted.model.seq4(set1,n,r2,replace=TRUE,ssp=ssp.OSMAC,y,X,stage1.weights)
          BETA[c,1+1,k,ixd,iss,ir1,] <- model1$coef # - b.MLE
          VAR[c,1+1,k,ixd,iss,ir1,] <- diag(model1$covH)
          VAR.MLE[c,1+1,k,ixd,iss,ir1,] <- diag(model1$cov)
          
          #3# OSMAC(s) S~X (sSUB)
          # OSSAT is optimal in theory, but b/c need pilot est, it may not always be 
          # better than OSMAC(s), which does not need pilot when sens spec high
          model2 <- weighted.model.seq4(set1,n,r2,replace=TRUE,ssp=ssp.OSMAC.S,y,X,stage1.weights)
          BETA[c,2+1,k,ixd,iss,ir1,] <- model2$coef # - b.MLE 
          VAR[c,2+1,k,ixd,iss,ir1,] <- diag(model2$covH)
          VAR.MLE[c,2+1,k,ixd,iss,ir1,] <- diag(model2$cov)
          
          #4# OSSAT proposed: (sAUG)
          model3 <- weighted.model.seq4(set1,n,r2,replace=TRUE,ssp=ssp.OSSAT,y,X,stage1.weights) 
          BETA[c,3+1,k,ixd,iss,ir1,] <- model3$coef # - b.MLE 
          VAR[c,3+1,k,ixd,iss,ir1,] <- diag(model3$covH)
          VAR.MLE[c,3+1,k,ixd,iss,ir1,] <- diag(model3$cov)
          
          #5# OSUMC (yMISS)
          model4 <- weighted.model.seq4(set1,n,r2,replace=TRUE,ssp=ssp.OSUMC,y,X,stage1.weights)
          BETA[c,4+1,k,ixd,iss,ir1,] <- model4$coef # - b.MLE
          VAR[c,4+1,k,ixd,iss,ir1,] <- diag(model4$covH)
          VAR.MLE[c,4+1,k,ixd,iss,ir1,] <- diag(model4$cov)

          #6# Uniform subsampling (Uniform)
          set1u <- sample(seq(1,n,1), size=r,replace=TRUE, prob=rep(1/n,n))
          model5 <- glm(y[set1u] ~ 0+X[set1u,], family=binomial(link="logit"))
          BETA[c,5+1,k,ixd,iss,ir1,] <- model5$coef # - b.MLE 
          VAR[c,5+1,k,ixd,iss,ir1,] <- diag(summary(model5)$cov.unscaled)

          #7# case-control sampling..
          set1cc <- sample(seq(1,n,1), size=r,replace=TRUE, prob=stage1.weights) 
          model6 <- weighted.model.set1(set1cc,n,y,X,stage1.weights)
          BETA[c,6+1,k,ixd,iss,ir1,] <- model6$coef # - b.MLE 
          VAR[c,6+1,k,ixd,iss,ir1,] <- diag(model6$covH)    # Cov(b) by inv Hessian
          VAR.MLE[c,6+1,k,ixd,iss,ir1,] <- diag(model6$cov) # Cov(b - bMLE)
 
        }      
        cat('>')
      }
}
}
}


iss = 1
emp.mse = apply((BETA[,,,,,,-1]-0.5)^2,c(2:6),mean,na.rm=T)  # c=10  ir=5  ixd=6  ir1=2
emp.mse.sd = apply((BETA[,,,,iss,,-1]-0.5)^2,c(2:6),sd)  # c=10  ir=5  ixd=6  ir1=2

# save(list=ls(), file="simu_new_rev3_final.Rdata")
# load('simu_new_rev3_final.Rdata')

# another name for X distributions...
X.dist1 = c('zeroMean','rareEvent','unequalVar','mixNormal','T3','Exp')

pdf("Figure2&3_r1_600.pdf",width=8,height=12)
par(mfrow=c(3,2))
ir1 = 2 # r1=600, Figure 2 and 3
for(ixd in c(1,3,2,4:6)){
  for(iss in 1:2){
    # MLE OSMAC(y) OSMAC(s) OSSAT OSUMC Unif CC 
    plot(emp.mse[2,,ixd,iss,ir1], xlab='r (total subsample size)', xaxt='n', ylab='MSE', 
         ylim=c(min(emp.mse[c(2:7),,ixd,iss,ir1],na.rm=T)*0.9, max(emp.mse[c(2:7),,ixd,iss,ir1],na.rm=T)*1.1),
         type='b', col='grey40',lwd=2, pch=1, 
         main=paste0('X~',X.dist1[ixd], ", sens/spec ", c('high','low')[iss], ', r1=',r1.size[ir1])) 
    lines(emp.mse[3,,ixd,iss,ir1], type='b', col='#E69F00',lwd=2, pch=3)
    lines(emp.mse[4,,ixd,iss,ir1], type='b', col='blue',lwd=2, pch=6)   
    lines(emp.mse[5,,ixd,iss,ir1], type='b', col='grey75',lwd=2, pch=0)
    lines(emp.mse[7,,ixd,iss,ir1], type='b', col='darkgreen',lwd=2, pch=7) 
    legend("topright", bty="n", lty=1, lwd=2, pch=c(1,3,6,0,7),
           col = c('grey40','#E69F00','blue','grey75', 'darkgreen'),
           legend=c('OSMAC(y)','OSMAC(s)','OSSAT','OSUMC','CC') )
    axis(side = 1, at=c(1:length(r.size)), labels = (r.size), las=1)
  }}
dev.off()


pdf("FigureS1&2_r1_200.pdf",width=8,height =12)
par(mfrow=c(3,2))
ir1 = 1 # r1=200, Figure S1 and S2
for(ixd in c(1,3,2,4:6)){
  for(iss in 1:2){
    # MLE OSMAC(y) OSMAC(s) OSSAT OSUMC Unif CC 
    plot(emp.mse[2,,ixd,iss,ir1], xlab='r (total subsample size)', xaxt='n', ylab='MSE', 
         ylim=c(min(emp.mse[c(2:7),,ixd,iss,ir1],na.rm=T)*0.9, max(emp.mse[c(2:7),,ixd,iss,ir1],na.rm=T)*1.1),
         type='b', col='grey40',lwd=2, pch=1, 
         main=paste0('X~',X.dist1[ixd], ", sens/spec ", c('high','low')[iss], ', r1=',r1.size[ir1])) 
    lines(emp.mse[3,,ixd,iss,ir1], type='b', col='#E69F00',lwd=2, pch=3)
    lines(emp.mse[4,,ixd,iss,ir1], type='b', col='blue',lwd=2, pch=6)   
    lines(emp.mse[5,,ixd,iss,ir1], type='b', col='grey75',lwd=2, pch=0)
    lines(emp.mse[7,,ixd,iss,ir1], type='b', col='darkgreen',lwd=2, pch=7) 
    legend("topright", bty="n", lty=1, lwd=2, pch=c(1,3,6,0,7),
           col = c('grey40','#E69F00','blue','grey75', 'darkgreen'),
           legend=c('OSMAC(y)','OSMAC(s)','OSSAT','OSUMC','CC') )
    axis(side = 1, at=c(1:length(r.size)), labels = (r.size), las=1)
  }}
dev.off()



 

## Figure 1: SSP concordance plots 
# SSP[,1:6,ixd,iss,ir1] = cbind(ssp.OSMAC, ssp.OSMAC.S, ssp.OSSAT, ssp.OSUMC, ssp.OSUMC.S, ssp.OSSAT)
pdf("Figure1_SSP_concordance_sesp_low_r1_600.pdf",width=8,height=12)
par(mfrow=c(3,2))
for(ixd in 1:6){
  ssp.OSMAC = SSP[,1,ixd,2,2]
  ssp.OSMAC.S = SSP[,2,ixd,2,2]
  ssp.OSSAT = SSP[,3,ixd,2,2]
  plot(log10(ssp.OSMAC), log10(ssp.OSMAC.S), ylim=range(c(log10(ssp.OSMAC.S),log10(ssp.OSSAT))),
       xlab='log10(SSP), OSMAC(y)',ylab='log10(SSP), OSMAC(s) or OSSAT',
       pch=3, col='#E69F00', cex=0.5, main=paste0('X~',X.dist1[ixd]))
  lines(log10(ssp.OSMAC), log10(ssp.OSSAT), type='p',pch=6, col='blue', cex=0.5)
  abline(0,1)
  legend("bottomright", pch=c(3,6), # bty="n", lty=1, 
         col = c('#E69F00','blue'),
         legend=c(paste0('OSMAC(s), corr=',round(cor(ssp.OSMAC, ssp.OSMAC.S,method = 'spearman'),3)), 
                  paste0('OSSAT, corr=',round(cor(ssp.OSMAC, ssp.OSSAT,method = 'spearman'),3))))
}
dev.off()
 