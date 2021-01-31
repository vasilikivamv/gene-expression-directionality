# The following code is based on the paper from C. Grazian & B. Liseo (2017)

rm(list=ls())

library(ggpubr)
library(gcmr)
library(gmm)
library(copula)
library(dplyr)
library(VineCopula)
library(copula)
library(ggplot2)
library(grid)



## "Quadratic" function--------------------------------------------------------------------------------------
Owen <- function(x, thre, gradiente)
{
  grad <- t(gradiente)
  par <- t(as.matrix(x))
  eps <- 1/thre
  
  z <- 1+par%*%grad
  ans <- z
  lo <- ((z)<eps)
  
  ans[ lo  ] <- log(eps) - 1.5 + 2*(z[lo])/eps - 0.5*(z[lo]/eps)^2
  ans[ !lo ] <- log( z[!lo] )
  -sum(ans)
}



## Computation of empirical likelihood-----------------------------------------------------------------------

##input: 
####### 1) SC --> individual contributions of the estimating function. 
####### It can be either an n x 1 vector (scalar parameter) or an n x p matrix 
####### (vector valued parameter)


##output: a list with elements
######## 1) w0 --> weights associated to each unit
######## 2) conv --> convergence of the algorithm (available only for vector valued parameter)
######## 3) elr --> empirical log-likelihood ratio
######## 4) el --> empirical log-likelihood


EL <- function(SC)
{
  n <- NROW(SC)
  p <- NCOL(SC)
  
  # Lagrange multiplier
  if(p!=1)
  {
    OBJ <- optim(par=rep(0,p), fn=Owen, thre=n, gradiente=SC, control=list(maxit=1e4))
    molt_lagr_0 <- OBJ$pa
    conv <- OBJ$conv
  }
  else
  {
    molt_lagr_0 <- optim(par=rep(0,p), fn=Owen, thre=n, gradiente=SC, method="Brent", lower=-1e1, upper=1e1)$pa
  }
  
  # weights
  w_emp_0 <- as.numeric( 1/(n*(1+molt_lagr_0%*%t(SC))) )
  
  if(p!=1)
  {
    list(w0=w_emp_0, conv=conv, elr=-2*sum(log(w_emp_0))-2*n*log(n), el=log(w_emp_0))
  }
  else
  {
    list(w0=w_emp_0, elr=-2*sum(log(w_emp_0))-2*n*log(n), el=log(w_emp_0) )
  }
}




## Load data--------------------------------------------------------------------------------------------------
genes1<- read.delim(file="D:/Gene_Data/New_data/Nkx2-1_Sftpa1.txt",
                    header = TRUE, sep = "\t")


## Scale the data---------------------------------------------------------------------------------------------

mydata <- data.frame(scale(genes1))



## Compute pseudo-observations for copula inference-----------------------------------------------------------

set.seed(1234)
udat = data.frame(pobs(mydata))

u <- udat[,1]        # marginal U ~ Uniform (0,1)   (??)
v <- udat[,2]        # marginal V ~ Uniform (0,1)

dat <- data.frame(u,v)


## Regression from V to U-------------------------------------------------------------------------------------
r12<-gcmr( u~v, data = dat, marginal = beta.marg(link = "logit"),
           cormat = arma.cormat(0, 0) )
summary(r12)


Er12<-exp(r12$estimate[1]+dat$v*r12$estimate[2])/
  (1+exp(r12$estimate[1]+dat$v*r12$estimate[2]))

vtou_rho2<-var(Er12)/var(dat$u)                  # 1.872321e-06


## Regression from U to V-------------------------------------------------------------------------------------
r21<-gcmr( v~u, data = dat, marginal = beta.marg(link = "logit"),
           cormat = arma.cormat(0, 0) )
summary(r21)


Er21<-exp(r21$estimate[1]+dat$u*r21$estimate[2])/
  (1+exp(r21$estimate[1]+dat$u*r21$estimate[2]))

utov_rho2<-var(Er21)/var(dat$v)                 # 0.000416715



## Bayesian Inference for Copula Directional Dependence-----------------------------------------------------

## Output: a sample of size S of  values approximately from the posterior distribution of rho.

## Direction U to V----------------------------------------------------------------------------------------

ELCOP_UtoV<- function(sample){               # Given observations from unknown copula
  
  S=1000
  
  rh<-runif(S, -1,1)                         # draw rho ~ Unif (-1,1) prior for Spearman's rho
  omega_UtoV<-rep(0,S)
  
  for (s in 1:S) {
    estim_UtoV = utov_rho2 - rh[s]
    
    omega_UtoV[s] <- exp(-EL(estim_UtoV)$elr)
  }
  
  psam_UtoV <- sample(rh, size=S, rep=T, prob=omega_UtoV)        # sample with replacement
  
  sintesi_UtoV<-c(quantile(psam_UtoV,.05), quantile(psam_UtoV,.5),  quantile(psam_UtoV,.95))
  
  return(as.numeric(sintesi_UtoV))
  
}

## Direction V to U------------------------------------------------------------------------------------------
ELCOP_VtoU<- function(sample){              # Given observations from unknown copula
  
  S=1000
  
  rh<-runif(S, -1,1)                        # draw rho ~ Unif (-1,1) prior for Spearman's rho
  omega_VtoU <- rep(0,S)
  for (s in 1:S) {
    estim_VtoU =  vtou_rho2 - rh[s]
    
    omega_VtoU[s] <- exp(-EL(estim_VtoU)$elr)
    
  }
  psam_VtoU<-sample(rh, size=S, rep=T, prob=omega_VtoU)    # sample with replacement
  
  sintesi_VtoU<-c(quantile(psam_VtoU,.05), quantile(psam_VtoU,.5),  quantile(psam_VtoU,.95))
  
  return(as.numeric(sintesi_VtoU)) 
  
}


## Final quantiles of the directional dependence-------------------------------------------------------------

set.seed(1234)
as.numeric(ELCOP_UtoV(dat))                               #  -0.469773943   0.009091824  0.347334895
as.numeric(ELCOP_VtoU(dat))                               #  -0.46310855   -0.01274599   0.41894777



## Plot of the posterior densities of the directional dependence---------------------------------------------
par(mfrow=c(1,1))
plot(density(psam_UtoV),col="indianred",main="Posterior densities of directional dependence Nkx2.1 and Sftpa1")
lines(density(psam_VtoU),col="dodgerblue")



