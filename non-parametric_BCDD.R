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
library(reshape2)


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
genes1<- read.delim(file="...Nkx2-1_Sftpa1.txt",
                    header = TRUE, sep = "\t")




## Compute pseudo-observations for copula inference-----------------------------------------------------------

set.seed(1234)
udat = data.frame(pobs(mydata))

u <- udat[,1]        # marginal U ~ Uniform (0,1)   
v <- udat[,2]        # marginal V ~ Uniform (0,1)

dat <- data.frame(u,v)


## Regression from V to U-------------------------------------------------------------------------------------
r12<-gcmr( u~v, data = dat, marginal = beta.marg(link = "logit"),
           cormat = arma.cormat(0, 0) )
summary(r12)


Er12<-exp(r12$estimate[1]+dat$v*r12$estimate[2])/
  (1+exp(r12$estimate[1]+dat$v*r12$estimate[2]))

vtou_rho2<-var(Er12)/var(dat$u)               


## Regression from U to V-------------------------------------------------------------------------------------
r21<-gcmr( v~u, data = dat, marginal = beta.marg(link = "logit"),
           cormat = arma.cormat(0, 0) )
summary(r21)


Er21<-exp(r21$estimate[1]+dat$u*r21$estimate[2])/
  (1+exp(r21$estimate[1]+dat$u*r21$estimate[2]))

utov_rho2<-var(Er21)/var(dat$v)               



## Non-parametric Bayesian Copula Directional Dependence-------------------------------------------------------------------------------------

## Output: a sample of size S of  values approximately from the posterior distribution of the directionality

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
  
  return(as.numeric(psam_UtoV))
  
  
  
  
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
  
  return(as.numeric(psam_VtoU))
  
}


## Final quantiles of the directional dependence-------------------------------------------------------------

set.seed(1234)
UtoV <- as.numeric(ELCOP_UtoV(dat))     # U to V                          
VtoU <- as.numeric(ELCOP_VtoU(dat))     # V to U                          


sintesi_UtoV<-c(quantile(UtoV,.05), quantile(UtoV,.5),  quantile(UtoV,.95))   # quantiles U to V
as.numeric(sintesi_UtoV)
mean(UtoV)  

sintesi_VtoU<-c(quantile(VtoU,.05), quantile((VtoU),.5),  quantile(VtoU,.95))  # quantiles V to U
as.numeric(sintesi_VtoU)
mean(VtoU)


rslt <- UtoV - VtoU     # difference


sintesi_diff <- c(quantile(rslt,.05), quantile(rslt,.5),  quantile(rslt,.95))    # quantiles of the difference
as.numeric(sintesi_diff)
mean(sintesi_diff)


table(UtoV > VtoU)  # logical count of the directions

## Plottinng------------------------------------------------------------------------------------------------------------------

UtoV <- as.data.frame(UtoV)
VtoU <- as.data.frame(VtoU)
rslt <- as.data.frame(rslt)
final3 <-  data.frame(UtoV, VtoU,rslt)
final2 <-  data.frame(UtoV, VtoU)
colnames(final2) <- c("Nkx2.1 to Sftpa1", "Sftpa1 to Nkx2.1")
colnames(final3) <- c("Nkx2.1 to Sftpa1", "Sftpa1 to Nkx2.1", "difference")




data3<- melt(final3)
data2<- melt(final2)
dens <- ggplot(data3,aes(x=value, fill=variable)) + geom_density(alpha=0.40)+ scale_fill_brewer(palette = "Set1")+ theme_light()
his <- ggplot(data2,aes(x=value, fill=variable)) + geom_histogram(position = "dodge",alpha=0.55)+ scale_fill_brewer(palette = "Set1")+ theme_light()

ggarrange(dens, his,ncol = 1, nrow = 2)

