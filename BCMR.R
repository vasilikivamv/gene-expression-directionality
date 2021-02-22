## Part of the following code comes from Lee & Kim (2018) paper
rm(list=ls())
library(gcmr)
library(VineCopula)
library(coda)

## Load data--------------------------------------------------------------------------------------------------
genes1<- read.delim(file="D:/Gene_Data/New_data/Nkx2-1_Sftpa1.txt",
                    header = TRUE, sep = "\t")

## Compute pseudo-observations for copula inference-----------------------------------------------------------

set.seed(1234)
udat = data.frame(pobs(genes1))

U <- udat[,1]        # marginal U ~ Uniform (0,1)   (??)
V <- udat[,2]        # marginal V ~ Uniform (0,1)

dat <- data.frame(U,V)



                        ##################################
                        #### ----- Classic GCMR ----- ####
                        ##################################



## Regression from V to U-------------------------------------------------------------------------------------
r12<-gcmr( U~V, data = dat, marginal = beta.marg(link = "logit"),
           cormat = arma.cormat(0, 0) )
summary(r12)


Er12<-exp(r12$estimate[1]+dat$V*r12$estimate[2])/
  (1+exp(r12$estimate[1]+dat$V*r12$estimate[2]))

vtou_rho2<-var(Er12)/var(dat$U)  ##  1.872321e-06


## Regression from U to V-------------------------------------------------------------------------------------
r21<-gcmr( V~U, data = dat, marginal = beta.marg(link = "logit"),
           cormat = arma.cormat(0, 0) )
summary(r21)


Er21<-exp(r21$estimate[1]+dat$U*r21$estimate[2])/
  (1+exp(r21$estimate[1]+dat$U*r21$estimate[2]))

utov_rho2<-var(Er21)/var(dat$V)  ## 0.000416715


## Result output is the difference and gives us the direction of influence------------------------------------

rslt <- utov_rho2 - vtou_rho2
rslt  ## 0.0004148427




                     ##########################################################
                     #### -----  Bayesian Copula Marginal Regression ----- ####
                     ##########################################################

## Compute the posterior densities---------------------------------------------------------------------------------

   ## From U to V ##
lpost_f_UtoV <- function(u, v, beta0, beta1, kappa)
{
  mu <- exp (u *beta1 + beta0) / (1 + exp (u *beta1 + beta0) )
  
  llik <- sum( dbeta (v, mu *kappa, (1-mu) *kappa), log=T)
  
  lprior <- sum (dnorm (c(beta1,beta0), 0, 10), log=T) +
    sum(dgamma (kappa, 1, 1, log=T) )
  
  lpost <- llik + lprior
  
  return(lpost)
}


## From V to U ##
lpost_f_VtoU <- function(u, v, beta0, beta1, kappa)
{
  mu <- exp (v *beta1 + beta0) / (1 + exp (v *beta1 + beta0) )
  
  llik <- sum( dbeta (u, mu *kappa, (1-mu) *kappa), log=T)
  
  lprior <- sum (dnorm (c(beta1,beta0), 0, 10), log=T) +
    sum(dgamma (kappa, 1, 1, log=T) )
  
  lpost <- llik + lprior
  
  return(lpost)
}

mh_UtoV = function(n_iter, beta0_init, beta1_init, kappa, u, v, cand_sd) {
  
  
  beta0_out = numeric(n_iter)              
  beta1_out = numeric(n_iter)             
  accpt = 0                                ## track the accepted values
  beta0_now = beta0_init                   ## initial value for beta0 coefficient
  beta1_now = beta1_init                  ## initial value for beta1 coefficient
  
  ## posterior density now
  lg_now = lpost_f_UtoV (beta0 = beta0_now, beta1 = beta1_now, kappa = kappa, u = U, v = V)
  
  
  
  for (i in 1:n_iter) {
    
    beta0_cand = rnorm (1, mean = beta0_now, sd = cand_sd)        # suggest candidate for beta0
    beta1_cand = rnorm (1, mean = beta1_now, sd = cand_sd)        # suggest candidate for beta1
    
    # calculate the posterior with the new candidates
    lg_cand = lpost_f_UtoV (beta0 = beta0_cand, beta1 = beta1_cand, u = U, v = V, kappa = kappa)
    
    lalpha = lg_cand - lg_now        # calculate alpha ratio in log-form
    alpha = exp(lalpha)              # transform back to the original scale
    
    u = runif(1)                     # generate a uniform number
    
    if (u < alpha) {                  # if u < alpha then
      beta0_now = beta0_cand          # accept the candidate for beta0
      beta1_now = beta1_cand          # accept the candidate for beta1
      accpt = accpt + 1
      lg_now = lg_cand
    }
    
    beta0_out[i] = beta0_now
    beta1_out[i] = beta1_now
    
    
  }
  
  list(beta0 = beta0_out, beta1 = beta1_out, accpt = accpt / n_iter)
}




## Metropolis-Hastings V to U--------------------------------------------------------------------

mh_VtoU = function(n_iter, beta0_init, beta1_init, kappa, u, v, cand_sd) {
  
  
  beta0_out = numeric(n_iter)
  beta1_out = numeric(n_iter)
  accpt = 0
  beta0_now = beta0_init
  beta1_now = beta1_init
  lg_now = lpost_f_VtoU (beta0 = beta0_now, beta1 = beta1_now, kappa = kappa, u = U, v = V)
  
  for (i in 1:n_iter) {
    
    beta0_cand = rnorm(1, mean = beta0_now, sd = cand_sd)
    beta1_cand = rnorm(1, mean = beta1_now, sd = cand_sd)
    lg_cand = lpost_f_VtoU (beta0 = beta0_cand, beta1 = beta1_cand, u = U, v = V, kappa = kappa)
    
    lalpha = lg_cand - lg_now
    alpha = exp(lalpha)
    
    u = runif(1)
    
    if (u < alpha) {
      beta0_now = beta0_cand
      beta1_now = beta1_cand
      accpt = accpt + 1
      lg_now = lg_cand
    }
    
    beta0_out[i] = beta0_now
    beta1_out[i] = beta1_now
    
    
  }
  
  list(beta0 = beta0_out, beta1 = beta1_out, accpt = accpt / n_iter)
}




## After 10000 iterations the results are----------------------------------------------------------------------
  

# U to V
set.seed(1234)
post_UtoV = mh_UtoV (n_iter = 10000, beta0_init = 0, beta1_init = 0, kappa = 1, u=U, v=V , cand_sd = 0.28)
str(post_UtoV)

# V to U
set.seed(1234)
post_VtoU = mh_VtoU (n_iter = 10000, beta0_init = 0, beta1_init = 0, kappa = 1, u=U, v=V , cand_sd = 0.28)
str(post_VtoU)


## The traceplots of the regression coefficients from U to V--------------------------------------------------

traceplot(as.mcmc(post_UtoV$beta0), main = "Traceplot of beta0 (U to V)")
traceplot(as.mcmc(post_UtoV$beta1),main = "Traceplot of beta1 (U to V)")


## The traceplots of the regression coefficients from V to U--------------------------------------------------
  

traceplot(as.mcmc(post_VtoU$beta0),main = "Traceplot of beta0 (V to U)")
traceplot(as.mcmc(post_VtoU$beta1),main = "Traceplot of beta1 (V to U)")


## Histograms------------------------------------------------------------------------------------------------

par(mfrow= c(2,2))
hist(as.mcmc(post_UtoV$beta0), freq = FALSE, main = "Histogram of beta0 (U to V)",
     xlab = "beta0", col="dodgerblue")
points(r21$estimate[1], 0.0, pch = 19, col = "indianred")
lines(density(as.mcmc(post_UtoV$beta0)), lty = 5)

hist(as.mcmc(post_UtoV$beta1), freq = FALSE, main = "Histogram of beta1 (U to V)", 
     xlab = "beta1",col = "dodgerblue")
points(r21$estimate[2], 0.0, pch = 19, col = "indianred")
lines(density(as.mcmc(post_UtoV$beta1)), lty = 5)


hist(as.mcmc(post_VtoU$beta0), freq = FALSE, main = "Histogram of beta0 (V to U)",
     xlab = "beta0", col="dodgerblue")
points(r12$estimate[1], 0.0, pch = 19, col = "indianred")
lines(density(as.mcmc(post_VtoU$beta0)), lty = 5)

hist(as.mcmc(post_VtoU$beta1), freq = FALSE, main = "Histogram of beta1 (V to U)", 
     xlab = "beta1",col = "dodgerblue")
points(r12$estimate[2], 0.0, pch = 19, col = "indianred")
lines(density(as.mcmc(post_VtoU$beta1)), lty = 5)



## Plot comparison for the posterior densities---------------------------------------------------------------


plot(density(post_VtoU$beta0), xlim = c(-2 , 2), 
     xlab = "beta0",main = "Posterior density of beta0 ", col = "dodgerblue")
polygon(density(post_VtoU$beta0), col="dodgerblue", border="black")

lines(density(post_UtoV$beta0), xlim = c(-2 , 2), col = "indianred")
polygon(density(post_UtoV$beta0), col="indianred", border="blue")

grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)


plot(density(post_VtoU$beta1), xlim = c(-2 , 2), 
     xlab = "beta1",main = "Posterior density of beta1", col = "dodgerblue")
polygon(density(post_VtoU$beta1), col="dodgerblue", border="black")

lines(density(post_UtoV$beta1), xlim = c(-2 , 2),  col = "indianred")
polygon(density(post_UtoV$beta1), col="indianred", border="blue")

grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)



## Multiple chains-----------------------------------------------------------------------------------------

set.seed(1234)


## U to V


nsim = 10000
post1 = mh_UtoV(n_iter=nsim, beta0_init = 5, beta1_init= 5, cand_sd= 0.4, kappa = 1, u=U, v=V)
post1$accpt  #acceptance rate

post2 = mh_UtoV(n_iter=nsim, beta0_init = 1, beta1_init= 1, cand_sd= 0.45, kappa = 1, u=U, v=V)
post2$accpt

post3 = mh_UtoV(n_iter=nsim, beta0_init = 0 , beta1_init= 0, cand_sd= 0.5, kappa = 1, u=U, v=V)
post3$accpt


post4 = mh_UtoV(n_iter=nsim, beta0_init = 0, beta1_init=0, cand_sd=0.35, kappa = 1, u=U, v=V)
post4$accpt

post5 = mh_UtoV(n_iter=nsim, beta0_init = 0, beta1_init=0, cand_sd=0.3, kappa = 1, u=U, v=V)
post5$accpt



## V to U
set.seed(1234)

nsim = 10000
post6 = mh_VtoU(n_iter=nsim, beta0_init = 5, beta1_init= 5, cand_sd= 0.4, kappa = 1, u=U, v=V)
post6$accpt

post7 = mh_VtoU(n_iter=nsim, beta0_init = 1, beta1_init= 1, cand_sd= 0.45, kappa = 1, u=U, v=V)
post7$accpt

post8 = mh_VtoU(n_iter=nsim, beta0_init = 0 , beta1_init= 0, cand_sd= 0.5, kappa = 1, u=U, v=V)
post8$accpt


post9 = mh_VtoU(n_iter=nsim, beta0_init = 0, beta1_init=0, cand_sd=0.35, kappa = 1, u=U, v=V)
post9$accpt

post10 = mh_VtoU(n_iter=nsim, beta0_init = 0, beta1_init=0, cand_sd=0.3, kappa = 1, u=U, v=V)
post10$accpt



## Traceplots of the direction U to V-----------------------------------------------------------------------

pmc0 = mcmc.list(as.mcmc(post1$beta0), as.mcmc(post2$beta0), 
                 as.mcmc(post3$beta0), as.mcmc(post4$beta0), as.mcmc(post5$beta0))


coda::traceplot(pmc0, main="Traceplot of beta0")

pmc1 = mcmc.list(as.mcmc(post1$beta1), as.mcmc(post2$beta1), 
                 as.mcmc(post3$beta1), as.mcmc(post4$beta1), as.mcmc(post5$beta1))

coda::traceplot(pmc1, main="Traceplot of beta1")




## Traceplots of the direction V to U-------------------------------------------------------------------------

pmc2 = mcmc.list(as.mcmc(post6$beta0), as.mcmc(post7$beta0), 
                 as.mcmc(post8$beta0), as.mcmc(post9$beta0), as.mcmc(post10$beta0))


coda::traceplot(pmc2, main="Traceplot of beta0")

pmc3 = mcmc.list(as.mcmc(post6$beta1), as.mcmc(post7$beta1), 
                 as.mcmc(post8$beta1), as.mcmc(post9$beta1), as.mcmc(post10$beta1))

coda::traceplot(pmc3, main="Traceplot of beta1")



## Estimations (mean, sd, quantiles)--------------------------------------------------------------------------


## U to V

post_UtoV$beta0_keep = post_UtoV$beta0[-c(1:1000)] # discard early iterations
post_UtoV$beta1_keep = post_UtoV$beta1[-c(1:1000)]
summary(as.mcmc(post_UtoV$beta0_keep))
summary(as.mcmc(post_UtoV$beta1_keep))


## V to U

post_VtoU$beta0_keep = post_VtoU$beta0[-c(1:1000)]
post_VtoU$beta1_keep = post_VtoU$beta1[-c(1:1000)]
summary(as.mcmc(post_VtoU$beta0_keep))
summary(as.mcmc(post_VtoU$beta1_keep))
