## Part of the following code comes from Lee & Kim (2018) paper
rm(list=ls())
library(gcmr)
library(VineCopula)
library(coda)
library(reshape2)
library(ggplot2)
library(ggpubr)

## Load data---------------------------------------------------------------------------------------------------------------------------------------------------------------------
genes1<- read.delim(file="...Nkx2-1_Sftpa1.txt",
                    header = TRUE, sep = "\t")

## Compute pseudo-observations for copula inference------------------------------------------------------------------------------------------------------------------------------

set.seed(1234)
udat = pobs(genes1)

U <- udat[,1]        # marginal U ~ Uniform (0,1)  
V <- udat[,2]        # marginal V ~ Uniform (0,1)

dat <- data.frame(U,V)



                        ##################################
                       #### ----- Frequentist GCMR ----- ####
                        ##################################



## Regression from V to U-------------------------------------------------------------------------------------
r12<-gcmr( U~V, data = dat, marginal = beta.marg(link = "logit"),
           cormat = arma.cormat(0, 0) )
summary(r12)


Er12<-exp(r12$estimate[1]+dat$V*r12$estimate[2])/
  (1+exp(r12$estimate[1]+dat$V*r12$estimate[2]))

vtou_rho2<-var(Er12)/var(dat$U)  


## Regression from U to V-------------------------------------------------------------------------------------
r21<-gcmr( V~U, data = dat, marginal = beta.marg(link = "logit"),
           cormat = arma.cormat(0, 0) )
summary(r21)


Er21<-exp(r21$estimate[1]+dat$U*r21$estimate[2])/
  (1+exp(r21$estimate[1]+dat$U*r21$estimate[2]))

utov_rho2<-var(Er21)/var(dat$V) 


## Result output is the difference and gives us the direction of influence-------------------------------------------------------------------------------------------------------

rslt <- utov_rho2 - vtou_rho2
rslt  


                     ########################################################################
                     #### -----  Bayesian Parametric Copula Directional Dependence ----- ####
                     ########################################################################



---------------------------------------------- ## Compute the posterior densities ##---------------------------------------------------------------------------------------------

## From U to V ------------------------------------------------------------------------------------------------------------------------------------------------------------------


# This function defines the posterior density of beta0 and beta1 parameters from U to V
#
#
# This function will be used at the Metropolis-Hastings algorithm
#
# u An input vector
# v An input vector
# beta0 parameter of the mean mu of Beta distribution
# beta1 parameter of the mean mu of Beta distribution
# kappa precision parameter of Beta distribution
#
# return posterior samples of beta0 and beta1 for direction U to V
# 
# 
#
lpost_f_UtoV <- function(u, v, beta0, beta1, kappa)
{
  mu <- exp (u *beta1 + beta0) / (1 + exp (u *beta1 + beta0) )
  
  llik <- sum( dbeta (v, mu *kappa, (1-mu) *kappa), log=T)
  
  lprior <- sum (dnorm (c(beta1,beta0), 0, 10), log=T) +
    sum(dgamma (kappa, 1, 1, log=T) )
  
  lpost <- llik + lprior
  
  return(lpost)
}


## Metropolis-Hastings

#  Metropolis- Hastings algorithm to calculate posteriors beta0 and beta1 parameters from U to V
#
#
#  n_iter Number of iterations of the chain
#  beta0_init initial value of the chain
#  beta1_init initial value of the chain
#  kappa value for the precision parameter of Beta distribution
#  u A vector
#  v A vector
#  cand_sd candidate value (proposed) standard deviation
# 
#  return posterior samples of beta0 and {beta1 for direction U to V
# 
#

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



## From V to U-------------------------------------------------------------------------------------------------------------------------------------------------------------------


# This function defines the posterior density of beta0 and beta1 parameters from V to U
#
#
# This function will be used at the Metropolis-Hastings algorithm
#
# u An input vector
# v An input vector
# beta0 parameter of the mean mu of Beta distribution
# beta1 parameter of the mean mu of Beta distribution
# kappa precision parameter of Beta distribution
#
# return posterior samples of beta0 and beta1 for direction V to U
# 
# 
#

lpost_f_VtoU <- function(u, v, beta0, beta1, kappa)
{
  mu <- exp (v *beta1 + beta0) / (1 + exp (v *beta1 + beta0) )
  
  llik <- sum( dbeta (u, mu *kappa, (1-mu) *kappa), log=T)
  
  lprior <- sum (dnorm (c(beta1,beta0), 0, 10), log=T) +
    sum(dgamma (kappa, 1, 1, log=T) )
  
  lpost <- llik + lprior
  
  return(lpost)
}



## Metropolis-Hastings

#  Metropolis- Hastings algorithm to calculate posteriors beta0 and beta1 parameters from V to U
#
#
#  n_iter Number of iterations of the chain
#  beta0_init initial value of the chain
#  beta1_init initial value of the chain
#  kappa value for the precision parameter of Beta distribution
#  u A vector
#  v A vector
#  cand_sd candidate value (proposed) standard deviation
# 
#  return posterior samples of beta0 and {beta1 for direction V to U
# 
#

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




## Running 10000 simulations-----------------------------------------------------------------------------------------------------------------------------------------------------
  

# U to V
set.seed(1234)
post_UtoV = mh_UtoV (n_iter = 10000, beta0_init = 0, beta1_init = 0, kappa = 1, u=U, v=V , cand_sd = 0.28)
str(post_UtoV)

# V to U
set.seed(1234)
post_VtoU = mh_VtoU (n_iter = 10000, beta0_init = 0, beta1_init = 0, kappa = 1, u=U, v=V , cand_sd = 0.28)
str(post_VtoU)


## The traceplots of the regression coefficients from U to V--------------------------------------------------------------------------------------------------------------------
traceplot(as.mcmc(post_UtoV$beta0), main = "Traceplot of beta0 (U to V)")
traceplot(as.mcmc(post_UtoV$beta1),main = "Traceplot of beta1 (U to V)")


## The traceplots of the regression coefficients from V to U---------------------------------------------------------------------------------------------------------------------
traceplot(as.mcmc(post_VtoU$beta0),main = "Traceplot of beta0 (V to U)")
traceplot(as.mcmc(post_VtoU$beta1),main = "Traceplot of beta1 (V to U)")


## Histograms of coefficients----------------------------------------------------------------------------------------------------------------------------------------------------

# U to V
par(mfrow= c(2,2))
hist(as.mcmc(post_UtoV$beta0), freq = FALSE, main = "Histogram of beta0 (U to V)",
     xlab = "beta0", col="dodgerblue")
points(r21$estimate[1], 0.0, pch = 19, col = "indianred")
lines(density(as.mcmc(post_UtoV$beta0)), lty = 5)

hist(as.mcmc(post_UtoV$beta1), freq = FALSE, main = "Histogram of beta1 (U to V)", 
     xlab = "beta1",col = "dodgerblue")
points(r21$estimate[2], 0.0, pch = 19, col = "indianred")
lines(density(as.mcmc(post_UtoV$beta1)), lty = 5)

# V to U
hist(as.mcmc(post_VtoU$beta0), freq = FALSE, main = "Histogram of beta0 (V to U)",
     xlab = "beta0", col="dodgerblue")
points(r12$estimate[1], 0.0, pch = 19, col = "indianred")
lines(density(as.mcmc(post_VtoU$beta0)), lty = 5)

hist(as.mcmc(post_VtoU$beta1), freq = FALSE, main = "Histogram of beta1 (V to U)", 
     xlab = "beta1",col = "dodgerblue")
points(r12$estimate[2], 0.0, pch = 19, col = "indianred")
lines(density(as.mcmc(post_VtoU$beta1)), lty = 5)



## Comparison of the posterior densities of regression coefficients--------------------------------------------------------------------------------------------------------------


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




## Estimations (mean, sd, quantiles)---------------------------------------------------------------------------------------------------------------------------------------------


## U to V coefficients

post_UtoV$beta0_keep = post_UtoV$beta0[-c(1:1000)] # discard early iterations
post_UtoV$beta1_keep = post_UtoV$beta1[-c(1:1000)] # discard early iterations
summary(as.mcmc(post_UtoV$beta0_keep))
summary(as.mcmc(post_UtoV$beta1_keep))


## V to U coefficients

post_VtoU$beta0_keep = post_VtoU$beta0[-c(1:1000)] # discard early iterations
post_VtoU$beta1_keep = post_VtoU$beta1[-c(1:1000)] # discard early iterations
summary(as.mcmc(post_VtoU$beta0_keep))
summary(as.mcmc(post_VtoU$beta1_keep))




## Compute Directionality--------------------------------------------------------------------------------------------------------------------------------------------------------



# U to V
for (i in 1:9000) {
  Er21<-exp(post_UtoV$beta0_keep[i]+dat$U*post_UtoV$beta1_keep[i])/
    (1+exp(post_UtoV$beta0_keep[i]+dat$U*post_UtoV$beta1_keep[i]))
  
  utov_rho2[i] <-var(Er21)/var(dat$V)
}


# V to U
for (i in 1:9000) {
  
  Er12 <-exp(post_VtoU$beta0_keep[i]+dat$V*post_VtoU$beta1_keep[i])/
    (1+exp(post_VtoU$beta0_keep[i]+dat$V*post_VtoU$beta1_keep[i]))
  
  vtou_rho2[i] <-var(Er12)/var(dat$U)
}


# difference
rslt <- utov_rho2 - vtou_rho2


# summarize the results
final2 <-  data.frame(utov_rho2,vtou_rho2)
final_results <- data.frame(utov_rho2,vtou_rho2, rslt)
colnames(final2) <- c("U to V", "V to U")
colnames(final_results) <- c("U to V", "V to U", "difference")



# plot
#data3<- melt(final_results)
data2<- melt(final2)
dens <- ggplot(data2,aes(x=value, fill=variable)) + geom_density(alpha=0.40)+ scale_fill_brewer(palette = "Set1")+ theme_light()
his <- ggplot(data2,aes(x=value, fill=variable)) + geom_histogram(position = "dodge",alpha=0.55)+ scale_fill_brewer(palette = "Set1")+ theme_light()
ggarrange(dens, his,ncol = 1, nrow = 2)



# logical --> check how many samples come from each direction
table(utov_rho2 > vtou_rho2)


# posterior quantiles U ---> V
sintesi_UtoV<-c(quantile(utov_rho2,.05), quantile(utov_rho2,.5),  quantile(utov_rho2,.95))
as.numeric(round(sintesi_UtoV, digits=7))
mean(utov_rho2)

# posterior quantiles V ---> U
sintesi_VtoU<-c(quantile(vtou_rho2,.05), quantile(vtou_rho2,.5),  quantile(vtou_rho2,.95))
as.numeric(round(sintesi_VtoU, digits=7))
mean(vtou_rho2)


# posterior quantiles of their difference
sintesi_diff <- c(quantile(rslt,.05), quantile(rslt,.5),  quantile(rslt,.95))
as.numeric(round(sintesi_diff, digits=7))
mean(sintesi_diff)
