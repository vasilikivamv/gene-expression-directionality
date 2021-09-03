## Load libraries and functions

library(gcmr)
library(VineCopula)
library(coda)
library(reshape2)
library(ggplot2)
library(ggpubr)
source("Bayesian_CDD.R")

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
                        ### ----- Frequentist CDD ----- ##
                        ##################################


# Part of the following code comes from Lee & Kim (2018) paper
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
                     ####      -----  Bayesian Copula Directional Dependence -----       ####
                     ########################################################################



 ## Compute the posterior densities ##---------------------------------------------------------------------------------------------

## From U to V ------------------------------------------------------------------------------------------------------------------------------------------------------------------


post_UtoV <- Bayesian_CDD_UtoV(n_iter = 10000, beta0_init = 0,
                               beta1_init = 0, cand_sd = 0.23, kappa = 0.25, u=U, v=V)



## From V to U-------------------------------------------------------------------------------------------------------------------------------------------------------------------


post_VtoU <- Bayesian_CDD_VtoU(n_iter = 10000, beta0_init = 0,
                               beta1_init = 0, cand_sd = 0.22, kappa = 0.25, u=U, v=V)



## The traceplots of the regression coefficients from U to V--------------------------------------------------------------------------------------------------------------------

traceplot(as.mcmc(post_UtoV$beta0), main = "Traceplot of beta0 (U to V)")
traceplot(as.mcmc(post_UtoV$beta1),main = "Traceplot of beta1 (U to V)")


## The traceplots of the regression coefficients from V to U---------------------------------------------------------------------------------------------------------------------

traceplot(as.mcmc(post_VtoU$beta0),main = "Traceplot of beta0 (V to U)")
traceplot(as.mcmc(post_VtoU$beta1),main = "Traceplot of beta1 (V to U)")



## Estimations of coefficients (mean, sd, quantiles)------------------------------------------------------------------------------------------------------------------------------


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
#data3<- melt(final_results)  uncomment if want to plot the posterior distribution of the difference between the directionality
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
mean(rslt)
