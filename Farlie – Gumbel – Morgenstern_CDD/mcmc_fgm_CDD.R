# (Preliminary work)
# This code transforms the FGM Copula Directional Dependence from the paper of
# "Kim, J. M., Jung, Y. S., Sungur, E. A., Han, K. H., Park, C., & Sohn, I. (2008). 
# A copula method for modeling directional dependence of genes. BMC bioinformatics, 9(1), 1-12."
# in the Bayesian setting.

# Libraries

library(gcmr)
library(VineCopula)
library(coda)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(truncnorm)


#   Transform data to standard uniform d.f. through its empirical d.f. (Copyright (C) 2018 Namgil Lee & Jong-Min Kim)
 
empiric_df<-function(data,x) {  

  dat<-sort(data)
  
  if(min(dat)>0) a<-0 else a<-floor(min(dat)/100)*100
  if(max(dat)<0) b<-0 else b<-ceiling(max(dat)/100)*100
  
  for(j in 1:length(x))
  {
    if(x[j]<a) x[j]<-a
    if(x[j]>b) x[j]<-b
  }
  
  dat<-c(a,dat,b)
  n<-length(dat)
  p<-c(rep(0,(n-1)))
  q<-c(rep(0,(n-1)))
  
  for(i in 2:(n-2))
  {
    p[i]<-(dat[i]+dat[i+1])/2
    q[i]<-(i-1)/(n-2)
  }
  p[1]<-a
  p[n-1]<-b
  q[1]<-0
  q[n-1]<-1
  approx(p,q,xout=c(x),ties = mean)$y
}



# priors - likelihood - posterior
FGM_copula <- function(u, v, a, b,theta){
  
  llik <- sum (1 + theta*(1-u*(1-a))*((1-u)^(a-1))*(1-v*(1+b))*(1-v)^(b-1), log = T)  # likelihood
  
  
  lprior <-sum (dgamma (c(a,b), 1, 1), log=T) +   # priors
    sum(dunif (theta, -1, 1, log=T) )
  
  lpost <- llik + lprior                # posterior
  
  return(lpost)
} 


# MCMC
metropolis_algorithm <- function(n_iter, a_init, b_init,theta_init, u, v, cand_sd, sd_theta){
  
  accpt = 0
  b_out = numeric(n_iter)              
  a_out = numeric(n_iter)
  theta_out = numeric(n_iter)
  b_now = b_init                  
  a_now = a_init
  theta_now = theta_init
  
  for (i in 1:n_iter) {
    
    a_cand = rtruncnorm (1, a = 1, b= 3,mean = a_now, sd = cand_sd)   # suggest with truncated normal
    b_cand = rtruncnorm (1, a = 1, b=3 , mean = b_now, sd = cand_sd)
    cand_theta =  min(((a_cand+1)/(a_cand-1))^(a_cand-1),((b_cand+1)/(b_cand-1))^(b_cand-1))
    theta_cand = rtruncnorm (1, a = 0, b=cand_theta, mean = theta_now, sd = sd_theta)
    
    
    # calculate the posterior with the new candidates
    lg_now = FGM_copula (a = a_now, b = b_now, theta = theta_now, u = U, v = V)
    lg_cand = FGM_copula (a = a_cand, b = b_cand, theta = theta_cand, u = U, v = V)
    
    
    lalpha = lg_cand - lg_now + dtruncnorm(a_cand) +dtruncnorm(b_cand) +
      dtruncnorm(theta_cand) + dtruncnorm(a_now)+dtruncnorm(b_now)+dtruncnorm(theta_now)
    alpha = exp(lalpha)   
    
    r = runif(1)                    
    
    if (r <  alpha) {             
      a_now = a_cand       
      b_now = b_cand
      theta_now = theta_cand
      lg_now = lg_cand
      accpt = accpt + 1
    }
    
    a_out[i] = a_now
    b_out[i] = b_now
    theta_out[i] = theta_now
    
    
    
  }
  
  list(a = a_out, b = b_out, theta = theta_out, accpt = accpt / n_iter)
  
}


# input data
net1<- read.delim(file="C:/R_Directory/DREAM5_data/Network1/input data/net1_expression_data.tsv",
                  header = TRUE, sep = "\t")
gold<- read.delim(file="C:/R_Directory/DREAM5_data/Network1/gold standard/DREAM5_NetworkInference_GoldStandard_Network1.tsv",
                  header = FALSE, sep = "\t")
interactions <- gold %>% subset(V3==1) %>% select(V1,V2)

# select a gene pair
j=20
genes <- data.frame(net1 %>% select_if(interactions$V1[j]==colnames(net1)),
                    net1 %>% select_if(interactions$V2[j]==colnames(net1)))


n <- nrow(genes)
p <- ncol(genes)

Emp.index <- matrix(rep(0,n*p),n,p)
for(i in 1:p)
  Emp.index[,i] <- empiric_df(genes[,i],genes[,i])

U <- Emp.index[,1]        # marginal U ~ Uniform (0,1)  
V <- Emp.index[,2]        # marginal V ~ Uniform (0,1)
dat <- data.frame(U,V)
u=U
v=V

udat = data.frame(u,v)



# Run the MCMC
set.seed(1234)
post <- metropolis_algorithm(n_iter = 10000, a_init = 1,b_init =1, cand_sd=0.4, sd_theta =0.8,
                             theta_init = 0, u=U, v=V)

str(post)


traceplot(as.mcmc(post$a), main = "Traceplot of a")
traceplot(as.mcmc(post$b), main = "Traceplot of b")
traceplot(as.mcmc(post$theta), main = "Traceplot of theta")



post$a_keep = post$a[-c(1:2000)] # discard early iterations
post$b_keep = post$b[-c(1:2000)] # discard early iterations
post$theta_keep = post$theta[-c(1:2000)]
summary(as.mcmc(post$a_keep))
summary(as.mcmc(post$b_keep))
summary(as.mcmc(post$theta_keep))



# Calculate the directional dependence from the posterior samples
rho.vtou = c()
rho.utov=c()


for (i in 1:length(post$a_keep)) {
  
  rho.utov[i] = (12 * post$theta_keep[i]^2)* ((beta(2,post$b_keep[i]+1))^2)*
    (beta(1,2*post$a_keep[i]-1)-
       2*(1+post$a_keep[i])*beta(2,2*post$a_keep[i]-1)+
       (1+post$a_keep[i])^2*beta(3,2*post$a_keep[i]-1))
  
  
  rho.vtou[i] = (12 * post$theta_keep[i]^2)*((beta(2,post$a_keep[i]+1))^2)*
    (beta(1,2*post$b_keep[i]-1)-
       2*(1+post$b_keep[i])*beta(2,2*post$b_keep[i]-1)+
       (1+post$b_keep[i])^2*beta(3,2*post$b_keep[i]-1))
  
  
  
}


# plot posterior distributions of CDD
final2 <-  data.frame(rho.utov,rho.vtou)
data2<- melt(final2)
dens <- ggplot(data2,aes(x=value, fill=variable)) + geom_density(alpha=0.40)+ scale_fill_brewer(palette = "Set1")+ theme_light()
his <- ggplot(data2,aes(x=value, fill=variable)) + geom_histogram(bins = 500,position = "dodge",alpha=0.55)+ scale_fill_brewer(palette = "Set1")+ theme_light()
ggarrange(dens, his,ncol = 1, nrow = 2)


# logical to find the proposed direction in terms of samples
table(rho.utov >rho.vtou)












### Convergence diagnostics

## Multiple chains-----------------------------------------------------------------------------------------


# Run the MCMC
post1 <- metropolis_algorithm(n_iter = 10000, a_init = 1,b_init =1, cand_sd=0.8, sd_theta = 0.9,
                              theta_init = 0, u=U, v=V)
str(post1)
post2 <- metropolis_algorithm(n_iter = 10000, a_init = 0,b_init =0, cand_sd=0.8, sd_theta = 0.9,
                              theta_init = 0, u=U, v=V)
str(post2)


post3 = metropolis_algorithm(n_iter = 10000, a_init = 10,b_init =10, cand_sd=0.8, sd_theta = 0.9,
                             theta_init = 1, u=U, v=V)
str(post3)

post4 = metropolis_algorithm(n_iter = 10000, a_init = 0,b_init =0, cand_sd=0.8, sd_theta = 0.9,
                             theta_init = -0.5, u=U, v=V)
str(post4)
post5 = metropolis_algorithm(n_iter = 10000, a_init =-2,b_init =-2, cand_sd=0.8, sd_theta = 0.9,
                             theta_init = -1, u=U, v=V)
str(post5)




pmca = mcmc.list(as.mcmc(post1$a), as.mcmc(post2$a), 
                 as.mcmc(post3$a), as.mcmc(post4$a), as.mcmc(post5$a))


pmcb = mcmc.list(as.mcmc(post1$b), as.mcmc(post2$b), 
                 as.mcmc(post3$b), as.mcmc(post4$b), as.mcmc(post5$b))


pmctheta = mcmc.list(as.mcmc(post1$theta), as.mcmc(post2$theta), 
                     as.mcmc(post3$theta), as.mcmc(post4$theta), as.mcmc(post5$theta))


coda::traceplot(pmctheta, main="Traceplot of beta0 U to V")


##  U to V------------------------------------------------------------------------------------
coda::autocorr.plot(as.mcmc(post$a))
coda::autocorr.plot(as.mcmc(post$b))
coda::autocorr.plot(as.mcmc(post$theta))

coda::autocorr.diag(as.mcmc(post$a))
coda::autocorr.diag(as.mcmc(post$theta))


coda:: autocorr.plot(as.mcmc(post$theta),lag.max=100)

# The Monte Carlo effective sample size 
effectiveSize(as.mcmc(post$a_keep))
effectiveSize(as.mcmc(post$b_keep))
effectiveSize(as.mcmc(post$theta_keep))


coda::gelman.diag(pmctheta)
coda::gelman.diag(pmca)
coda::gelman.diag(pmcb)

coda::gelman.plot(pmca)
coda::gelman.plot(pmcb)
coda::gelman.plot(pmctheta)



