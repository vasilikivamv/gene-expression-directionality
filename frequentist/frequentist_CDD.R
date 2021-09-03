## Part of the following code is based on the paper of Lee & Kim (2018)

library(gcmr)
library(VineCopula)

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


## Result output is the difference and gives us the direction of influence------------------------------------

rslt <- utov_rho2 - vtou_rho2
rslt



