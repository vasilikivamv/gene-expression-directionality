## Multiple chains-----------------------------------------------------------------------------------------

set.seed(1234)


## U to V


nsim = 10000
post1 = mh_UtoV(n_iter=nsim, beta0_init =1 , beta1_init= 1, cand_sd= 0.22, kappa = 0.25, u=U, v=V)
post1$accpt  #acceptance rate

post2 = mh_UtoV(n_iter=nsim, beta0_init = -1 , beta1_init= -1, cand_sd= 0.22, kappa = 0.25, u=U, v=V)
post2$accpt

post3 = mh_UtoV(n_iter=nsim, beta0_init = 0 , beta1_init= 0, cand_sd= 0.4, kappa = 0.25, u=U, v=V)
post3$accpt


post4 = mh_UtoV(n_iter=nsim, beta0_init = r21$estimate[1] , beta1_init= r21$estimate[2], cand_sd=0.22, kappa = 0.5, u=U, v=V)
post4$accpt

post5 = mh_UtoV(n_iter=nsim, beta0_init = 0.2 , beta1_init= 0.2, cand_sd=0.22, kappa = 0.5, u=U, v=V)
post5$accpt

post6 = mh_UtoV(n_iter=nsim, beta0_init = 3 , beta1_init= 3, cand_sd=0.22, kappa = 0.5, u=U, v=V)
post6$accpt

post7 = mh_UtoV(n_iter=nsim, beta0_init = 5 , beta1_init= 5, cand_sd=0.22, kappa = 0.5, u=U, v=V)
post7$accpt

post8 = mh_UtoV(n_iter=nsim, beta0_init = -3 , beta1_init= -3, cand_sd=0.22, kappa = 0.5, u=U, v=V)
post8$accpt

post9 = mh_UtoV(n_iter=nsim, beta0_init = -5 , beta1_init= -5, cand_sd=0.22, kappa = 0.5, u=U, v=V)
post9$accpt

post10 = mh_UtoV(n_iter=nsim, beta0_init = 0.5 , beta1_init= 0.5, cand_sd=0.22, kappa = 0.5, u=U, v=V)
post10$accpt

## V to U
set.seed(1234)

nsim = 10000
post11 = mh_VtoU(n_iter=nsim, beta0_init =1 , beta1_init= 1, cand_sd= 0.22, kappa = 0.25, u=U, v=V)
post11$accpt  #acceptance rate

post12 = mh_VtoU(n_iter=nsim, beta0_init = -1 , beta1_init= -1, cand_sd= 0.22, kappa = 0.25, u=U, v=V)
post12$accpt

post13 = mh_VtoU(n_iter=nsim, beta0_init = 0 , beta1_init= 0, cand_sd= 0.4, kappa = 0.25, u=U, v=V)
post13$accpt


post14 = mh_VtoU(n_iter=nsim, beta0_init = r21$estimate[1] , beta1_init= r21$estimate[2], cand_sd=0.22, kappa = 0.5, u=U, v=V)
post14$accpt

post15 = mh_VtoU(n_iter=nsim, beta0_init = 0.2 , beta1_init= 0.2, cand_sd=0.22, kappa = 0.5, u=U, v=V)
post15$accpt

post16 = mh_VtoU(n_iter=nsim, beta0_init = 3 , beta1_init= 3, cand_sd=0.22, kappa = 0.5, u=U, v=V)
post16$accpt

post17 = mh_VtoU(n_iter=nsim, beta0_init = 5 , beta1_init= 5, cand_sd=0.22, kappa = 0.5, u=U, v=V)
post17$accpt

post18 = mh_VtoU(n_iter=nsim, beta0_init = -3 , beta1_init= -3, cand_sd=0.22, kappa = 0.5, u=U, v=V)
post18$accpt

post19 = mh_VtoU(n_iter=nsim, beta0_init = -5 , beta1_init= -5, cand_sd=0.22, kappa = 0.5, u=U, v=V)
post19$accpt

post20 = mh_VtoU(n_iter=nsim, beta0_init = 0.5 , beta1_init= 0.5, cand_sd=0.22, kappa = 0.5, u=U, v=V)
post20$accpt


## Traceplots of the direction U to V-----------------------------------------------------------------------

pmc0 = mcmc.list(as.mcmc(post1$beta0), as.mcmc(post2$beta0), 
                 as.mcmc(post3$beta0), as.mcmc(post4$beta0), as.mcmc(post5$beta0),
                 as.mcmc(post6$beta0),
                 as.mcmc(post7$beta0),
                 as.mcmc(post8$beta0),
                 as.mcmc(post9$beta0),
                 as.mcmc(post10$beta0))


coda::traceplot(pmc0, main="Traceplot of beta0 U to V")

pmc1 = mcmc.list(as.mcmc(post1$beta1), as.mcmc(post2$beta1), 
                 as.mcmc(post3$beta1), as.mcmc(post4$beta1), as.mcmc(post5$beta1),
                 as.mcmc(post6$beta1),
                 as.mcmc(post7$beta1),
                 as.mcmc(post8$beta1),
                 as.mcmc(post9$beta1),
                 as.mcmc(post10$beta1))

coda::traceplot(pmc1, main="Traceplot of beta1 U to V")




## Traceplots of the direction V to U-----------------------------------------------------------------------

pmc00 = mcmc.list(as.mcmc(post11$beta0), as.mcmc(post12$beta0), 
                 as.mcmc(post13$beta0), as.mcmc(post14$beta0), as.mcmc(post15$beta0),
                 as.mcmc(post16$beta0),
                 as.mcmc(post17$beta0),
                 as.mcmc(post18$beta0),
                 as.mcmc(post19$beta0),
                 as.mcmc(post20$beta0))


coda::traceplot(pmc00, main="Traceplot of beta0 V to U")

pmc11 = mcmc.list(as.mcmc(post11$beta1), as.mcmc(post12$beta1), 
                  as.mcmc(post13$beta1), as.mcmc(post14$beta1), as.mcmc(post15$beta1),
                  as.mcmc(post16$beta1),
                  as.mcmc(post17$beta1),
                  as.mcmc(post18$beta1),
                  as.mcmc(post19$beta1),
                  as.mcmc(post20$beta1))

coda::traceplot(pmc11, main="Traceplot of beta1 V to U")







# One major difference between the two chains we've looked at is the level of autocorrelation in each.
# Autocorrelation is a number between (-1,1)and which measures how linearly dependent the current value of the
# chain is on past values (called lags). We can see this with an autocorrelation plot:
# 
# Autocorrelation is important because it tells us how much information is available in our Markov chain.
# Sampling 1000 iterations from a highly correlated Markov chain yields less information about
#the stationary distribution than would obtain from 1000 samples
# independently
# drawn from the stationary distribution.



##  U to V------------------------------------------------------------------------------------
coda::autocorr.plot(as.mcmc(post_UtoV$beta0_keep))
coda::autocorr.plot(as.mcmc(post_UtoV$beta1_keep))

coda::autocorr.diag(as.mcmc(post_UtoV$beta0_keep))
coda::autocorr.diag(as.mcmc(post_UtoV$beta1_keep))

## V to U-------------------------------------------------------------------------------------------
coda::autocorr.plot(as.mcmc(post_VtoU$beta0_keep))
coda::autocorr.plot(as.mcmc(post_VtoU$beta1_keep))

coda::autocorr.diag(as.mcmc(post_VtoU$beta0_keep))
coda::autocorr.diag(as.mcmc(post_VtoU$beta1_keep))




# The Monte Carlo effective sample size is how many independent
# samples from the stationary distribution you would have to draw to have equivalent
# information in your Markov chain
coda:: autocorr.plot(as.mcmc(post_UtoV$beta0_keep),lag.max=100)


effectiveSize(as.mcmc(post_UtoV$beta0_keep))




coda::gelman.diag(pmc1)
# Potential scale reduction factors:
#   
#   Point est. Upper C.I.
# [1,]       1.01       1.01

coda::gelman.diag(pmc0)
# Potential scale reduction factors:
#   
#   Point est. Upper C.I.
# [1,]          1       1.01

coda::gelman.diag(pmc11)
# Potential scale reduction factors:
#   
#   Point est. Upper C.I.
# [1,]       1.04       1.08
coda::gelman.diag(pmc00)
# Potential scale reduction factors:
#   
#   Point est. Upper C.I.
# [1,]       1.12       1.24
# 




coda::gelman.plot(pmc1)
coda::gelman.plot(pmc0)
coda::gelman.plot(pmc11)
coda::gelman.plot(pmc00)
