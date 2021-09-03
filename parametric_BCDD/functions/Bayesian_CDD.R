
Bayesian_CDD_UtoV <- function(n_iter, beta0_init, beta1_init, kappa,
                         u, v , cand_sd){
  ## Compute the posterior densities 
  
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
  
  set.seed(1234)
  post_UtoV = mh_UtoV (n_iter = n_iter, beta0_init = beta0_init, beta1_init = beta1_init,
                       kappa = kappa, u=u, v=v ,
                       cand_sd = cand_sd)
  str(post_UtoV)
  return(post_UtoV)
  
}




Bayesian_CDD_VtoU <- function(n_iter, beta0_init, beta1_init, kappa,
                              u, v , cand_sd){
  
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
  
  set.seed(1234)
  post_VtoU = mh_VtoU (n_iter = n_iter, beta0_init = beta0_init, beta1_init = beta1_init,
                       kappa = kappa, u=u, v=v ,
                       cand_sd = cand_sd)
  str(post_VtoU)
  return(post_VtoU)
  
}
