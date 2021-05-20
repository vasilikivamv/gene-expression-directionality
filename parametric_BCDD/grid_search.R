# Grid search--------------------------------------------------------------------------------------------
set.seed(1234)

# Creates a set of combinations of parameters for grid search
hyper_grid <- expand.grid(
  beta0_init = runif(5,-1,1),
  beta1_init = runif(5,-1,1), 
  cand_sd = c(0.08,0.18,0.28,0.38,0.5),                               
  n_iter = 10000,                       
  kappa = runif(10,0,1),
  accept = NA)


for(i in seq_len(nrow(hyper_grid))) {

  post = mh_UtoV(n_iter = hyper_grid$n_iter[i],
                    beta0_init = hyper_grid$beta0_init[i],
                    beta1_init = hyper_grid$beta1_init[i],
                    cand_sd = hyper_grid$cand_sd[i],
                    kappa = hyper_grid$kappa[i],
                    u=U,
                    v=V)
  hyper_grid$accept[i] <- post$accpt
}


library(dplyr)
# the top 100 models arranged in terms of acceptance rate
hyper_grid %>%
  arrange(accept) %>%
  head(100)

sum ( (0.22 <= hyper_grid$accept) & (hyper_grid$accept <= 0.25) ) # 153 model in this acceptance range

# contains the values of the paramters of the best models
best <- subset(hyper_grid, (0.22 <= hyper_grid$accept) & (hyper_grid$accept <= 0.25))

