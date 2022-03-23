# Data: https://www.synapse.org/#!Synapse:syn2787209/wiki/70350

# Load libraries
library(tidyverse)
library(VineCopula)
library(kdecopula)
library(copula)

# Load data and gold standard
net1 <- read.delim(file=".../DREAM5/Network1/input data/net1_expression_data.tsv",
                 header = TRUE, sep = "\t")
gold <- read.delim(file=".../DREAM5_data/Network1/gold standard/DREAM5_NetworkInference_GoldStandard_Network1.tsv",
                  header = FALSE, sep = "\t")

# data frame with known interactions
interactions <- gold %>% subset(V3==1) %>% select(V1,V2)

# extract a pair of genes (here choose a value for j in ( 1, nrow(interactions) ) )
genes <- data.frame(net1 %>% select_if (interactions$V1[j] == colnames(net1)),
                      net1 %>% select_if (interactions$V2[j] == colnames(net1)))
  
set.seed(1234) # reproducibility

# pseudo-observations to be fed in copula functions
udat <- pobs(genes)  
U <- udat[,1]        # marginal U ~ Uniform (0,1)  
V <- udat[,2]        # marginal V ~ Uniform (0,1)

# fit a nonparametrc copula through kernel estimates
fit <- kdecop(udat)

# store results
VgivenU <- matrix(nrow = nrow(udat), ncol = 1)

# Calculation of copula partial derivatives wrt U=u ---> V|U=u
for (i in 1:nrow(udat)) {
  
  # function to calculate partial derivative
    givenU_fun <- function(v){
      hkdecop(c(udat[i,1],v), fit, cond.var = 1)
    }
  
  # copula regression E[V|U=u]
   VgivenU [i,1] <- 1 - integrate(Vectorize(givenU_fun), lower=0, upper=1)$value
}
  
# store results  
UgivenV <- matrix(nrow = nrow(udat), ncol = 1)
  
# Calculation of partial derivatives wrt V=v ---> U|V=v
for (i in 1:nrow(udat)) {
    
    # function to calculate partial derivative
    givenV_fun <- function(u){
      hkdecop( c(u, udat[i,2]), fit, cond.var = 2)
    }
    
     # copula regression E[U|V=v]
    UgivenV[i,1] <- 1 - integrate(Vectorize(BB1GivenV_fun), lower=0, upper=1)$value
}
  
# U to V
utov <- var( VgivenU[,1] )/ var(V) 
# V to U
vtou <- var( UgivenV[,1] )/ var(U)

# result is their difference
rslt <- utov - vtou
  

## ! If utov > vtou ----> The drectional dependence is from U to V
