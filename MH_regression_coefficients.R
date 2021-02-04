## "Quadratic" function------------------------------------------------------------------------------------------
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

## Empirical Likelihood------------------------------------------------------------------------------------------
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



## Metropolis algorithm for the regression coefficients----------------------------------------------------------

easyMCMC = function(niter, startval, proposalsd){
  x = rep(0,niter)
  x[1] = startval     
  for(i in 2:niter){
    currentx = x[i-1]
    proposedx = rnorm(1,mean=currentx,sd=proposalsd)      # proposal distribution Normal
    A =  exp(-(EL(proposedx)$elr)) / exp(-(EL(currentx)$elr))      # target distribution the Empirical Likelihood
    if(runif(1)<A){
      x[i] = proposedx                      # accept move with probabily min(1,A)
    } else {
      x[i] = currentx                       # otherwise "reject" move, and stay where we are
    }
  }
  return(x)
}


## Run 3 times-------------------------------------------------------------------------------------------------------


z1=easyMCMC(10000,-2,0.1)
z2=easyMCMC(10000,-2,0.1)
z3=easyMCMC(10000,-2,0.1)

## Plot the chains------------------------------------------------------------------------------------------------

plot(z1,type="l")
lines(z2,col="dodgerblue")
lines(z3,col="green")
abline(h=0.05774,col="indianred",lwd=1.5)    # true value from GCMR for Intercept ??0 = 0.05774

