compute.spectrum <- function(mom)
{
  n <- length(mom$a)
  tau <- pp <- delta <- D <- rep(0,length(mom$q))
  xx <- log(mom$a)
  for (jj in 1:length(mom$q))
  {
    yy <- log(mom$mom[,jj])
    fin <- is.finite(yy)
    
    fit.lin <- lm(yy[fin] ~ xx[fin])
    fit.quad <- lm(yy[fin] ~ xx[fin] + I(xx[fin]^2))
    
    pp[jj] <- anova(fit.quad)[2,5]
    delta[jj] <- AIC(fit.quad) + (2*4*5/(n-4-1)) - AIC(fit.lin) - (2*3*4/(n-3-1))
    
    tau[jj] <- fit.lin$coefficients[2]
  }
  
  D <- tau/(1-mom$q)
  
  if (sum(mom$q == 1) > 0)
  {
    if (exists("H", where = mom))
    {
      yy <- mom$H
      fin <- is.finite(yy)
      fit.lin <- lm(yy[fin] ~ xx[fin])
      D[mom$q == 1] <- fit.lin$coefficients[2]
    } else
    {
      uno <- which(mom$q == 1)
      D[uno] <- (D[uno-1] + D[uno+1])/2
    }
  }
  
  alfa <- tau
  f <- tau
  alfa <- -derv(tau,mom$q[2]-mom$q[1])
  f <- mom$q*alfa + tau
  
  return(list(q = mom$q, tau = tau, D = D, alfa = alfa, f = f, pp = pp, delta = delta))
}