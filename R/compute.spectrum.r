#' A function to compute global multifractal spectrum
#' @param mom moments object returned by compute.moments.lin()
#' @param weigh whether to weigh points with scale or not

compute.spectrum <- function(mom, weigh = F)
{
  n <- length(mom$a)
  tau <- pp <- delta <- Dq <- rep(0, length(mom$q))
  xx <- log(mom$a)
  if (weigh) w <- xx - min(xx) + 1 else w <- NULL
  for (jj in 1:length(mom$q))
  {
    yy <- log(mom$qD[,jj])
    fin <- is.finite(yy)
    
    fit.lin <- lm(yy[fin] ~ xx[fin], weights = w)
    fit.quad <- lm(yy[fin] ~ xx[fin] + I(xx[fin]^2), weights = w)
    
    pp[jj] <- anova(fit.quad)[2,5]
    delta[jj] <- AIC(fit.quad) + (2*4*5/(n-4-1)) - AIC(fit.lin) - (2*3*4/(n-3-1))
    
    Dq[jj] <- fit.lin$coefficients[2]
  }
  
  tau <- Dq*(1-mom$q)

  alfa <- tau
  f <- tau
  alfa <- -derv(tau, mom$q[2] - mom$q[1])
  f <- mom$q*alfa + tau
  
  return(list(q = mom$q, tau = tau, Dq = Dq, alfa = alfa, f = f, pp = pp, delta = delta))
}