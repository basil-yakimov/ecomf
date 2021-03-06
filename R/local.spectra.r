#' A function to compute local multifractal spectra
#' @param inp input object returned by compute.moments.lin()
#' @param sc vector of sample scales
#' @param smooth smoothing parameter

local.spectra <- function(inp, sc = 0, w = 1, smooth = 1)
{
  if (exists("a", inp)) inp$A <- inp$a
  if (length(w) == 1) w <- rep(1, length(inp$a))
  
  q <- inp$q
  
  if (sc[1] == 0)
  {
    fA <- log(max(inp$A))
    sA <- log(min(inp$A))
    dA <- fA - sA
    sc <- seq(sA + 0.05*dA, fA - 0.05*dA, length = 4)
  } else
  {
    sc <- log(sc)
  }
  
  tau <- rep(0,length(q)*length(sc))
  dim(tau) <- c(length(q),length(sc))
  alfa <- f <- Dq <- tau
  
  for (jj in 1:length(q))
  {
    xx <- log(inp$A)
    yy <- log(inp$qD[,jj])
    
    fin <- is.finite(yy)
    
    spl <- smooth.spline(x = xx[fin], y = yy[fin], w = w[fin], spar = smooth, tol = .0001)
    
    der <- predict(spl, x = sc, deriv = 1)
    
    
    Dq[jj,] <- der$y
    tau[jj,] <- der$y*(1-q[jj])
  }

  for (ii in 1:length(sc))
  {
    alfa[,ii] <- -derv(tau[,ii],.1)
    f[,ii] <- q*alfa[,ii] + tau[,ii]
  }

  return(list(a = sc, q = q, tau = tau, Dq = Dq, alfa = alfa, f = f))
}