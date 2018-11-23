plot.mom <- function(mom, q = 0, Mq = F, spar = 0.5)
{
  fA <- log(max(mom$a))
  sA <- log(min(mom$a))
  dA <- fA - sA
  sc <- seq(sA + 0.05*dA, fA - 0.05*dA, length = 4)
  
  ind <- mom$q == q
  
  xx <- log(mom$a)
  yy <- if (Mq) log(mom$mom[, ind]) else log(mom$qD[, ind])
  
  lab <- if (Mq) bquote(M[.(q)]) else bquote(''^ ~ .(q) ~ D)
  plot(exp(xx), exp(yy), pch = 19, col = rgb(.7, .9, .7, .5), log = "xy",
       xlab = "L", ylab = lab)
  
  spl <- smooth.spline(x = xx, y = yy, w = xx - min(xx) + 1, spar = spar, tol = .0001)
  
  xxx <- seq(min(xx), max(xx), length = 100)
  
  lines(exp(xxx), exp(predict(spl, xxx)$y), lwd = 2, col = "red")
  
  abline(v = exp(sc), col = c("grey", "red", "blue", "black"))
}

plot.pp <- function(mf)
{
  plot(mf$q, mf$pp, xlab = "q", ylab = "p")
  abline(h = 0.05, lty = 2)
  points(mf$q[mf$pp < 0.05], mf$pp[mf$pp < 0.05], pch = 19, col = "red")
}

plot.delta <- function(mf)
{
  plot(mf$q, mf$delta, xlab = "q", ylab = expression(delta[AIC]))
  abline(h = 0, lty = 2)
  points(mf$q[mf$delta < 0], mf$delta[mf$delta < 0], pch = 19, col = "red")
}

plot.spectra <- function(mfl, xlim = NA, ylim = NA, title = "")
{
  if (sum(is.na(xlim)) > 0) xlim <- range(0, mfl$alfa, 1)
  if (sum(is.na(ylim)) > 0) ylim <- range(0, mfl$f)
  
  plot(mfl$alfa[,1], mfl$f[,1], type = "o", pch = 21, bg = "white", 
       ylim = ylim, xlim = xlim, 
       xlab = "a", ylab = "f", main = title)
  points(mfl$alfa[,2], mfl$f[,2], type = "o", pch = 22, bg = "red")
  points(mfl$alfa[,3], mfl$f[,3], type = "o", pch = 24, bg = "blue")
  points(mfl$alfa[,4], mfl$f[,4], type = "o", pch = 21, bg = "black")
  abline(h = 0, v = 0)
  abline(v = 1, lty = 2)
}