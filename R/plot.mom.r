plot.mom <- function(mom, q = 0, Mq = F)
{
  fA <- log(max(mom$a))
  sA <- log(min(mom$a))
  dA <- fA - sA
  sc <- seq(sA + 0.05*dA, fA - 0.05*dA, length = 4)
  
  ind <- mom$q == q
  
  xx <- log(mom$a)
  yy <- if (Mq) log(mom$mom[, ind]) else log(mom$qD[, ind])
  
  lab <- if (Mq) bquote(''^ ~ .(q) ~ D) else bquote(''^ ~ .(q) ~ M)
  plot(exp(xx), exp(yy), pch = 19, col = rgb(.7, .9, .7, .5), log = "xy",
       xlab = "L", ylab = lab)
  
  spl <- smooth.spline(x = xx, y = yy, w = xx - min(xx) + 1, spar = 0.3, tol = .0001)
  
  xxx <- seq(min(xx), max(xx), length = 100)
  
  lines(exp(xxx), exp(predict(spl, xxx)$y), lwd = 2, col = "red")
  
  abline(v = exp(sc), col = c("grey", "red", "blue", "black"))
}