mom <- function(p, q)
{
  p <- p[p > 0]
  p <- p/sum(p)
  m <- rep(0, length(q))
  for (ii in 1:length(q))
  {
    m[ii] <- sum(p^q[ii])
  }
  return(m)
}