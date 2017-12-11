
shannon <- function(p)
{
  p <- p[p > 0]
  p <- p/sum(p)
  return(-sum(p*log(p)))
}