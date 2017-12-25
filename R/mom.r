#' A function to compute moments
#' @param p abundance vector
#' @param q vector of orders
#' @keywords moment
#' @example mom(c(0.5, 0.3, 0.1, 0.1), 0:5)

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