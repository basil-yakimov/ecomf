#' A function to compute Shannon index
#' @param p abundance vector
#' @keywords Shannon
#' @example shannon(c(0.5, 0.3, 0.1, 0.1))

shannon <- function(p)
{
  p <- p[p > 0]
  p <- p/sum(p)
  return(-sum(p*log(p)))
}