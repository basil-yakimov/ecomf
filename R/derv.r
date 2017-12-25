#' A function to compute local derivative
#' @param vec vector
#' @param h interval size
#' @keywords derivative
#' @example derv(vec = seq(0,5, by = 0.1)^ 2, h = 0.1)

derv <- function(vec,h)
{
  n <- length(vec)
  f <- rep(0,n)
  f[2:(n-1)] <- (vec[3:n] - vec[1:(n-2)]) / (2*h)
  f[1] <- (-3*vec[1] + 4*vec[2] - vec[3]) / (2*h)
  f[n] <- (vec[n-2] - 4*vec[n-1] + 3*vec[n]) / (2*h)
  return(f)
}
