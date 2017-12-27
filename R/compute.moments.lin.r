#' A function to compute moments with a linear transect
#' @param ab vector of abundances
#' @param cv vector of coverages
#' @param ht vector of heights
#' @param q vector of orders

compute.moments.lin <- function(ab, q)
{
  ab <- as.matrix(ab)
  n <- dim(ab)[1]  # nr of samples
  nn <- sum(1:n)  # nr of cells
  
  m <- qD <- matrix(0, nrow = nn, ncol = length(q))
  a <- rep(0, nn)
  H <- rep(0, nn)
  pmin <- rep(0, nn)
  
  counter <- 1
  for (ii in 1:(n-1))
  {
    for (jj in 1:(n-ii+1))
    {
      aa <- ab[jj:(jj+ii-1), ]

      if (is.matrix(aa))
      {
        aa <- colSums(aa)
      }
      if (sum(aa) > 0)
      {
        p <- aa/sum(aa)
        
        m[counter,] <- mom(p,q)
        a[counter] <- ii
        H[counter] <- shannon(p)
        pmin[counter] <- min(p[p>0])
        
        counter <- counter + 1
      }
    }
  }
  
  aa <- colSums(ab)
  p <- aa/sum(aa)
  
  m[counter,] <- mom(p,q)
  a[counter] <- n
  H[counter] <- shannon(p)
  pmin[counter] <- min(p[p>0])
  
  qD <- m ^ (1/(1-matrix(q, nrow = dim(m)[1], ncol = dim(m)[2], byrow = T)))
  qD[, q == 1] <- exp(H)
  
  return(list(mom = m[1:counter,], H = H[1:counter], qD = qD[1:counter,], 
              a = a[1:counter], q = q, pmin = pmin[1:counter]))
}