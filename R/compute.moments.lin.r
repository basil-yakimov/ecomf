#' A function to compute moments with a linear transect
#' @param ab vector of abundances
#' @param q vector of orders
#' @param x vector of sample positions

compute.moments.lin <- function(ab, q, x = 0)
{
  
  ab <- as.matrix(ab)
  n <- dim(ab)[1]  # nr of samples
  nn <- sum(1:n)  # nr of cells
  regular <- F
  if (length(x) != n) {
    x <- 1:n
    regular <- T
  }
  
  m <- qD <- matrix(0, nrow = nn, ncol = length(q))
  a <- rep(0, nn)
  H <- rep(0, nn)
  pmin <- nmin <- rep(0, nn)
  
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
        a[counter] <- x[jj+ii-1] - x[jj]
        H[counter] <- shannon(p)
        pmin[counter] <- min(p[p > 0])
        nmin[counter] <- min(aa[aa > 0])
        counter <- counter + 1
      }
    }
  }
  
  aa <- colSums(ab)
  p <- aa/sum(aa)
  
  m[counter,] <- mom(p,q)
  a[counter] <- x[n] - x[1]
  H[counter] <- shannon(p)
  pmin[counter] <- min(p[p > 0])
  nmin[counter] <- min(aa[aa > 0])
  
  qD <- m ^ (1/(1-matrix(q, nrow = dim(m)[1], ncol = dim(m)[2], byrow = T)))
  qD[, q == 1] <- exp(H)
  
  
  if(regular)
  {
    index <- (1:counter)
    a <- a + 1
  }  else index <- (1:counter)[a[1:counter] > 0]
  
  return(list(mom = m[index,], H = H[index], qD = qD[index,], 
              a = abs(a[index]), q = q, pmin = pmin[index], nmin = nmin[index]))
}