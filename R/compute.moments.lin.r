#' A function to compute moments with a linear transect
#' @param ab vector of abundances
#' @param cv vector of coverages
#' @param ht vector of heights
#' @param q vector of orders

compute.moments.lin <- function(ab, cv, ht, q)
{
  ab <- as.matrix(ab)
  cv <- as.matrix(cv)
  ht <- as.matrix(ht)
  n <- dim(ab)[1]  # nr of samples
  if (dim(cv)[1] != n | dim(ht)[1] != n) error("data dimentions must be equal")
  nn <- sum(1:n)  # nr of cells
  
  m <- matrix(0, nrow = nn, ncol = length(q))
  a <- rep(0, nn)
  H <- rep(0, nn)
  pmin <- rep(0, nn)
  
  counter <- 1
  for (ii in 1:(n-1))
  {
    for (jj in 1:(n-ii+1))
    {
      aa <- ab[jj:(jj+ii-1), ]
      cc <- cv[jj:(jj+ii-1), ]
      hh <- ht[jj:(jj+ii-1), ]
      
      if (is.matrix(aa))
      {
        aa <- colSums(aa)
        cc <- colSums(cc)
        hh <- colSums(hh)
      }
      if (sum(aa) > 0)
      {
        p <- ( aa/sum(aa) + cc/sum(cc) + hh/sum(hh) ) / 3
        
        
        m[counter,] <- mom(p,q)
        a[counter] <- ii
        H[counter] <- shannon(p)
        pmin[counter] <- min(p[p>0])
        
        counter <- counter + 1
      }
    }
  }
  
  aa <- colSums(ab)
  cc <- colSums(cv)
  hh <- colSums(ht)
  p <- ( aa/sum(aa) + cc/sum(cc) + hh/sum(hh) ) / 3
  
  m[counter,] <- mom(p,q)
  a[counter] <- n
  H[counter] <- shannon(p)
  pmin[counter] <- min(p[p>0])
  
  return(list(mom = m[1:counter,], H = H[1:counter], a = a[1:counter], q = q, pmin = pmin[1:counter]))
}