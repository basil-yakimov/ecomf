#' A function to compute moments with a linear transect
#' @param ab matrix of abundances
#' @param clust matrix of hierarchical clustering
#' @param area vector of areas for clustering levels
#' @param q vector of orders

compute.moments.tri <- function(ab, clust, area, q)
{
  # todo: check clust validity
  
  ab <- as.matrix(ab)
  n <- dim(ab)[1]  # nr of samples
  nlev <- dim(clust)[2] # nr of hierarchy levels
  nn <- sum(apply(clust, 2, max)) # nr of cells
  
  m <- qD <- matrix(0, nrow = nn, ncol = length(q))
  a <- rep(0, nn)
  H <- rep(0, nn)
  pmin <- nmin <- rep(0, nn)
  
  counter <- 1 # current cell
  
  for (lev in 1:nlev)
  {
    for (ii in 1:max(clust[, lev], na.rm = T))
    {
      p <- ab[clust[, lev] == ii, ]
      if (is.matrix(p)) p <- colSums(p)
      if (sum(p) > 0)
      {
        m[counter, ] <- mom(p, q)
        a[counter] <- area[lev]
        H[counter] <- shannon(p)
        
        nmin[counter] <- min(p[p > 0])
        pmin[counter] <- nmin[counter]/sum(p)
        counter <- counter + 1
      }
    }
  }
  
  qD <- m ^ (1/(1-matrix(q, nrow = dim(m)[1], ncol = dim(m)[2], byrow = T)))
  qD[, q == 1] <- exp(H)
  
  counter <- counter - 1
  
  return(list(mom = m[1:counter,], H = H[1:counter], qD = qD[1:counter,], 
              a = a[1:counter], q = q, pmin = pmin[1:counter], nmin = nmin[1:counter]))
}
