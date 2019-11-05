#' A function to compute moments with a linear transect
#' @param ab matrix of abundances
#' @param dist matrix of functional dissimilarities
#' @param clust matrix of hierarchical clustering
#' @param area vector of areas for clustering levels
#' @param q vector of orders

compute.fd.moments.tri <- function(ab, dist, clust, area, q)
{
  # todo: check clust validity
  
  ab <- as.matrix(ab)
  n <- dim(ab)[1]  # nr of samples
  nlev <- dim(clust)[2] # nr of hierarchy levels
  nn <- sum(apply(clust, 2, max)) # nr of cells
  
  dist <- dist[colnames(ab), colnames(ab)]
  
  fd.mom.line <- function(x)
  {
    x <- x/sum(x)
    dd <- dist[x > 0, x > 0]
    x <- x[x > 0]
    rao <- x %*% dd %*% x
    
    FD.mom <- rep(0, length(q))
    
    if (length(x) == 1) return(FD.mom)
    
    for (ii in 1:length(q))
    {
      FD.mom[ii] <- (x ^ q[ii] %*% dd %*% x^q[ii]) / rao^q[ii]
    }
    
    return(sqrt(FD.mom))
  }
  
  fd.shannon <- function(x)
  {
    x <- x/sum(x)
    dd <- dist[x > 0, x > 0]
    x <- x[x > 0]
    if (length(x) == 1) return(0)
    
    rao <- x %*% dd %*% x
    temp <- x %*% t(x)
    
    return(- sum(dd * temp/rao[1] * log(temp/rao[1])))
  }
  
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
        m[counter, ] <- fd.mom.line(p)
        a[counter] <- area[lev]
        H[counter] <- fd.shannon(p)
        
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
