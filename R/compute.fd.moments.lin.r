#' A function to compute moments with a linear transect
#' @param ab vector of abundances
#' @param dist matrix of functional dissimilarities
#' @param q vector of orders
#' @param x vector of sample positions

compute.fd.moments.lin <- function(ab, dist, q, x = 0)
{
  ab <- as.matrix(ab)
  n <- dim(ab)[1]  # nr of samples
  nn <- sum(1:n)  # nr of cells
  regular <- F
  if (length(x) != n) {
    x <- 1:n
    regular <- T
  }
  
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
        
        m[counter,] <- fd.mom.line(p)
        a[counter] <- x[jj+ii-1] - x[jj]
        H[counter] <- fd.shannon(p)
        counter <- counter + 1
      }
    }
  }
  
  aa <- colSums(ab)
  p <- aa/sum(aa)
  
  m[counter,] <- fd.mom.line(p)
  a[counter] <- x[n] - x[1]
  H[counter] <- fd.shannon(p)

  qD <- m ^ (1/(1-matrix(q, nrow = dim(m)[1], ncol = dim(m)[2], byrow = T)))
  qD[, q == 1] <- sqrt(exp(H))
  
  
  if(regular)
  {
    index <- (1:counter)
    a <- a + 1
  }  else index <- (1:counter)[a[1:counter] > 0]
  
  return(list(mom = m[index,], H = H[index], qD = qD[index,], 
              a = abs(a[index]), q = q))
}




compute.FqD.lin <- function(ab, dist, q, x = 0)
{
  ab <- as.matrix(ab)
  n <- dim(ab)[1]  # nr of samples
  nn <- sum(1:n)  # nr of cells
  regular <- F
  if (length(x) != n) {
    x <- 1:n
    regular <- T
  }
  
  dist <- dist[colnames(ab), colnames(ab)]
  
  FqD <- function(x)
  {
    x <- x/sum(x)
    dd <- dist[x > 0, x > 0]
    x <- x[x > 0]
    rao <- x %*% dd %*% x
    temp <- x %*% t(x)
    
    qD <- rep(0, length(q))
    
    if (length(x) == 1) return(qDQ)
    
    for (ii in 1:length(q))
    {
      qD[ii] <- sqrt(( (x^q[ii] %*% dd %*% x^q[ii]) / rao^q[ii] ) ^ (1/(1-q[ii])) / rao)
    }
    
    qD[q == 1] <- sqrt(exp(- sum(dd * temp/rao[1] * log(temp/rao[1]))) / rao)
    
    return(qD)
  }
  
  qD <- matrix(0, nrow = nn, ncol = length(q))
  a <- rep(0, nn)

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
        
        qD[counter,] <- FqD(p)
        a[counter] <- x[jj+ii-1] - x[jj]
        counter <- counter + 1
      }
    }
  }
  
  aa <- colSums(ab)
  p <- aa/sum(aa)
  
  qD[counter,] <- FqD(p)
  a[counter] <- x[n] - x[1]

  if(regular)
  {
    index <- (1:counter)
    a <- a + 1
  }  else index <- (1:counter)[a[1:counter] > 0]
  
  return(list(qD = qD[index,], a = abs(a[index]), q = q))
}