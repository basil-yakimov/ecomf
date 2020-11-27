#' A function to compute moments with a linear transect
#' @param ab vector of abundances
#' @param tree phylogenetic tree
#' @param q vector of orders
#' @param x vector of sample positions

compute.pd.moments.lin <- function(ab, tree, q, x = 0)
{
  ab <- as.matrix(ab)
  n <- dim(ab)[1]  # nr of samples
  nn <- sum(1:n)  # nr of cells
  regular <- F
  if (length(x) != n) {
    x <- 1:n
    regular <- T
  }
  
  require(geiger)
  
  pd.mom.line <- function(x)
  {
    x <- x/sum(x)
    if (sum(x > 0) > 1)
    {
      tmp.tree <- treedata(tree, x[x > 0], warnings = F)$phy
      branches <- matrix(NA, nrow = nrow(tmp.tree$edge), ncol = 4)
      branches[, 1:2] <- tmp.tree$edge
      branches[, 3] <- tmp.tree$edge.length
      for (ii in 1:nrow(branches))
      {
        leaves.node <- tips(tmp.tree, branches[ii, 2])
        branches[ii, 4] <- sum(x[leaves.node], na.rm = T)
      }
      TT <- max(branching.times(tmp.tree))
      PD.mom <- rep(0, length(q))
      for (ii in 1:length(q))
      {
        PD.mom[ii] <- sum(branches[, 3]*(branches[, 4]/TT)^q[ii])
      }
      
      return(PD.mom)
    }
    else
    {
      return(rep(0, length(q)))
    }
  }
  
  pd.shannon <- function(x)
  {
    x <- x/sum(x)
    if (sum(x > 0) > 1)
    {
      tmp.tree <- treedata(tree, x[x > 0], warnings = F)$phy
      branches <- matrix(NA, nrow = nrow(tmp.tree$edge), ncol = 4)
      branches[, 1:2] <- tmp.tree$edge
      branches[, 3] <- tmp.tree$edge.length
      for (ii in 1:nrow(branches))
      {
        leaves.node <- tips(tmp.tree, branches[ii, 2])
        branches[ii, 4] <- sum(x[leaves.node], na.rm = T)
      }
      TT <- max(branching.times(tmp.tree))
      H <- -sum(branches[,3]*branches[,4]/TT*log(branches[,4]/TT))
      return(H)
    }
    else
    {
      return(0)
    }
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
        
        m[counter,] <- pd.mom.line(p)
        a[counter] <- x[jj+ii-1] - x[jj]
        H[counter] <- pd.shannon(p)
        counter <- counter + 1
      }
    }
  }
  
  aa <- colSums(ab)
  p <- aa/sum(aa)
  
  m[counter,] <- pd.mom.line(p)
  a[counter] <- x[n] - x[1]
  H[counter] <- pd.shannon(p)

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




compute.PqD.lin <- function(ab, dist, q, x = 0)
{
  ab <- as.matrix(ab)
  n <- dim(ab)[1]  # nr of samples
  nn <- sum(1:n)  # nr of cells
  regular <- F
  if (length(x) != n) {
    x <- 1:n
    regular <- T
  }
  
  require(geiger)
  
  PqD <- function(x)
  {
    x <- x/sum(x)
    if (sum(x > 0) > 1)
    {
      tmp.tree <- treedata(tree, x[x > 0], warnings = F)$phy
      branches <- matrix(NA, nrow = nrow(tmp.tree$edge), ncol = 4)
      branches[, 1:2] <- tmp.tree$edge
      branches[, 3] <- tmp.tree$edge.length
      for (ii in 1:nrow(branches))
      {
        leaves.node <- tips(tmp.tree, branches[ii, 2])
        branches[ii, 4] <- sum(x[leaves.node], na.rm = T)
      }
      TT <- max(branching.times(tmp.tree))
      qD <- rep(0, length(q))
      for (ii in 1:length(q))
      {
        qD[ii] <- sum(branches[, 3]*(branches[, 4]/TT)^q[ii]) ^ (1/(1-q[ii])) / TT
      }
      qD[q == 1] <- exp(-sum(branches[,3]*branches[,4]/TT*log(branches[,4]/TT))) / TT
      
      return(qD)
    }
    else return(rep(0, length(q)))
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
        
        qD[counter,] <- PqD(p)
        a[counter] <- x[jj+ii-1] - x[jj]
        counter <- counter + 1
      }
    }
  }
  
  aa <- colSums(ab)
  p <- aa/sum(aa)
  
  qD[counter,] <- PqD(p)
  a[counter] <- x[n] - x[1]

  if(regular)
  {
    index <- (1:counter)
    a <- a + 1
  }  else index <- (1:counter)[a[1:counter] > 0]
  
  return(list(qD = qD[index,], a = abs(a[index]), q = q))
}