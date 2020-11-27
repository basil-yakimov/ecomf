#' A function to compute moments with a linear transect
#' @param ab matrix of abundances
#' @param tree phylogenetic tree
#' @param clust matrix of hierarchical clustering
#' @param area vector of areas for clustering levels
#' @param q vector of orders

compute.pd.moments.tri <- function(ab, tree, clust, area, q)
{
  # todo: check clust validity
  
  ab <- as.matrix(ab)
  n <- dim(ab)[1]  # nr of samples
  nlev <- dim(clust)[2] # nr of hierarchy levels
  nn <- sum(apply(clust, 2, max)) # nr of cells
  
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
        m[counter, ] <- pd.mom.line(p)
        a[counter] <- area[lev]
        H[counter] <- pd.shannon(p)
        
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
