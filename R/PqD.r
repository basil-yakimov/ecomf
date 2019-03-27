#' A function to compute phylogenetic diversity indices
#' @param dat matrix of abundances
#' @param tree phylogenetic tree
#' @param q vector of orders
#' @param hill transformation to phylogenetic hill numbers

PqD <- function(dat, tree, q = 0:2, hill = T)
{
  require(geiger)
  
  PqD.line <- function(x)
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
        qD[ii] <- sum(branches[, 3]*(branches[, 4]/TT)^q[ii]) ^ (1/(1-q[ii]))
      }
      qD[q == 1] <- exp(-sum(branches[,3]*branches[,4]/TT*log(branches[,4]/TT)))
      
      if (hill) qD <- qD/TT
      
      return(qD)
    }
    else return(rep(0, length(q)))
  }
  
  dat <- as.matrix(dat)
  nr <- nrow(dat)
  nc <- ncol(dat)
  if (nc == 1) {
    dat <- t(dat)
    nr <- nrow(dat)
    nc <- ncol(dat)
  }
  
  out <- t(apply(dat, 1, PqD.line))
  
  if (length(q) == 1) out <- as.vector(out)
  
  return(out)
}