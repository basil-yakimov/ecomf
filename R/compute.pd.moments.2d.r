#' A function to compute phylogenetic moments from a set of points in 2 dimensions
#' @param ab matrix of abundances
#' @param coords 2-column matrix of point coordinates
#' @param tree phylogenetic tree
#' @param itNum number of starting points inside the hull
#' @param print wheather to print point locations map
#' @param r vector of radii
#' @param q vector of orders

compute.pd.moments.2d <- function(ab, coords, tree, itNum = 1000, print = T, r = NULL, q = NULL)
{
  require(sp)
  require(rgeos)
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
  
  sp <- SpatialPoints(coords) 
  buf <- gBuffer(sp, width = 5)
  hull <- gConvexHull(buf)
  env <- gEnvelope(buf)
  
  cent <- gCentroid(sp)
  hullL <- as(hull, "SpatialLines")
  d <- gDistance(cent, hullL)
  
  if (print)
  {
    plot(gEnvelope(buf))
    plot(sp, col = 'black', pch = 19, add = T)
    plot(buf, col = '#00FF0009', lty = 2, add = TRUE)
    plot(hull, col = '#FFFF0009', add = TRUE)
    plot(cent, add = T, col = "red", pch = 19, cex = 2)
    plot(gBuffer(cent, width = d), add = T, col = "#03FF0909")
  }
  
  xr <- env@bbox[1,]
  yr <- env@bbox[2,]
  
  #---#
  
  if (is.null(r))
  {
    dd <- dist(coords, diag = T, upper = T)
    dd <- as.matrix(dd)
    diag(dd) <- rep(1000000, nrow(coords))
    
    dmax <- apply(dd, MARGIN = 1, max)
    dmin <- apply(dd, MARGIN = 1, min)
    
    r <- exp(seq(log(1.5*dmin), log(0.75*dmax), length = 20))
  }
  
  if (is.null(q)) q <- seq(-3, 3, by = 0.1)
  
  mom <- NULL
  A <- sam <- pmin <- nmin <- H <- 0
  counter <- 1
  
  #--------------------------------------#
  
  for (uu in 1:itNum)
  {
    if (uu %% 100 == 0) print(uu)
    nptx <- runif(1, xr[1], xr[2])
    npty <- runif(1, yr[1], yr[2])
    
    newp <- SpatialPoints(t(as.matrix(c(nptx, npty))))
    
    if (!is.na(over(newp, hull)))
    {
      for (ii in r)
      {
        newb <- gBuffer(newp, width = ii)
        isOver <- over(sp, newb)
        sites <- which(!is.na(isOver))
        sam[counter] <- length(sites)
        if (sam[counter] > 2)
        {
          h <- gConvexHull(sp[sites])
          A[counter] <- gArea(h)
          set <- colSums(ab[sites,])
          mom <- rbind(mom, pd.mom.line(set))
          nmin[counter] <- min(set[set > 0])
          pmin[counter] <- nmin[counter]/sum(set)
          H[counter] <- pd.shannon(set)
          counter <- counter + 1
        }
      }
    }
  }
  
  qD <- mom ^ (1/(1-matrix(q, nrow = dim(mom)[1], ncol = dim(mom)[2], byrow = T)))
  qD[, q == 1] <- exp(H)
  
  return(list(mom = mom, H = H, qD = qD, a = A, q = q, pmin = pmin, nmin = nmin, sam = sam))
}