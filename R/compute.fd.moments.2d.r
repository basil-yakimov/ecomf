#' A function to compute functional moments from a set of points in 2 dimensions
#' @param ab matrix of abundances
#' @param coords 2-column matrix of point coordinates
#' @param dist matrix of functional dissimilarities
#' @param itNum number of starting points inside the hull
#' @param print wheather to print point locations map
#' @param r vector of radii
#' @param q vector of orders

compute.fd.moments.2d <- function(ab, coords, dist, itNum = 1000, print = T, r = NULL, q = NULL)
{
  require(sp)
  require(rgeos)

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
          mom <- rbind(mom, fd.mom.line(set))
          nmin[counter] <- min(set[set > 0])
          pmin[counter] <- nmin[counter]/sum(set)
          H[counter] <- fd.shannon(set)
          counter <- counter + 1
        }
      }
    }
  }
  
  qD <- mom ^ (1/(1-matrix(q, nrow = dim(mom)[1], ncol = dim(mom)[2], byrow = T)))
  qD[, q == 1] <- sqrt(exp(H))
  
  return(list(mom = mom, H = H, qD = qD, a = A, q = q, pmin = pmin, nmin = nmin, sam = sam))
}


compute.FqD.2d <- function(ab, coords, dist, itNum = 1000, print = T, r = NULL, q = NULL)
{
  require(sp)
  require(rgeos)
  
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
  
  qD <- NULL
  A <- sam <- pmin <- nmin <- 0
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
          qD <- rbind(qD, FqD(set))
          nmin[counter] <- min(set[set > 0])
          pmin[counter] <- nmin[counter]/sum(set)
          counter <- counter + 1
        }
      }
    }
  }
  
  return(list(qD = qD, a = A, q = q, pmin = pmin, nmin = nmin, sam = sam))
}