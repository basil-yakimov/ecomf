#' A function to compute moments with a rectangular lattice of samples
#' @param ab 3d-array of abundances
#' @param q vector of orders

compute.moments.lat <- function(ab, q)
{
    L1 <- dim(ab)[1]
    L2 <- dim(ab)[2]
    sizer <- min(L1, L2)
    sp <- dim(ab)[3]
    p <- rep(0, sp)
    siz <- 1:sizer
    
    nn <- sum((L1 - siz + 1)*(L2 - siz + 1))  # nr of cells
    
    m <- qD <- matrix(0, nrow = nn, ncol = length(q))
    a <- rep(0, nn)
    H <- rep(0, nn)
    pmin <- nmin <- rep(0, nn)
    
    counter <- 1
    for (uu in siz)
    {
      for (ii in 1:(L1 - uu + 1)) {
        for (jj in 1:(L2 - uu + 1)) {
          aa <- apply(ab[ii:(ii+uu-1), jj:(jj+uu-1), , drop = F], 3, sum)
          if (sum(aa) > 0)
          {
            p <- aa/sum(aa)
            
            m[counter,] <- mom(p,q)
            a[counter] <- uu
            H[counter] <- shannon(p)
            pmin[counter] <- min(p[p > 0])
            nmin[counter] <- min(aa[aa > 0])
            counter <- counter + 1
          }
        }
      }
    }
    
    qD <- m ^ (1/(1-matrix(q, nrow = dim(m)[1], ncol = dim(m)[2], byrow = T)))
    qD[, q == 1] <- exp(H)
    
    return(list(mom = m, H = H, qD = qD, 
                a = a^2, q = q, pmin = pmin, nmin = nmin))
}
