#' A function to compute moments with a rectangular lattice of samples
#' @param ab 3d-array of abundances
#' @param q vector of orders

computeMoms2 <- function(ab, q)
{
    L1 <- dim(ab)[1]
    L2 <- dim(ab)[2]
    sizer <- min(L1, L2)
    sp <- dim(ab)[3]
    p <- rep(0, sp)
    siz <- 1:sizer
    
    
    
    moment <- rep(0,length(q)*sizer)
    dim(moment) <- c(sizer, length(q))
    for (aa in sizer:1)
    {
        cur <- L - siz[aa] + 1
        sqNum <- 0
        for (ii in 0:(cur-1))
        {
            for (jj in 0:(cur-1))
            {
                p <- rep(0, sp)
                sqAb <- 0
                for (iii in (ii+1) : (ii+siz[aa]) )
                {
                    for (jjj in (jj+1) : (jj+siz[aa]) )
                    {
                        p <- p + n[iii,jjj,]
                    }
                }
                p <- p[p > 0]/sum(p)
                for (uu in 1:length(q))
                {
                    if (log) moment[aa,uu] <- moment[aa,uu] + log(sum(p^q[uu]))
                    else moment[aa,uu] <- moment[aa,uu] + sum(p^q[uu])
                }
                sqNum <- sqNum + 1
            }
        }
        moment[aa,] <- moment[aa,]/sqNum
    }
    if (log) moment <- exp(moment)
    return(list(mom = moment, A = siz^2, q = q))
}

