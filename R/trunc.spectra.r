#' A function to truncate anomalous "horns" from multifractal spectra
#' @param mf spectra object returned by local.spectra() or compute.spectrum()

trunc.spectra <- function(mf)
{
  if (is.matrix(mf$f))
  {
    n <- ncol(mf$f)
    len <- nrow(mf$f)
    start <- which(mf$q == 0)
    for (ii in 1:n)
    {
      pos <- start - 1
      while (mf$f[pos, ii] <= mf$f[pos+1, ii] & pos != 1) pos <- pos - 1
      if (pos != 1) 
      {
        mf$f[1:pos, ii] <- NA
        mf$alfa[1:pos, ii] <- NA
      }
      pos <- start + 1
      while (mf$f[pos, ii] <= mf$f[pos-1, ii] & pos != len) pos <- pos + 1
      if (pos != len) 
      {
        mf$f[pos:len, ii] <- NA
        mf$alfa[pos:len, ii] <- NA
      }
    }
  } else
  {
    len <- length(mf$f)
    start <- which(mf$q == 0)
    pos <- start - 1
    while (mf$f[pos] <= mf$f[pos+1] & pos != 1) pos <- pos - 1
    if (pos != 1) 
    {
      mf$f[1:pos] <- NA
      mf$alfa[1:pos] <- NA
    }
    pos <- start + 1
    while (mf$f[pos] <= mf$f[pos-1] & pos != len) pos <- pos + 1
    if (pos != len) 
    {
      mf$f[pos:len] <- NA
      mf$alfa[pos:len] <- NA
    }
  }
  return(mf)
}