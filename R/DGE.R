
.as.manta <- function(x){

  if(F){
  ## UNFINISHED
  align2manta()
  #.DGEList.as.manta()

    x <- new("DGEList", list(samples = samples, counts = counts))
 }

  
}




.as.DGEList <- function(x){

  #names(x$samples) <- sub('condition','group',names(x$samples))
  if(F){
    if('common.dispersion' %in% names(x))
      cd <- x$common.dispersion
    else
      cd <- 1e-16
    
    new('DGEList', list(counts=x$counts, samples=x$samples, common.dispersion=cd))
  }

  stop('as.DGEList has been depricated')
}




