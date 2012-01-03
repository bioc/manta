makeSampleDF <- function(counts, group=factor(colnames(counts)), lib.size=colSums(counts)){

  if (ncol(counts) != length(group)) 
    stop("Length of 'group' must equal number of columns in 'counts'")

  if (ncol(counts) != length(lib.size)) 
    stop("Length of 'lib.size' must equal number of columns in 'counts'")

  group <- as.factor(group)
  
  samples <- data.frame(group = group, lib.size = lib.size)
  row.names(samples) <- colnames(counts)
  
  return(samples)
}

.geomean.cmn.lib.size <- function(obj){
   prod(obj$samples$lib.size)^(1/ncol(obj$counts))
}





.setLibrarySizes <- function(obj, lib.sizes){

  ## rebuild the sample dataframe using the new specified lib sizes
  obj$samples <- data.frame(group = obj$samples$group, lib.size = lib.sizes)

  ## now reset the common library size via a geometric mean (from edgeR's 'equalizeLibSizes')
  obj$common.lib.size <- .geomean.cmn.lib.size(obj)

  return(obj)

}
