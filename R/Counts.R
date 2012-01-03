
collapseRepliCounts <- function(x, pair=nv(levels(x$samples$group)[1:2] , c('ref','obs'))){

  if(! class(x) %in% c('manta','DGEList'))
    stop('x must be of class manta or DGEList')
  
  if('pseudo.alt' %in% names(x)){
    cts <- x$pseudo.alt
  }else{
    warning("You're trying to collapse replicates, yet no psudo.alt exists in your object. \
               Consider running estimateCommonDisp first!.  Using raw counts instead.")
    cts <- x$counts
  }
  
  .checkCondPair(x, pair)
  
  if(dim(cts)[2] != 2){
    cg.f <- x$samples$group  
    o <- apply(cts[,cg.f==pair['obs']], 1, sum, na.rm=T)
    r <- apply(cts[,cg.f==pair['ref']], 1, sum, na.rm=T)
  }else{
    o <- cts[,as.character(pair['obs'])]
    r <- cts[,as.character(pair['ref'])]
  }
  
  xor <- data.frame(obs=o,ref=r)
  names(xor) <- as.character(pair[c('obs','ref')])
  
  return(xor)
}


