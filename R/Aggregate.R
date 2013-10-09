

.aggBinCounts <- function(x, pair=nv(levels(x$samples$group)[1:2] , c('ref','obs')), breaks=c(-1,0,1,9,99, max(x[,pair])), clmn.labs=c('0','(0,1]','(1-10]','(10,100]','(100,max]')){
  out <- rbind.data.frame(
                          hist(x[,pair[1]], breaks=breaks, plot=F)$counts,
                          hist(x[,pair[2]], breaks=breaks, plot=F)$counts
                          )
  names(out) <- clmn.labs
  rownames(out) <- as.character(pair)
  return(out)
}


.aggDESigCumDist <- function(x, pair=nv(levels(x$samples$group)[1:2] , c('ref','obs')), breaks=c(0,.001,.01, .05, 1), labels, nr=0, fc.clmn='R', sig.clmn='fdr'){
  out <- rbind.data.frame(
                          cumsum(hist(x[,sig.clmn], breaks=breaks, plot=F, right=F)$counts),
                          cumsum(hist(subset(x, x[,fc.clmn] <  nr)[,sig.clmn], breaks=breaks, plot=F, right=F)$counts),
                          cumsum(hist(subset(x, x[,fc.clmn] >= nr)[,sig.clmn], breaks=breaks, plot=F, right=F)$counts) 
                          )
  names(out) <- c('<=.001','<=.01','<=.05','<=1')
  rownames(out) <- c('total',pair)
  return(out)
}