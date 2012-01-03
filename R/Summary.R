





summary.manta <- function(object, ...){

  if(!'meta' %in% names(object))
    stop('could not find a meta slot in your object')

  lapply(object$meta.sum, function(ms)  leghead(ms, nrow(ms)))


  #print('calculating abundance shifts')
  #sapply(names(obj$manta), function(ml) mantaMethod(obj, ml, pair))
}


.getCondFactor <- function(v, pair, cond.lev='ref')
  pair[as.integer(1:length(v) %in% grep(pair[names(pair)==cond.lev], names(v)))+1]


normfact2absTMM <- function(x, pair, f=nv(x$samples, 'norm.factors'), sums=colSums(x$counts))
  log2(f[pair['obs']] * do.call('/', as.list(by(sums, .getCondFactor(sums, pair, 'ref'), mean)[1:2]))) 




meta2counts <- function(obj, meta.lev, rm.sum=TRUE, meta.subset=NULL){

  counts <- do.call(rbind, obj$meta[[meta.lev]])
  genes <- rep(names(obj$meta[[meta.lev]]), sapply(obj$meta[[meta.lev]], nrow))
  meta <- unlist(sapply(obj$meta[[meta.lev]], rownames))

  df <- cbind.data.frame(genes, meta, counts)
  rownames(df) <- NULL

  if(rm.sum)
    df$sum <- NULL

  if(!is.null(meta.subset)){
    df <- subset(df, meta == meta.subset)
    df$meta <- NULL
    rownames(df) <- df$genes
    df$genes <- NULL
  }  
  
  return(df)  

}
