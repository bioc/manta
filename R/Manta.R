

require(methods)
setClass("manta", contains=c("DGEList"), representation("list"))

manta <- function(counts, samples=makeSampleDF(counts), genes=NULL, meta=NULL, meta.sum=NULL, weights=NULL, norm=TRUE, disp=TRUE, ...){

  if(ncol(counts) < 2)
    stop('must have at least two columns in manta$counts')
  
  # begin stolen edgeR's DGEList validation bit
  counts <- as.matrix(counts)
  if (ncol(counts) > 0 & is.null(colnames(counts))) 
    colnames(counts) <- paste("sample", 1:ncol(counts), sep = ".")
  
  obj <- new('manta', list(counts=counts, samples=samples), ...)  
  
  
  if(!is.null(genes)) {
    genes <- as.data.frame(genes, stringsAsFactors = FALSE)
    if(nrow(genes) != nrow(counts)){
      stop("counts and genes have different nrows")
    }
    rownames(genes) <- rownames(counts)
    obj$genes <- genes
  }
  ## end stolen edgeR DGEList validation bit

  #if(is.null(common.lib.size))
    obj$common.lib.size <- .geomean.cmn.lib.size(obj)  # move to top as param default?
  #else
  #  obj$common.lib.size <- common.lib.size

  if(norm)
    obj <- calcNormFactors(obj)
  if(disp)
    obj <- estimateCommonDisp(obj)

  if(!is.null(weights)){
    if(!all(dim(weights) == dim(obj$counts)))
      stop('weights df must be the same dims as counts')
    obj$weights <- weights
  }
  ## check meta now too
  if(!is.null(meta)){
    if (all(sapply(meta, length) != nrow(counts)))      
      warning("counts rows and meta lists have different lengths")
    if (max(sapply(meta, length)) != nrow(counts))
      stop("at least one of the meta lists must equal the counts rows")
    obj$meta <- meta
  }             

  ## check meta summary tables now too
  if(!is.null(meta.sum)){
    if (length(meta.sum) != length(obj$meta))
      stop("meta and meta.sum lists have different lengths")    
    obj$meta.sum <- meta.sum
  }             
  return(obj)
}
 
