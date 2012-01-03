
tableMetas <- function(x, agg.clmn, meta.clmns, cond.clmn=NULL, count.clmns=NULL){

  x[,c(meta.clmns, cond.clmn)] <- lapply(x[,c(meta.clmns, cond.clmn)], factor)
  # warning message: "provided #### variables to replace 1 variables"
  #(when cond.clmn is NULL and only one meta column aka 'genus_sp')

  meta <- lapply(meta.clmns, function(mt)
                 by(x, x[,agg.clmn], function(g)
                    sstable(g, c(mt, cond.clmn) ,count.clmns)
                    )
                 )
  names(meta) <- meta.clmns
  return(meta)
}



tableMetaSums <- function(x, meta.clmns, cond.clmn=NULL, count.clmns=NULL){

  x[,c(meta.clmns, cond.clmn)] <- lapply(x[,c(meta.clmns, cond.clmn)], factor)
  
  meta.sum <- lapply(meta.clmns, function(y)
                     sstable(x,c(y,cond.clmn),count.clmns)
                     )
  names(meta.sum) <- meta.clmns
  return(meta.sum)
}

.meta2metasum <- function(meta, pair){

  for(i in names(meta)) meta[[i]]$name <- rownames(meta[[i]])
  tmp <- do.call(rbind.data.frame, meta)
  out <- aggregate(tmp[,c(as.character(pair), 'sum')], list(tmp$name), sum)
  rownames(out) <- out$Group.1
  out$Group.1 <- NULL

  return( out[rev(order(out$sum)),] )
}

counts2manta <- function(x, annotation, a.merge.clmn, agg.clmn, gene.clmns=NULL, meta.clmns=NULL, ...){
		
  if(ncol(x)>2)
    stop('number of columns of x cannot currently exceed 2')
  pair <- colnames(x)  
				
  x <- merge(annotation, x, by.x=a.merge.clmn, by.y='row.names')
  rownames(x) <- x[,a.merge.clmn]
  x <- x[!duplicated(x),]

  is.long <- sapply(annotation[,agg.clmn], nchar) > 255
  if(any(is.long)){
    warning('annotations were found to be over the 256 char limit... truncating')
    annotation[is.long, agg] <- substr(annotation[is.long, agg], 1, 250)
  }
  		
  counts.rg <- regroup(x[,pair],
                         old=annotation[,a.merge.clmn],
                         new=annotation[,agg.clmn],
                         clmns=pair, funcs=rep('sum',length(pair)), combine=FALSE)

  samples <- makeSampleDF(counts.rg)
  genes <- x[match(rownames(counts.rg), x[,agg.clmn]), names(x)[names(x) %in% gene.clmns] ]

  if(is.null(meta.clmns)){
    manta(counts= counts.rg, samples = samples, genes = genes)
  }else{
    manta(counts = counts.rg, samples = samples, genes = genes[, gene.clmns],
          meta = tableMetas(x, agg.clmn=agg.clmn, meta.clmns=meta.clmns, count.clmns=pair),
          meta.sum = tableMetaSums(x,        meta.clmns=meta.clmns, count.clmns=pair), ...
          )
  }
}

                   
align2manta <- function(x, cond.clmn, agg.clmn, gene.clmns, meta.clmns, weight.clmn=NULL, tag.clmn=NULL, ...){

  all.clmns <- c(tag.clmn, cond.clmn, agg.clmn, gene.clmns, meta.clmns)
  
  ## input validation
  if(ncol(x) < 2)
    stop('not enough columns for counts')

  if(!all(sapply(x[,all.clmns], class)=='character'))
    stop("all 'clmns' must be of class 'character'")
  
  if(!all(all.clmns %in% colnames(x)))
    stop("some of your 'clmn' names are not in the column names of your input alignment table (x)")
  
  if(!agg.clmn %in% gene.clmns)
    stop('agg.clmn must be in the gene.clmns')

  count.clmn <- NULL 
  ## check for duplicate tags (if clmn specified) in order to dis-count 
  if(!is.null(tag.clmn)){
    ntags <- length(unique(x[,tag.clmn]))
    if(ntags < nrow(x)){
      warning('number of unique tags is less than table length. discounting...')
      discounts <- 1/table(x[,tag.clmn])
      x$count <- discounts[match(x[,tag.clmn],names(discounts))]
      if(ntags != as.integer(sum(x$count)))  #this check is pretty much totally unnecessary
        stop('something went wrong with the discounting: unique read count != sum of partial counts')
      count.clmn <- 'count'
    }
    if(ntags == nrow(x))
      warning('your tag clmn is unnecessary: number of unique tags == number of rows of input alignment')
  }

  ## cross-tabulate the actual counts
  if(is.null(count.clmn)){   #if no duplicate tags
    counts <- table(x[,agg.clmn], x[,cond.clmn])
  }else{    # if duplicate tags (use the dis-count column)
    counts <- as.table(by(x$count, list(x[,agg.clmn], x[,cond.clmn]), sum))
    counts[is.na(counts)] <- 0
    counts <- round(counts) ## unfortunate but necessary
  }

  ## subset the input table down (same length as counts) into a gene lookup table 
  genes <- x[match(rownames(counts), x[,agg.clmn]), gene.clmns]

  ## generate a meta-level list of gene lists of meta vs treatment cross-tabulation count sub tables (used in pie charts of raPlot)
  meta <- tableMetas(x, agg.clmn, meta.clmns, cond.clmn, count.clmn)
  
  ## now generate the gene-agnostic taxonomic count summary tables
  meta.sum <- tableMetaSums(x, meta.clmns, cond.clmn, count.clmn)
  
  ## add weights
  if(!is.null(weight.clmn)){
    weights <- generateWeights(x, weight.clmn, agg.clmn, cond.clmn, tag.clmn)
  }else{
    weights <- NULL 
  }
  
  manta(counts=counts, samples=makeSampleDF(counts), weights=weights, genes=genes, meta=meta, meta.sum=meta.sum, ...)
}


