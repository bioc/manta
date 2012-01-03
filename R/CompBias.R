
.compbias <- function(x, pair=nv(levels(x$samples$group)[1:2], c('ref','obs')), meta.lev='phylum',  meta.lev.lim=min(5,nrow(x$meta.sum[[meta.lev]]))){

  .checkMetaLev(x, meta.lev)
  .checkCondPair(x, pair)
  
  sub.tax.nms <- rownames(x$meta.sum[[meta.lev]])
  RAy <- raPlot(x$counts[,pair[c('ref','obs')]], jitter=FALSE, plot=FALSE)
  
  taxa.cts <- sapply(x$meta[[meta.lev]], nv,'sum')
  
  R.dist <- list()

  taxa <- sub.tax.nms[1:meta.lev.lim]
  for(tax in taxa){
    tax.cts <- sapply(taxa.cts, '[', tax)
    R.dist[[tax]] <- RAy$R[!is.na(tax.cts)] #rep( , tax.cts[!is.na(tax.cts)])
  }

  return(R.dist)

}


compbiasPlot <- function(x, pair=nv(levels(x$samples$group)[1:2] , c('ref','obs')), meta.lev='phylum',  meta.lev.lim=nrow(x$meta.sum[[meta.lev]]), breaks=10, xlim=c(0,.5), type='box', col=NULL,...){

  if(is.null(col)){
    if('color' %in% names(x$meta.sum[[meta.lev]]))
      col <- nv(x$meta.sum[[meta.lev]],'color')
    else
      col <- rainbow(meta.lev.lim)
  }
  R.dist <- .compbias(x, pair=pair, meta.lev=meta.lev, meta.lev.lim=meta.lev.lim)
  
  
  taxa <- names(R.dist)
  if(type == 'dist'){
    ylims <- range(RAy$R)
    plot(0,0, pch='', ylim=ylims, xlim=xlim, ylab='R',xlab='density')
    sapply(taxa, function(tx) with(hist(R.dist[[tx]], plot=FALSE, breaks=breaks),
                                   lines(density, mids, col=col[tx], ylim=ylims, xlim=xlim,  ...)))
    if(legend)
      legend(x='bottomright', legend=taxa, fill=col)
  }else if(type=='box'){
    boxplot(R.dist,  col=col[taxa], ...)
  }else if(type=='violin'){
    violins(R.dist,  col=col[taxa], ...)
  }


}



compbiasTest <- function(x, pair=nv(levels(x$samples$group)[1:2] , c('ref','obs')), meta.lev='phylum',  meta.lev.lim=min(10,nrow(x$meta.sum[[meta.lev]]))){

  R.dist <- .compbias(x, pair=pair, meta.lev=meta.lev, meta.lev.lim=meta.lev.lim)

  R <- unlist(R.dist)
  subtaxa <- rep(names(R.dist), sapply(R.dist,length))

  anova(lm(R~subtaxa))
}

