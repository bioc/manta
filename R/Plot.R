plot.manta <- manta.ra <- function(x, uniques=4, pair=nv(levels(x$samples$group)[1:2] , c('ref','obs')),
                       nr=0, alpha = 0.01, normalize=FALSE, 
                       meta.level=names(x$meta.sum)[1], meta.lgnd.lim=6, lgd.pos='topright', lgd.cex=.75, lgd.trunk=FALSE, pie.lwd=1,
                       annot=NULL, vrb.axlabs=TRUE, jitter=.43, border='black',
                       rex=2, flat=FALSE, tail=.5, arms=.5, spine=1, ...){
  
  ## try using missing()
  .checkMetaLev(x, meta.level)
  .checkCondPair(x, pair)
  
  if('weights' %in% names(x)){
    jit.wgts <- x$weights[,pair[c('ref','obs')]]
  }else{
    jit.wgts <- NULL
  }

  
  if(meta.level != FALSE & 'meta.sum' %in% names(x)){
    if(nrow(x$meta.sum[[meta.level]])==1){
      ## disable the pies
      warning('Only one meta sub-level. Defaulting to flat pie-less raPlot')
      flat <- TRUE
      meta.level <- FALSE
    }
  }

  if('nr' %in% names(x) & nr==0){
    message('using normalization ratio [nr] from within manta object')
    nr <- x$nr
  }
  
  ab <- collapseRepliCounts(x, pair)  
  rownames(ab) <- rownames(x$counts)
 
  
  RAy <- raPlot(ab[,as.character(pair[c('ref','obs')])], uniques=uniques, normalize=normalize,  
         nr=nr, alpha = alpha, jitter=jitter, jit.wgts= jit.wgts, 
         rex=rex, flat=flat, tail=tail, arms=arms, spine=spine, border=border,...)
  
  
  if(vrb.axlabs)
    raAddAxLabs(pair, normalize, line=2)
    
  ## now use these to overlay the pie charts
  if('meta' %in% names(x) & meta.level!=FALSE){
    if(!jitter)
      message("you should try setting 'jitter=TRUE' to show the hidden pies")
    
    tables <- lapply(x$meta[[meta.level]], nv, 'sum')
    colors <- nv(leghead(x$meta.sum[[meta.level]], n=meta.lgnd.lim, na.col='white', other.col='gray', na.name='unclassified'),'color')
    
    par(new=TRUE)
    with(RAy[names(tables),], 
       pies(x = tables, x0 = A, y0 = R,
         radii = sizes * rex, #xlim=xlims, ylim=ylims,
          color.table=colors, border=border, lwd=pie.lwd)
    )
    
    if(!is.null(lgd.pos) & as.logical(meta.lgnd.lim)){
      if(lgd.trunk){
        abrv.gspecies <- function(x){
          if(length(grep('sp\\.', x))== 0 & length(grep('spp\\.', x)) == 0){
            a <- strsplit(x,' ')[[1]]
            x <- paste(substr(a[1],1,1),a[2], sep='.')
          }
          x
        }
        names(colors) <- sapply(names(colors), abrv.gspecies)
      }
      legend(lgd.pos, inset=.02, fill=colors, border=gray(.3), legend=names(colors), cex=lgd.cex)
    }
  }

  ## now annotate with text over the top of the new pies
  if(!is.null(annot) & length(annot)!=0){
    if(!is.character(annot))
      stop('annot must be a character vector of names that match the count table rownames')
    an.subset <- match(annot, names(RAy$R))
    with(RAy, text(A[an.subset], R[an.subset], annot, col=names(annot)))
  }
 
  with(RAy, invisible(data.frame(A=A, R=R, sizes=sizes*rex, row.names=rownames(RAy))))
}



