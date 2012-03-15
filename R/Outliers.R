outGenes <- function(x, n=50, p=.05, FC=1, A.pct=.05, uk.filter=NULL, method='BH',  verbose=TRUE){

  is.DGE.test <- class(x) == 'DGEExact'
  if(is.DGE.test)
    x <- x$table

  #if(is.

  
  orig.nms <- names(x)
  
  if(is.DGE.test)
    names(x) <- sub('logCPM','A', sub('logFC','R', names(x)))


  
  x$adj.p <- p.adjust(x$PValue, method=method)

  if(!is.null(uk.filter))
    x <- x[apply(sapply(uk.filter, function(f) !grepl(f,rownames(x))),1, all),]

  x.in <- x
  x <- subset(x, PValue <= p  &  (R > FC | R < -FC))

  if(n <= 0)
    n <- nrow(x)
  
  n.a <- round(A.pct * n)
  n.r <- n - n.a

  r.up <- subset(x, R > 0)
  r.dn <- subset(x, R < 0)

  n.r.up <- min(nrow(r.up), n.r/2)
  n.r.dn <- min(nrow(r.dn), n.r/2)
  
  a <- x.in[rev(order(x.in$A))[1:n.a],]
    
  r.up <- r.up[rev(order(r.up$R)), ][1:n.r.up,]
  r.dn <- r.dn[    order(r.dn$R) , ][1:n.r.dn,]

  if(verbose){

    clms <- grep('p.value', names(x), invert=TRUE)
    
    cat('obs/ref (up-regulated)', sep='\n')
    print(r.up[,clms]); cat('\n')
    
    cat('ref/obs (down-regulated)', sep='\n')
    print(r.dn[,clms]); cat('\n')
    
    cat('high abundance (house-keeping)', sep='\n')
    print(a[,clms]) ; cat('\n')   
  }


  
  out <- rbind.data.frame(r.up, r.dn, a)

  names(out) <- c(orig.nms,'adj.p')
  out$otlr <- c(rep('up',nrow(r.up)),rep('dn',nrow(r.dn)), rep('rt',nrow(a)))

  invisible(out)

}


