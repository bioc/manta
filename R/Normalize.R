
nr <- function(obj, meta.lev, pair=nv(levels(x$samples$group)[1:2], c('ref','obs'))){

  .checkMetaLev(obj, meta.lev)

  nrms <- numeric()
  nams <- rownames(obj$meta.sum[[meta.lev]])
  ## for each meta sub level, calculate the normalization factors
  for(meta.sub.lev in nams){
    cts <- meta2counts(obj, meta.lev, meta.subset=meta.sub.lev)
    f <- calcNormFactors(as.matrix(cts))

    if(!is.null(pair)){
      f <- c(f, nv( normfact2absTMM(x=cts, pair=pair, f=f, sums=colSums(cts)), 'nr'))
      ## normal odds ratios
    }
    nrms <- rbind(nrms, f)
  }
  rownames(nrms) <- nams

  return(nrms)
}

.getCondGroups <- function(v, pair)
  sapply(v,function(ve) pair[sapply(pair, function(c) length(grep(c,ve))>0)])


.abs.nf <- function(x, f, gf=x$samples$group){
  s <- by(apply(x$counts, 2, sum), gf, mean)
  f <- f * do.call('/',as.list(s))^-1
  return(f)
}

nf2nr <- function(x, pair=nv(levels(x$samples$group)[1:2], c('ref','obs')), method='nf', absolute=TRUE){

  .checkCondPair(x, pair)
  
  cnds <- pair[c('obs','ref')]
  if(class(x) == 'DGEList')
    gf <- nv(x$samples,'group')
  else if (class(x) == 'manta')
    gf <- nv(x$samples, 'group')
  else  # use only if also remove the samples and counts dependancies below
    gf <- .getCondGroups(colnames(x$counts), pair)

  nf <- switch(method,
               nf= sapply(cnds,
                 function(cnd)
                 geomean(nv(x$samples,'norm.factors')[gf==cnd])),
               ctagg = calcNormFactors(
                 sapply(cnds,
                        function(cnd)
                        rowSums(x$counts[,gf==cnd])))
               )
  f <- do.call('/', as.list(nf))
  
  if(absolute)
    f <- .abs.nf(x, f, gf)
  nr <- log2(f)
  if(is.na(nr) | nr==0 | abs(signif(nr,7)) == 3.203427e-16)  #may not be appropriate
    nr <- NA
  return(nr)
}



#nr2nf <- function(nr){
#
#}

.wtd.var <- function (x, weights = NULL, normwt = FALSE, na.rm = TRUE)
{
    if (!length(weights)) {
        if (na.rm)
            x <- x[!is.na(x)]
        return(var(x))
    }
    if (na.rm) {
        s <- !is.na(x + weights)
        x <- x[s]
        weights <- weights[s]
    }
    if (normwt)
        weights <- weights * length(x)/sum(weights)
    xbar <- sum(weights * x)/sum(weights)
    sum(weights * ((x - xbar)^2))/(sum(weights) - 1)
}



.normalize <- function(x, meta.lev, pair=nv(levels(x$samples$group)[1:2], c('ref','obs')), method=c("mtd.balanced",'mtd.balnorm',"mtd.normal","mtd.raw","TMM", "RLE", "quantile","mean","median"), expnt=NA, uniques){

  .checkCondPair(x, pair)

  nr <- NA
  vr <- sd <- NA #variance
  method <- match.arg(method)
  
  if(class(x)!='manta')
    stop('x must be of class manta')

  raw.methods <- c('mean', 'median')
  is.mtd <- length(grep('mtd',method))>0
  if(is.mtd | method %in% raw.methods){

    rawR <- .rawR(x, pair, uniques)
    if(!is.null(rawR)){
      if(is.mtd){        
        ##if(nrow(x$meta.sum[[meta.lev]])==1)
        ##  warning('only one possible level weighting: reverting to an unweighted mean')
        mtd.method <- strsplit(method,'\\.')[[1]][2]
        if(is.na(mtd.method))
          mtd.method <- 'balanced'
        mtd <- .mtd(x, meta.lev=meta.lev, pair=pair, method=mtd.method, expnt=expnt)
        Rmtd <- nerge(list(rawR=rawR, mtd=mtd))
        
        nr <- with(Rmtd, weighted.mean( rawR, mtd))
        vr <- with(Rmtd, wtd.var(rawR, mtd)) # Hmisc
      }else{
        nr <- eval(call(method, rawR))
        if(method=='mean')
          vr <- var(rawR)
      }
    }
  }else if(method %in% c("TMM", "RLE", "quantile")){
    tmp.nf <- try(calcNormFactors(as.matrix(x$counts),method=method)) #, refColumn=pair['ref']
    if(class(tmp.nf)=='try-error'){
      nr <- NA
    }else{
      x$samples$norm.factors <- tmp.nf
      nr <- nf2nr(x, pair=pair)
    }
    
    if(!is.na(nr))
      if(nr < 10e-12 & nr > -10e-12) #TMM and quantile sometimes return funny almost zero values
        nr <- NA
    #if(method=='TMM')
    #  vr <- .calcTMMvar(x, pair)
  }

  if(!is.na(vr))
    sd <- sqrt(vr/sqrt(nrow(x$counts)))
  attributes(nr) <- list(sd=sd) 
  return(nr)
}


.blncd.dvrst.ct <- function(s, bal.denom=c('max','mean')){
  bal.denom = match.arg(bal.denom)
  sum(s/eval(call(bal.denom, s, na.rm=T)))
}


.calcLibSizeNormFact <- function(x){
  ct.sum <- apply(x$counts, 2, sum)
  ct.sum/(sum(ct.sum)/2)
}


.rawR <- function(x, pair, uniques)
  raPlot(x$counts[,pair[c('obs','ref')]], jitter=FALSE, uniques=uniques, plot=FALSE)$R





.mtd <- function(x, meta.lev, pair, method=c('raw','balanced'), expnt=NA, xpnt.atr=FALSE){ #, sbst=1:nrow(x$counts)){ 
  ## currently only works with no replicates only ...
  ## ... and assumes counts columns are in the same order as first two meta taxa columns

  method <- match.arg(method)
  dvrst <- sapply(x$meta[[meta.lev]],
                  function(gt){  
                    switch(method,
                           balanced=.blncd.dvrst.ct(gt$sum, bal.denom='max'),                                           
                           raw = nrow(gt),
                           'ERROR: switch statement couldnt find your method'
                           )
                    
                  }
                  )

  
  if(nrow(x$meta.sum[[meta.lev]])==1){
    warning("only one meta sub level found. expnt has no effect")
    expnt <- 1
  }
  if(is.na(expnt))
    expnt <- .calcMTDexpnt(dvrst)

  if(xpnt.atr){
    ret <- dvrst^expnt
    attributes(ret) <- list(expnt=expnt, names=names(dvrst))  # disqualifies as vector and can't be used in nerge
    return(ret)
  }else{
    return(dvrst^expnt)
  }
}


.calcMTDexpnt <- function(dvrst, share=.5, pt.pct=.01, max.expnt=20){

    for(expnt in max.expnt:1){

      #sq <- seq(.1,1,.1) #sapply(sq, function(s)
      ct <- sum(cumsum(rev(sort(dvrst^expnt)))/sum(dvrst^expnt) < share)
      pct <- ct/length(dvrst)

      if(pct > pt.pct & ct > 10)
       break
    }
    return(expnt)
}

#meta.lev='genus_sp'; method='raw'; expnt=2; uniques=TRUE; nr=NA; sd=NA; waa=NA; mr=TRUE; ahist=TRUE
.MTDheatplot <- function(x, meta.lev='genus_sp', method='raw', expnt=NA, uniques=TRUE, nr=NA, sd=NA, waa=NA, mr=TRUE, ahist=TRUE, colrng=c('yellow','red'), pair, ...){
  
  .checkCondPair(x, pair)

  #RA <- plot(x, pch='',meta.lev=0,flat=T,rex=1, pair=pair, spine=0, uniques=uniques, main=meta.lev, border=NULL);
  RA <- raPlot(x$counts[,as.character(pair[c('obs','ref')])], jitter=.43, pch=20, rex=1,  spine=0, uniques=uniques, col=gray(.93), ...)
  if(is.null(RA$R))
    stop('uniques == FALSE left no points in RA')
  mtd <- .mtd(x, meta.lev=meta.lev, method=method, pair=pair, expnt=expnt) #, sbst=names(RA$R))

  RA.df <- nerge(list(R=RA$R, A=RA$A,  mtd=mtd))
        
  RA.df <- RA.df[order(RA.df$mtd),]
  
  cols <- colorRampPalette(colrng)(20)
  col.mtd <- log(RA.df$mtd, 3) #, expnt) 
  with(RA.df, points(A,R,col=cols[round(col.mtd/max(col.mtd)*19)+1]))

  if(!is.null(nr)){
    if(is.na(nr))
      nr <- with(RA.df, weighted.mean(R, mtd))
    abline(h=nr)
    waa <- with(RA.df, weighted.mean(A, mtd))
    if(is.na(sd))
      sd <- with(RA.df, sqrt(Hmisc::wtd.var(R, mtd)/sqrt(nrow(x$counts))))
    if(!is.null(sd))
      segments(x0=waa,x1=waa, y0=nr-sd, y1=nr+sd)
    if(mr)
      with(RA.df,abline(h=mean(R), col='gray'))
  }
  if(ahist){
    old.par <- par()
    plt <-  par('plt')
    plt[4] <- plt[3]+.025
    lims <- usr2lims()
    par(fig=plt, mar=rep(0,4), new=TRUE)
    require(plotrix)
    hst <- with(RA.df, weighted.hist(A,mtd,breaks=40, plot=F))
    plot(hst$mids, hst$counts, type='l', yaxt='n',xlim=lims$x, col='orange', bty='n', xaxt='n')
    par(old.par)
  }
  
  invisible(list(df=RA.df, nr=nr, sd=sd, waa=waa))
}



.calcTMMvar <- function (obj, pair, logratioTrim = 0.3, sumTrim = 0.05, Acutoff = -1e+10){
  
  .checkCondPair(x, pair)

  obs <- obj$counts[,pair['obs']]
  ref <- obj$counts[,pair['ref']]

  ## -- lifted verbatim from edgeR -- #
  nO <- sum(obs)
  nR <- sum(ref)
  logR <- log2((obs/nO)/(ref/nR))
  absE <- (log2(obs/nO) + log2(ref/nR))/2
  v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  n <- sum(fin)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >= loS & rank(absE) <= hiS)
  ## -- end lift -- #
  
  wv <- 2^wtd.var(logR[keep], 1/v[keep])
  log2(.abs.nf(obj,wv))
}
