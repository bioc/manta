generateWeights <- function(x, w.clmn,  agg.clmn, cond.clmn, ct.clmns=NULL){
  
  wgts <- x[,w.clmn]
  agg <- x[,agg.clmn]
  cond <- x[,cond.clmn]
  if(!is.null(ct.clmns))
    wgts <- wgts * x[,ct.clmns]
  
  w <- wgts
  avg <- mean(wgts, na.rm=TRUE)
  
  ## negative log transformation for e-values
  if(avg < 0){
    w <- -log(wgts, 10)  #assumes an evalue cut off of e < 0
    w[wgts == 0] <- 400  # for exact matches
  }
  
  sums <- tab2df(as.table(by(w, list(agg, cond), sum, na.rm=T)))
  sums[is.na(sums)] <- 0
  
  avgs <- apply(sums, 2, mean)
  return(sums/avgs)

}