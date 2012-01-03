.checkCondPair <- function(x, pair){

  if(!(is.character(pair) | is.factor(pair)))
    stop('pair must be a character or factor')
    
  if(!all(names(pair) %in% c('ref','obs')))
    stop("pair must be a named vector with names 'ref' and 'obs'")
    
  if(length(pair) != 2)
    stop('pair should be of length 2')
    
  if(!all(as.character(pair) %in% levels(x$samples$group)))
    stop('elements in pair should exist in groups of x')

}


.checkMetaLev <- function(x, meta.lev) {
   
  if(!'meta' %in% names(x))
    stop('could not find a meta slot in your object')

  if(!'meta.sum' %in% names(x))
    stop('could not find a meta.sum slot in your object')

  if(is.integer(meta.lev) | is.numeric(meta.lev)){
     if(meta.lev > length(x$meta.sum))
       stop('meta.lev must be an integer less than or equal to the length of x$meta.sum')
  }else{ 
    if(is.character(meta.lev)){
      if(!meta.lev %in% names(x$meta.sum))
        stop(paste('meta.lev',meta.lev,'must be present in meta(.sum)', paste(names(meta.sum), collapse=', ')))
    }else{
      stop('meta.lev must be an integer or character string corresponding to a taxinomic level in meta or meta.sum')
    }
  }
  
}

