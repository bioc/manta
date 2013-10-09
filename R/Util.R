cmdArgsToVariables <- function(){
  cmd.arg.pairs.eq <- sub("--","",commandArgs()[grep("=",commandArgs())])
  cmd.arg.pairs <- strsplit(cmd.arg.pairs.eq,"=")
  keys <- sapply(cmd.arg.pairs,"[",1)
  values <- sapply(cmd.arg.pairs,"[",2)
  ## convert the reserved logical types properly
  rsvd <- list('TRUE', 'FALSE', 'NA','NULL')  #NULL gets converted to NA
  if(length(keys) > 0 & length(values) > 0){ 	
    for(i in seq(along=values)){
      value <- values[i]
      for(rsv in rsvd){
      	if (values[i] == tolower(rsv) | values[i] == rsv)
          value <- as.logical(rsv)
      }
      assign(keys[i],value,envir=.GlobalEnv)      
    }
  }else{
    warning("Couldn't find any command line arguments. No variables to assign.")
  }
}
