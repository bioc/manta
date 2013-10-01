## Run a R batch command line version of plotRA

## run a command like the following:
##  R CMD BATCH --input.file.path=myinput.tab manta.R manta.OUT
## the tab delimted input file should have the same format as the dataframe input to plotRA:
##  'treatment, group, subgroup, weight'  (the last two columns are optional)


## convert the command line paramters into variables
library(manta)
cmdArgsToVariables() 

if(exists('input.file.path')){ 

  plot.ra.tab <- read.csv(input.file.path, as.is=T, row.names=1)

  ## clean out rows that don't have all non-missing values
  plot.ra.tab <- subset(plot.ra.tab, apply(!sapply(plot.ra.tab, is.na), 1, all))

  
  #pdf(paste(input.file.path,'.pdf',sep=''), width=9, height=12)
  png(paste(base.plot.file.name,'.png',sep=''),height=2560/2/72-5, width=1600/72)
  
  ppgoRA <- makeRA(x=plot.ra.tab, use.weights=TRUE, epsilon=.5, normalize='auto')

  func.tab.out <- ppgoRA$groups    
  write.table(func.tab.out, paste(input.file.path, '-counts+stats.tab', sep=''), sep='\t', quote=F)
  
  
  if(!exists('txt.sig'))
    txt.sig <- .05
    
  plotRA(ppgoRA, txt.sig=txt.sig, pie.sig=1, cex=.3)
  dev.off()

  

}


