


pplacer2manta <- function(dir='code/R-mbrust/manta/inst/extdata/pplacer',
                          group.pattern='_([[:alpha:]]+)_',
                          groups=c('coastal','costal','DCM','surface','upwelling'),
                          uk.name='unknown', ...){

  gene.names <- list.files(dir)

  require(RSQLite)
  drv <- dbDriver("SQLite") 
  con <- dbConnect(drv, paste(dir, gene.names[1], 'taxtable.db',sep='/'))
  ranks <- unique(dbReadTable(con, "ranks")$rank)
  dbDisconnect(con)

  ## begin building manta components
  counts <- matrix(0, nrow=length(gene.names), ncol=length(groups), dimnames=list(gene.names,groups))
  meta <- list()
  for(rank in  ranks[!(ranks %in% 'root')]){
    meta[[rank]] <- vector('list', length(gene.names))
    names(meta[[rank]]) <- gene.names
  }  

  message('Reading pplacer taxtable databases:\n')
  ## loop through gene directories, querying each SQLite database along the way
  for(gene in gene.names){
    cat(gene)
    cat(', ')
    con <- dbConnect(drv, paste(dir, gene, 'taxtable.db',sep='/'))
                                        #dbListTables(con)
    
    ## read in the raw placement table in order to get upper level taxinomic ranks
    genetree <- dbGetQuery(con, 'SELECT placeclass.placement_id, placeclass.rank, name, tax_name
                          FROM   placement_classifications placeclass, placement_names placenames, taxa
                          WHERE placeclass.placement_id = placenames.placement_id
                            AND placeclass.tax_id = taxa.tax_id
                          GROUP BY  placeclass.placement_id, placeclass.rank')

    ## read in best classifications in order to enable subtractions for unknown/unclassifiable taxinomic ranks
    bestclass <- dbReadTable(con, 'best_classifications')
    ranks <- unique(dbReadTable(con, "ranks")$rank)
    placenames <- dbReadTable(con, "placement_names")
    dbDisconnect(con)
    
    genetree2 <- merge(bestclass, placenames)
    genetree2$group <- m(group.pattern, genetree2$name)
    grouptotals <- table(genetree2$group)

    ## expand length of vector to group total length (filling with zeros)
    for( group in groups[!(groups %in% names(grouptotals))])
      grouptotals[group] <- 0
    
    grouptotals <- grouptotals[groups]  
    group.totals.sum <- cbind.data.frame(data.frame(as.list(grouptotals)), sum=sum(grouptotals))
    
    if(sum(group.totals.sum)==0)
      stop('group totals sum == 0.  Please fix or remove this gene')
    
    genetree$group <- m(group.pattern, genetree$name)
    genetree$name <- NULL

    ## populate counts table
    for(group in names(grouptotals))
      counts[gene, group] <- grouptotals[group]
    
    n <- length(unique(genetree$placement_id))

    ## perform a cross tabulation for each taxinomic rank
    for(this.rank in  ranks[!(ranks %in% 'root')]){
      #print(this.rank)
      genetree.sub <- subset(genetree, rank== this.rank)
      
      if(nrow(genetree.sub)!=0){
        
        tax.tab <- tab2df(sstable(genetree.sub, idx.clmns =c('tax_name', 'group')))
                               
        tax.tab.rnms <- rownames(tax.tab)

        ## fill out dataframe with zeros where there are missing groups
        for( clmn in names(group.totals.sum)[!(names(group.totals.sum) %in% names(tax.tab))])
          tax.tab[[clmn]] <- 0      
        tax.tab <- tax.tab[,names(group.totals.sum)]
        
        tax.tab <- rbind.data.frame(tax.tab, group.totals.sum - apply(tax.tab, 2, sum))
        rownames(tax.tab) <- c(tax.tab.rnms,uk.name)
        
      }else{
        tax.tab <- data.frame(group.totals.sum)
        rownames(tax.tab) <- uk.name
      }
      meta[[this.rank]][[gene]] <- tax.tab
    }  
  }
  
  ## finish by building last manta component
  meta.sum <- lapply(meta[!(names(meta) %in% 'root')], manta:::.meta2metasum, groups)
  
  manta(counts=counts, samples=makeSampleDF(counts), genes=data.frame(gene.names), meta=meta, meta.sum=meta.sum, ...)
  
}




metataxa2subcounts <- function(x, meta.lev='species', taxa.filter){
 t(sapply(x$meta[[meta.lev]], function(i) i[taxa.filter,]))
}