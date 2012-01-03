
# 1) seq_id
# 2) bit_score -- Bits of information in alignment to reference
# 3) abundance*seq_len -- Fractional read count
# 4) abundance -- Length normalized fractional read count
# 5) abundance/total_ab -- Fractional abundance of this reference out of all references (length normalized)
# 6) coverage -- Mean coverage of reference by aligned reads
# 7) coverage/abundance -- Mean read length of aligned reads
# 8) seq_len -- Length of reference sequence
# 9) gc_percent -- %GC content of reconstructed sequence (only calculated when FASTQ files provided, 0.0 otherwise)
# 10) ref_seq->seq_name -- Name of larger reference this sequence belongs to (used for multi-chromosome or multi-scaffold references)
# 11) ref_seq->seq_desc -- Description of above overarching reference sequence


readSeastar <- function(path, clmn.names=c('seq_id','bit_score','read_count','abundance','abundance_per_total_cov','coverage','coverage_per_abundance','seq_len','gc_percent','ref_seq_name','seq_desc'),
                              clmn.class=c('character',rep('numeric',6), 'integer', 'numeric', rep('character',2)),
                        name.clmn='ref_seq_name', ret.df=FALSE, ct.calc=expression(abundance*seq_len), header=FALSE){

  if(!name.clmn %in% clmn.names)
    stop('name.clmn must be one of the columns specified in clmn.names')
  df <- read.delim(path, col.names=clmn.names, colClasses=clmn.class, stringsAsFactors=F, header=header)
  rownames(df) <- df[,name.clmn]
  df$count <- with(df, eval(ct.calc))
  if(ret.df)
    df
  else
    nv(df, 'count')

}


                        
seastar2counts <- function(treat.paths, by.col='ref_seq_name', id.prefix=NA, all=TRUE, uniques=0, header=FALSE){ #lcl|

  treat.names <- names(treat.paths)
  
  ssdfs <- lapply(treat.names, function(x) readSeastar(treat.paths[x], header=header))
  names(ssdfs)<- treat.names
  
  df.merge <- nerge(ssdfs, all=T)
  
  if(!is.na(uniques))
    df.merge[is.na(df.merge)] <- uniques
  
  if(!is.na(id.prefix))
    rownames(df.merge) <- paste(id.prefix, rownames(df.merge), sep='')
  
  return(df.merge)
}  
