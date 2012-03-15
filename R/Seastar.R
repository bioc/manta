
# 1) seq_id = sequence_id
# 2) bit_score -- Bits of information in alignment to reference
# 3) read_count = abundance*seq_len -- Fractional read count
# 4) raw_abundance -- Length normalized fractional read count
# 5) fractional_abundance -- Fractional abundance of this reference out of all references (length normalized)
# 6) mean_coverage -- Mean coverage of reference by aligned reads
# 7) mean_read_length -- Mean read length of aligned reads
# 8) seq_len -- Length of reference sequence
# 9) gc_percent -- %GC content of reconstructed sequence (only calculated when FASTQ files provided, 0.0 otherwise)
# 10) ref_seq->seq_name -- Name of larger reference this sequence belongs to (used for multi-chromosome or multi-scaffold references)
# 11) ref_seq->seq_desc -- Description of above overarching reference sequence


readSeastar <- function(path, clmn.names=c('seq_id','bit_score','read_count','raw_abundance','fractional_abundance','mean_coverage','mean_read_length','seq_len','gc_percent','ref_seq_name','ref_seq_desc'),
                              clmn.class=c('character',rep('numeric',6), 'integer', 'numeric', rep('character',2)),
                        name.clmn='seq_id', ret.df=FALSE, ct.calc=expression(fractional_abundance*seq_len), header=FALSE, ...){

  if(!name.clmn %in% clmn.names)
    stop('name.clmn must be one of the columns specified in clmn.names')
  df <- read.delim(path, col.names=clmn.names, colClasses=clmn.class, stringsAsFactors=F, header=header, ...)
  rownames(df) <- df[,name.clmn]
  df$count <- with(df, eval(ct.calc))
  if(ret.df)
    df
  else
    nv(df, 'count')

}


                        
seastar2counts <- function(treat.paths, id.prefix=NA, all=TRUE, uniques=0,  ...){

  treat.names <- names(treat.paths)
  
  ssdfs <- lapply(treat.names, function(x) readSeastar(treat.paths[x], ...))
  names(ssdfs)<- treat.names
  
  df.merge <- nerge(ssdfs, all=all)
  
  if(!is.na(uniques))
    df.merge[is.na(df.merge)] <- uniques
  
  if(!is.na(id.prefix))
    rownames(df.merge) <- paste(id.prefix, rownames(df.merge), sep='')
  
  return(df.merge)
}  
