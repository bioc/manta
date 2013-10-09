\name{readSeastar}
\alias{readSeastar}
\title{Read SEAStAR output format}
\description{
This function reads Vaughn Iverson's SEAStAR tabular format output and appends headers to it.
}
\usage{
readSeastar(path,
clmn.names=c('seq_id','bit_score','read_count','raw_abundance','fractional_abundance','mean_coverage','mean_read_length','seq_len','gc_percent','ref_seq_name','ref_seq_desc'),
clmn.class = c("character", rep("numeric", 6), "integer", "numeric", rep("character", 2)),
name.clmn='seq_id', ret.df=FALSE, ret.clmn='read_count', ct.calc=expression(raw_abundance*seq_len), header = FALSE, ...)
}
\arguments{
\item{path}{Path to the seastar file}
\item{clmn.names}{The feild names for the absent header column}
\item{clmn.class}{The class types for each of the columns}
\item{name.clmn}{the seastar column to be used to name the rows of the dataframe}
\item{ret.df}{Return the entire dataframe rather than just the calculated count column}
\item{ret.clmn}{If ret.df is FALSE, return just this column as a vector. Can be 'count' as calculated by 'ct.calc'}
\item{ct.calc}{the equation to use to calculate the counts}
\item{header}{If the seastar tables have headers}
\item{...}{Additional parameters passed on to read.delim}
}
\value{
a SEAStAR formated matrix of per-reference/contig/gene stats (including counts)
}
\seealso{
seastar2counts
}
\examples{

conditions <- c('ambient','plusFe')
ss.names <- caroline::nv(paste('Pgranii-',conditions,'.seastar', sep=''), conditions)

ss.paths <- system.file("extdata",ss.names, package="manta")

df <- readSeastar(ss.paths[1])
  
}

