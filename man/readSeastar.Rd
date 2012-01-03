\name{readSeastar}
\alias{readSeastar}
\title{Read SEAStAR output format}
\description{
This function reads Vaughn Iverson's SEAStAR tabular format output and appends headers to it.
}
\usage{
readSeastar(path,
clmn.names=c('seq_id','bit_score','read_count','abundance','abundance_per_total_cov','coverage','coverage_per_abundance','seq_len','gc_percent','ref_seq_name','seq_desc'),
clmn.class = c("character", rep("numeric", 6), "integer", "numeric", rep("character", 2)),
name.clmn='ref_seq_name', ret.df=FALSE, ct.calc=expression(abundance*seq_len), header = FALSE)
}
\arguments{
\item{path}{path to the seastar file}
\item{clmn.names}{The feild names for the absent header column}
\item{clmn.class}{The class types for each of the columns}
\item{name.clmn}{the seastar column to be used to name the rows of the dataframe}
\item{ret.df}{return the entire dataframe rather than just the calculated count column}
\item{ct.calc}{the equation to use to calculate the counts}
\item{header}{if the seastar tables have headers}
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

