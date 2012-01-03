\name{in2manta}
\alias{in2manta}
\alias{align2manta}
\alias{counts2manta}
\alias{tableMetas}
\alias{tableMetaSums}
\title{Convert a count or alignment table into a MANTA object}
\description{
This function converts a table of alignment data (BLAST for example) where each record is a read and each column is some parameter of the blast(s).  
The function will perform a cross-tabulation of the annotated reads into count data using (at the very least) an aggregation index column and a condition column. 
Alternatively, the input can be pre-tabulated counts and a corresponding annotation table.  The two tableMeta helper functions (called by the x2manta functions) are also documented here.
}
\usage{
align2manta(x, cond.clmn, agg.clmn, gene.clmns, meta.clmns, weight.clmn=NULL, tag.clmn=NULL, ...)
counts2manta(x, annotation, a.merge.clmn, agg.clmn, gene.clmns=NULL, meta.clmns=NULL, ...)
tableMetaSums(x, meta.clmns, cond.clmn=NULL, count.clmns=NULL)
tableMetas(x, agg.clmn, meta.clmns, cond.clmn=NULL, count.clmns=NULL)
}
\arguments{
\item{x}{The alignment or counts table.}
\item{cond.clmn}{A string indicating which column contains the conditions.  (only two different levels two allowed)}
\item{agg.clmn}{A string indicating which column in the annotation table contains the aggregation index.}
\item{annotation}{An annotation table with genes and/or meta information.}
\item{a.merge.clmn}{A stringn indicating which column in the annotation table on which to merge.  These should correspond to the row names of the counts table.}
\item{gene.clmns}{A vector of strings indicating which column contains gene annotation information.}
\item{meta.clmns}{A vector of strings indicating which column contains meta/taxinomic information.}
\item{weight.clmn}{A string indicating which column contains weighting information.}
\item{count.clmns}{A string indicating which column in a pre-cross-tabulated meta annotation table contains counts.  These counts are multiplied by the result of tableMetas()'s cross-tabulated counts to genereate the meta tables.}
\item{tag.clmn}{A string indicating which column corresponds to individual read names.  Only necessary when a single read has multiple annotation records in the table.}
\item{...}{additional parameters passed along to manta()}
}
\value{
A MANTA object
}
\seealso{
manta 
}
\examples{

align.path <- system.file("extdata","PapaGO-BLAST.results-diatoms.tab", package="manta")
a <- read.delim(align.path, stringsAsFactors=FALSE)
x <- align2manta(a, cond.clmn='treatment', agg.clmn='what_def',
			gene.clmns=c('what_def','kid','pathway'), 
			meta.clmns=c('family','genus_sp'))

	

cts.path <- system.file("extdata","PapaGO-BWA.counts-diatoms.tab", package="manta")
cts <- read.delim(cts.path)
cts.annot.path <- system.file("extdata","PapaGO-BWA.annot-diatoms.tab", package="manta")
cts.annot <- read.delim(cts.annot.path, stringsAsFactors=FALSE)

x <- counts2manta(cts, annotation=cts.annot,
                    a.merge.clmn='query_seq', agg.clmn='what_def', meta.clmns=c('family','genus_sp'),
                    gene.clmns=c('what_def','kid','pathway'))



}

