\name{normfact2absTMM}
\alias{normfact2absTMM}
\title{Convert a count or alignment table into a MANTA object}
\description{
This function converts a table of alignment data (BLAST for example) where each record is a read and each column is some parameter of the blast(s).  
The function will perform a cross-tabulation of the annotated reads into count data using (at the very least) an aggregation index column and a condition column. 
Alternatively, the input can be pre-tabulated counts and a corresponding annotation table.  The two tableMeta helper functions (called by the x2manta functions) are also documented here.
}
\usage{
normfact2absTMM(x, pair, f=nv(x$samples, 'norm.factors'), sums=colSums(x$counts))
}

\arguments{
\item{x}{The manta object (can be NULL if f and sums are specified independently).}
\item{pair}{A named vector indicating the condition names and the corresponding the reference or observation status of each.}
\item{f}{The Normalization factors. By default uses the norm.factors column of the sample dataframe in x.}
\item{sums}{The column sums of the counts. By default uses the column sums of the x$counts table.}
}
\value{
a scalar normalization factor
}
\seealso{
manta, mantaMethod 
}
\examples{

conditions <- nv(factor(x=1:2, labels=c('ambient','plusFe')) ,c('ref','obs'))
manta.path <- system.file("extdata","PapaGO-BLAST.results-diatoms.Rdata", package="manta")
load(manta.path)

x$mm <- normfact2absTMM(x=x, pair=conditions)

}

