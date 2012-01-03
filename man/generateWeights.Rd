\name{generateWeights}
\alias{generateWeights}
\title{Generate Weights}
\description{
Because the manta plot uses integer count data, many of the points overlap and hide a large portion of the data.  
This function allows one to apply a weighting scheme to jitter points out from under each other both to show the density and expose the content of their pies (if applicable).
}
\usage{
generateWeights(x, w.clmn,  agg.clmn, cond.clmn, ct.clmns=NULL)
}
\arguments{
\item{x}{an alignment dataframe}
\item{w.clmn}{a string corresponding to a column with weight data}
\item{agg.clmn}{a string corresponding to a column with an aggregation index.}
\item{cond.clmn}{a string corresponding to a column with condition factor}
\item{ct.clmns}{a string corresponding to a column with count data}
}
\value{
a 2 by n weight matrix
}
\seealso{
manta
}
\examples{
align.path <- system.file("extdata","PapaGO-BLAST.results-diatoms.tab", package="manta")
a <- read.delim(align.path, stringsAsFactors=FALSE)
w <- generateWeights(a, 'what_e_value', 'what_def', 'treatment')

}

