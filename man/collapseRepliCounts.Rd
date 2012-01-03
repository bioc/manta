
\name{collapseRepliCounts}
\alias{collapseRepliCounts}
\title{Collapse multiple technical replicate count columns into two columns for plotting}
\description{
.
}
\usage{
collapseRepliCounts(x, pair=nv(levels(x$samples$group)[1:2] , c('ref','obs')))
}
\arguments{
\item{x}{a MANTA object.}
\item{pair}{the pairs indicating which columns to collapse}
}
\value{
a collapsed, two-column count table
}
\seealso{
DGEList, manta
}
\examples{


cts <- matrix(data=rnbinom(28,2,.4), ncol=4, nrow=7)

colnames(cts) <- apply(expand.grid(c('a','b'),1:2), 1, paste, collapse='_')

x <- manta(cts, makeSampleDF(cts, group=rep(c('a','b'),2)))
collapseRepliCounts(x, pair=c('a','b'))



}


