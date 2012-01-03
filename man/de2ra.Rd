\name{de2ra}
\alias{de2ra}
\title{create a pseudo-RA plot out of an edgeR exactTest output}
\description{
The edgeR package has a great function 'exactTest' for assessing differential expression of genes.  The output is somewhat already in the form of an RA plot. This function goes the extra step of bringing the condition unique points closer to the main body of the plot.
}
\usage{
de2ra(de, sig.cols=c('gray','black'), q=.05, method=p.adjust.methods, ...)
}
\arguments{
\item{de}{a list of named vectors.}
\item{sig.cols}{colors for specifying significance.}
\item{q}{the q value (or adjusted p) to use for differentiating significant from non significant differential expression.}
\item{method}{the p.adjust method}
\item{...}{additional parameters passed to 'plot'}
}
\seealso{
\link{raPlot}
}
\value{
an raPlot
}
\examples{

manta.path <- system.file("extdata","PapaGO-BLAST.results-diatoms.Rdata", package="manta")
load(manta.path)

de <- exactTest(x)
de2ra(de)

}
