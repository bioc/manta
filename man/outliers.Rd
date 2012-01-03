
\name{outGenes}
\alias{outGenes}
\title{find the most significant or highest fold change outlier genes}
\description{
.
}
\usage{
outGenes(x, n=50, p=1e-5, FC=1, A.pct=.05, uk.filter=NULL, method='BH',  verbose=TRUE)
}
\arguments{
\item{x}{a MANTA object.}
\item{n}{number of genes to return.}
\item{p}{adjusted p-value cut-off.}
\item{FC}{fold change cut-off.}
\item{A.pct}{percentage of genes that should be right most outliers.}
\item{uk.filter}{a character vector of unknown gene name paterns to filter out (eg 'lcl').}
\item{method}{multiple testing correction method from p.adjust}
\item{verbose}{should the function print these results (or just return the table)}
}
\value{
a table of the outlier genes
}
\seealso{
topTags
}
\examples{


manta.path <- system.file("extdata","PapaGO-BLAST.results-diatoms.Rdata", package="manta")
load(manta.path)

de <- exactTest(x)
outGenes(de)

}
