\name{nf2nr}
\alias{nf2nr}
\title{convert the normalization factors to a normalization line}
\description{
.
}
\usage{
nf2nr(x, pair, method='nf', absolute=TRUE)
}
\arguments{
\item{x}{a MANTA object.}
\item{pair}{The two dimentions (conditions) which you wish to convert.}
\item{method}{.}
\item{absolute}{.}

}
\value{
a scalar normalization ratio
}
\seealso{
DGEList, manta, calcNormFactors
}
\examples{

manta.path <- system.file("extdata","PapaGO-BLAST.results-diatoms.Rdata", package="manta")
load(manta.path)

nf2nr(x)


}


