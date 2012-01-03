\name{summary.manta}
\alias{summary.manta}
\title{Summarize a MANTA object}
\description{
Currently this funciton merely dumps the contents of the meta.sums tables to screen if available.
}
\usage{
summary.manta(object, ...)
}
\arguments{
\item{object}{a MANTA object.}
\item{...}{additinal arguments}
}
\value{
A MANTA summary printout.
}
\seealso{
manta
}
\examples{

manta.path <- system.file("extdata","PapaGO-BLAST.results-diatoms.Rdata", package="manta")
load(manta.path)
summary(x)

}

