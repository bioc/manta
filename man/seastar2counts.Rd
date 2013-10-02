\name{seastar2counts}
\alias{seastar2counts}
\title{Convert seastar output to count data}
\description{
This function performs a simple merge between two different SEAStAR tables.
}
\usage{
seastar2counts(treat.paths, id.prefix=NA, all=TRUE, uniques=0, ...)

}
\arguments{
\item{treat.paths}{a named vector of strings of paths to 2 seastar tabular files.}
\item{id.prefix}{prefix to use when naming the rownames of the merged table}
\item{all}{used in the same way as the underlying merge function's all parameter.}
\item{uniques}{determines the replacement value for genes which are unique to the opposite library.}
\item{...}{additional paramters passed on to readSeastar}
}
\value{
a named count matrix or vector
}
\seealso{
nerge, merge
}
\examples{

conditions <- c('ambient','plusFe')
ss.names <- nv(paste('Pgranii-',conditions,'.seastar', sep=''), conditions)

ss.paths <- nv(system.file("extdata",ss.names, package="manta"), conditions)

dfm <- seastar2counts(ss.paths)

}

