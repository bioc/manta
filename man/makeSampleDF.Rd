\name{makeSampleDF}
\alias{makeSampleDF}
\title{Make a Sample Dataframe for use in Initializing a MANTA object}
\description{
The sample dataframe contains a row for each sample with a factor indicating grouping and library sizes.
}
\usage{
makeSampleDF(counts, group=factor(colnames(counts)), lib.size=colSums(counts))
}
\arguments{
\item{counts}{first number}
\item{group}{a factor specifying which of each of count columns belong to each of the two conditions.}
\item{lib.size}{the sizes (cumulative counts) of each of the libraries.}
}
\value{
a sample dataframe
}
\seealso{
DGEList, manta, setLibrarySizes
}
\examples{
cts.path <- system.file("extdata","PapaGO-BWA.counts-diatoms.tab", package="manta")
cts <- read.delim(cts.path)
sdf <- makeSampleDF(counts=cts)

}

