\name{manta}
\alias{manta}
\title{Create a MANTA object}
\description{
The MANTA object contains counts, genes, library information just like a EdgeR's DGEList.  
Additionally, however, it contains 'meta' annotation (typically taxinomic classifications).  
This function converts all of listed component elements into a MANTA object.
}
\usage{
manta(counts, samples=makeSampleDF(counts), genes=NULL, meta=NULL, meta.sum=NULL, weights=NULL, norm=TRUE, disp=TRUE, ...)
}
\arguments{
\item{counts}{A numeric matrix containing the read counts.  Rows should be named by a unique gene identifier.}
\item{samples}{The experimental sample dataframe (nearly identical to the one in a DGEList object).}
\item{weights}{A numberic matrix of count weights for each count.  Should be the same dimentions as the count table.}
\item{genes}{A dataframe of genes annotations that have corresponding count data.}
\item{meta}{A taxinomic level list of gene lists of cross tabulations of taxinomic (meta) annotations for genes that have corresponding count data.}
\item{meta.sum}{A list of aggregated counts for (one per taxinomic level).}
\item{norm}{boolean specifiying if manta should automatically normalize using calcNormFactors().}
\item{disp}{boolean specifying if manta should automatically estimate the common dispersion via estimateCommonDisp().}
\item{...}{additional parameters passed to DGEList}
}
\value{
A MANTA object.
}
\seealso{
DGEList
}
\examples{

cts.path <- system.file("extdata","PapaGO-BWA.counts-diatoms.tab", package="manta")
cts <- read.delim(cts.path)
samples <- makeSampleDF(cts)

x <- manta(counts= cts, samples = samples)

}

