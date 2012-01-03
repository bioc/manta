\name{meta2counts}
\alias{meta2counts}
\title{Convert a manta object's meta slot data into counts}
\description{
This is a helper function for the mantaMethod
}
\usage{
meta2counts(obj, meta.lev, rm.sum=TRUE, meta.subset=NULL)
}
\arguments{
\item{obj}{A manta object.}
\item{meta.lev}{The taxinomic rank level (phylum, genus, etc) from which to pull counts.}
\item{rm.sum}{Whether or not to remove the overall 'sum' (over all conditions) column.}
\item{meta.subset}{Which sub rank level (the subset) for which to estimate the normalization. By default is does the overall normalization}
}
\value{
a count matrix
}
\seealso{
manta, mantaMethod
}
\examples{


manta.path <- system.file("extdata","PapaGO-BLAST.results-diatoms.Rdata", package="manta")
load(manta.path)

tab <- meta2counts(x, meta.lev='genus_sp', meta.subset='Pseudo-nitzschia granii')


}

