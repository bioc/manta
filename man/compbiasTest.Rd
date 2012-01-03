\name{compbiasTest}
\alias{compbiasTest}
\title{compositional bias test}
\description{
.
}
\usage{
compbiasTest(x, pair=nv(levels(x$samples$group)[1:2] , c('ref','obs')), meta.lev='phylum',  meta.lev.lim=min(10,nrow(x$meta.sum[[meta.lev]])))

}
\arguments{
\item{x}{a MANTA object.}
\item{pair}{named vector of the pair of conditions of interest}
\item{meta.lev}{which taxinomic level should this test be run at}
\item{meta.lev.lim}{how many underlying taxinomic levels should the analysis be limited to}
}
\value{
A DGEList object.
}
\seealso{
DGEList, manta
}
\examples{

manta.path <- system.file("extdata","PapaGO-BLAST.results-diatoms.Rdata", package="manta")
load(manta.path)
			
compbiasTest(x, meta.lev='genus_sp')

}


