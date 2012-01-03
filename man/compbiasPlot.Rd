\name{compbiasPlot}
\alias{compbiasPlot}
\title{plot the compositional bias for the sub-taxinomic rank}
\description{
.
}
\usage{
compbiasPlot(x, pair=nv(levels(x$samples$group)[1:2] , c('ref','obs')), meta.lev='phylum',  meta.lev.lim=nrow(x$meta.sum[[meta.lev]]), breaks=10, xlim=c(0,.5), type='box', col=NULL, ...)
}
\arguments{
\item{x}{a MANTA object.}
\item{pair}{named vector of the pair of conditions of interest}
\item{meta.lev}{the taxinomic rank to plot.}
\item{meta.lev.lim}{the taxinomic rank for diversity measure.}
\item{breaks}{breaks for the histogram.}
\item{xlim}{x axis limits.}
\item{type}{type of plot.}
\item{col}{colors}
\item{...}{additional params.}
}
\value{
a compositional bias plot
}
\seealso{
DGEList, manta
}
\examples{

manta.path <- system.file("extdata","PapaGO-BLAST.results-diatoms.Rdata", package="manta")
load(manta.path)
			
compbiasPlot(x, meta.lev='genus_sp')

}


