\name{metataxa2subcounts}
\alias{metataxa2subcounts}
\title{create a new (sub) count table out of a subcomponent of the metatranscriptome}
\description{
A simple and easy way to pull a subcomponent transcriptome out of the metatranscriptome.  
This one line function is useful for those wanting to just see the count data for one species or to create a new DGE or manta object on a subset of data.
}
\usage{
metataxa2subcounts(x, meta.lev='species', taxa.filter)

}
\arguments{
\item{x}{a manta object.}
\item{meta.lev}{the name of the taxinomic rank on which to subset.}
\item{taxa.filter}{the name of the taxonomy at the meta.lev rank on which to subset.}
}
\value{
a count table
}
\seealso{
manta
}
\examples{

load(system.file('extdata','PapaGO-BLAST.results-diatoms.Rdata',package='manta'))
metataxa2subcounts(x, meta.lev='species', taxa.filter='Pseudo-nitzschia granii')

}

