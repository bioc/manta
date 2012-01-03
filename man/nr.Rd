\name{nr}
\alias{nr}
\title{Print out all the normalization ratios for each subset in a specified taxinomic rank of a manta object}
\description{
.
}
\usage{
nr(obj, meta.lev, pair)
}
\arguments{
\item{obj}{A manta object.}
\item{meta.lev}{which taxinomic rank to use for subset normalization estimation.}
\item{pair}{A named vector of conditions factors}

}
\value{
table of normalization ratios
}
\seealso{
manta, meta2manta 
}
\examples{


manta.path <- system.file("extdata","PapaGO-BLAST.results-diatoms.Rdata", package="manta")
load(manta.path)
x <- calcNormFactors(x)

conditions <- factor(x=1:2, labels=c('ambient','plusFe')); names(conditions) <- c('ref','obs')
			
nr(x, meta.lev='genus_sp', pair=conditions)			

}

