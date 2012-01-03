\name{plot.manta}
\alias{manta.ra}
\alias{plot.manta}
\title{Plot a MAnTA object}
\description{
A MANTA RAy plot is designed for visualizing comparative meta-transcriptomics, 
but can be used in any case where fold change data assessment calls for displaying additional meta information.
}
\usage{
manta.ra(x, uniques=4, pair=nv(levels(x$samples$group)[1:2] , c('ref','obs')),
                       nr=0, alpha = 0.01, normalize=FALSE, 
                       meta.level=names(x$meta.sum)[1], meta.lgnd.lim=6, lgd.pos='topright', lgd.cex=.75, lgd.trunk=FALSE, pie.lwd=1,
                       annot=NULL, vrb.axlabs=TRUE, jitter=.43, border='black',
                       rex=2, flat=FALSE, tail=.5, arms=.5, spine=1, ...)
}
\arguments{
\item{x}{a MAnTA object.}
\item{uniques}{a boolean specifying whether or not to plot the library-unique genes (those with zero counts in one or the other library).}
\item{pair}{a named vector indicating the condition names and the corresponding the reference or observation status of each.}
\item{nr}{a numeric value indicating the asymptotic normalization ratio line.}
\item{alpha}{a statistical significance value.}
\item{normalize}{A boolean specifying whether or not to normalize the counts into proportions.}
\item{meta.level}{a number or string specifies which taxonomic level (which meta list element) should be used to generate the pie charts.  0 or FALSE indicates no pies should be drawn.}
\item{meta.lgnd.lim}{an integer specifying the top n most common taxonomic levels to use in the pies and legend.}
\item{lgd.pos}{a string specifying the general position of the legend (see the legend function documentation).}
\item{lgd.trunk}{a boolean specifying if the legend should trunkate (for example) Genus species to G.species}
\item{lgd.cex}{a numeric value specifying the character expansion for the legend text.}
\item{pie.lwd}{the line thickness of the pie polygon border.}
\item{annot}{a named vector of strings that match which points should be annotated in the plot. names indicate the colors of the text.}
\item{vrb.axlabs}{a boolean indicating if verbose axis labels should be use to spell out exactly how the axes are calculated.}
\item{jitter}{whether or not or how much to jitter the counts into surrounding, non-overlapping space.}
\item{rex}{a numeric value specifying the radial expansion of the plot points.}
\item{flat}{a boolean for the radial expansion of points as a function of both R and A axes.}
\item{tail}{a numeric or boolean value indicating the line thickness of the two trailing curved significance lines of the RAy.}
\item{arms}{a numeric or boolean value indicating the line thickness of the two leading straight separator lines of the RAy.}
\item{spine}{a numeric or boolean value indicating the line thickness of the normalization line (whose y position is specified by mm).}
\item{border}{a vector of strings used to color the borders of the points.}
\item{...}{other parameters passed to plot.}

}
\value{
A MAnTA RAy plot.
}
\seealso{
maPlot, plotMA, raPlot, pies
}
\examples{

manta.path <- system.file("extdata","PapaGO-BLAST.results-diatoms.Rdata", package="manta")
load(manta.path)
plot(x, meta.lev='genus_sp')

}

