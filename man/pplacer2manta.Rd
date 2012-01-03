\name{pplacer2manta}
\alias{pplacer2manta}
\title{convert a pplacer taxinomic placement repository to a MANTA object}
\description{
This function creates a single manta object by traversing a directory of directories of pplacer SQLite taxonomy databases (where each database called "taxtable.db" resides in a sub folder named by a single gene locus).
}
\usage{
pplacer2manta(dir, group.pattern='_([[:alpha:]]+)_', 
 		groups=c('coastal','costal','DCM','surface','upwelling'),
		uk.name='unknown', ...)

}
\arguments{
\item{dir}{directory of gene directories each with a single gene taxonomic read placement database.}
\item{group.pattern}{the regular expression used to parse out the name of the different conditions/groups from the read names.}
\item{groups}{the actual list of groups.}
\item{uk.name}{the rowname to appear in each othe meta tabulation tables for the unplaceable reads.}
\item{...}{additional parameters passed along to manta()}
}
\value{
A MANTA object
}
\seealso{
manta, in2manta
}
\examples{

KOG.SQLite.repo <- system.file('extdata','pplacer',package='manta')
pplacer2manta(dir=KOG.SQLite.repo, 
 	      groups=c('coastal','costal','DCM','surface','upwelling'),
 	      norm=FALSE, disp=FALSE
 	     )

}

