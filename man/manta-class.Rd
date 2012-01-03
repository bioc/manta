\name{manta-class}
\alias{manta-class}
\docType{class}
\title{Microbial Assemblage Normalized Transcript Analysis - class}

\description{
A simple list-based class for storing read counts from digital gene expression technologies and other important information for the analysis of (meta)transcriptomic data.
}

\section{Slots/List Components}{
Objects of this class contain (at least) the following list components:
  
  \code{counts}: numeric matrix containing the read counts.
  
  \code{samples}: data.frame containing the library size and group labels.
  
Additionally the class should contain the following meta/taxinomic information list components:
    
  \code{meta}: numeric matrix containing the read counts.
    
  \code{meta.sum}: data.frame containing the library size and group labels.
  
Also, the list could also contain further optional information:

  \code{genes}: data.frame containing further gene annotation information for each row in counts.

}

\seealso{ \code{\link{manta}} }

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{manta} objects also have a \code{show} method.
}

\keyword{classes}
