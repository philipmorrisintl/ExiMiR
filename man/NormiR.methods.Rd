\name{NormiR.methods}
\alias{NormiR.bgcorrect.methods}
\alias{NormiR.normalize.methods}
\alias{NormiR.pmcorrect.methods}
\alias{NormiR.summary.methods}
\alias{NormiR.spikein.args}

\title{
ExiMiR functions for enumerating the names of low-level normalization methods or arguments 
}

\description{
These functions enumerate the names of methods or arguments of the low-level functions for
miRNA raw data normalization (i.e. background correction, (spike-in probe-based) normalization,
PM corrrection, summarization). They take into account how the input \code{AffyBatch} object
was created as well as the underlying array type.
}

\usage{
NormiR.bgcorrect.methods(object)
NormiR.normalize.methods(object)
NormiR.pmcorrect.methods(object)
NormiR.summary.methods()
NormiR.spikein.args()
}

\arguments{
\item{object}{
An \code{AffyBatch} object.
}
}

\value{
List of strings containing the names of the methods or arguments available for the input \code{AffyBatch} object.
}

\details{
See accompanying vignette.  
}

\author{
Sylvain.Gubian, Alain.Sewer, PMP SA
}

\seealso{
\code{NormiR},
\code{bg.correct.miR},
\code{norm.miR},
\code{summarize.miR}.
}

