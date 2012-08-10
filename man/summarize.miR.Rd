\name{summarize.miR}
\alias{summarize.miR}

\title{
ExiMiR low-level function for summarization.
}

\description{
This function performs summarization on an \code{AffyBatch} object using a
GAL or CDF annotation environment and generates an \code{ExpressionSet} object
containing the results.
}

\usage{
summarize.miR(abatch,
	      pmcorrect.method='pmonly',
	      pmcorrect.param=list(),
	      summary.method='medianpolish',
	      summary.param=list(),
	      summary.subset=NULL)
}

\arguments{
\item{abatch}{
	An \code{AffyBatch} object.
}
\item{pmcorrect.method}{
	Character vector. It contains the name of the PM correction method, which is enabled
	only when numerical values are available for the MM probes of the input
	\code{AffyBatch} object. Running \code{NormiR.pmcorrect.methods(abatch)} indicates
	which other methods can be chosen instead of the default one \code{pmonly}.
}
\item{pmcorrect.param}{
	A \code{R} list of optional parameters for the selected \code{pmcorrect.method}, as specified
	in the documentation of the original function \code{pmcorrect} function of the \code{affy} package.
}
\item{summary.method}{
	Character vector. It contains the name of the summarization method. Running
        \code{NormiR.summary.methods()} indicates which other methods can be chosen instead
	of the default one \code{medianpolish}.
}
\item{summary.param}{
	A \code{R} list of optional parameters for the selected \code{summary.method}, as specified
	in the documentation of the original \code{AffyBatch} method \code{computeExprSet} contained in
	the \code{affy} package.
}
\item{summary.subset}{
	A \code{R} list of probe set identifiers. When left to its default \code{NULL} value, the
	summarized expression values are computed for all probe sets available on the array.
}	
}

\value{
An \code{ExpressionSet} containing the summarized expression data.
}

\examples{
data(galenv)
data(GSE20122)
eset <- summarize.miR(GSE20122,
                      summary.method="medianpolish")
}

\author{
Sylvain.Gubian, Alain Sewer, PMP SA
}

\seealso{
\code{NormiR.pmcorrect.methods},
\code{NormiR.summary.methods},
\code{NormiR}.
}

