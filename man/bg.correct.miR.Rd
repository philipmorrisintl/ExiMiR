\name{bg.correct.miR}
\alias{bg.correct.miR}

\title{
ExiMiR low-level function for background correction
}
\description{
This function performs background correction on an \code{AffyBatch} object.

The methods supported by \code{bg.correct.miR} are provided by the \code{affy}
or \code{limma} packages, depending on whether the input \code{AffyBatch} object
has been created with \code{ReadAffy} or \code{ReadExi}/\code{createAB}, respectively.

}

\usage{
bg.correct.miR(abatch,
	       bgcorrect.method='auto',
	       bgcorrect.param=list(),
	       verbose=FALSE)
}

\arguments{
\item{abatch}{
	An \code{AffyBatch} object.
}
\item{bgcorrect.method}{
	Character vector. It contains the name of the background correction method. Running
	\code{NormiR.bgcorrect.methods(abatch)} indicates which methods can be used, depending
        on the raw data contained in the \code{abatch} object. The \code{auto} option
	corresponds to the default choice of applying \code{rma} for single-channel arrays
	and \code{normexp} for dual-channel arrays.
}
\item{bgcorrect.param}{
	A \code{R} list containing the parameters required by the selected background correction method,
	as specified in the documentation of the original functions \code{bg.correct} of the
	\code{affy} package or \code{backgroundCorrect} of the \code{limma} package. As an illustration
	the parameters of the \code{normexp} method of the \code{limma} package are given below.
	\describe{
	\item{normexp.method}{
		Character vector. The variant of the \code{normexp} method, matching exactly
		the argument \code{normexp.method} of the \code{backgroundCorrect} function. 
	}
	\item{offset}{
		Numeric value to add to intensities. It matches exactly the argument
 		\code{offset} of the \code{backgroundCorrect} function. 
	}
	}
}
\item{verbose}{
	Logical. The default value is \code{TRUE}. The details of the function execution are displayed on the console.
}
}

\details{
See accompanying vignette.  
}

\value{
An \code{AffyBatch} object containing the background-corrected raw expression data.
}

\examples{
data(galenv)
data(GSE20122)
NormiR.bgcorrect.methods(GSE20122)
GSE20122.bgcorrected <- bg.correct.miR(GSE20122,
                                       bgcorrect.method='normexp',
                                       bgcorrect.param=list(offset=50))
}

\author{
Sylvain.Gubian, Alain Sewer, PMP SA
}

\seealso{
\code{NormiR.bgcorrect.methods},
\code{NormiR}.
}