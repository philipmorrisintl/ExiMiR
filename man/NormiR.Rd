\name{NormiR}
\alias{NormiR}

\title{
ExiMiR high-level function for miRNA raw data normalization.
}

\description{
This function applies a standard raw data normalization pipeline (i.e. background
correction, normalization, PM correction if needed, and summarization) on the
input \code{AffyBatch} object and returns the result in an \code{ExpressionSet} object. 

The methods supported by \code{NormiR} for the background correction are provided by
the \code{affy} or \code{limma} packages, depending on whether the input \code{AffyBatch}
object has been created with \code{ReadAffy} or \code{ReadExi}/\code{createAB},
respectively.

By default, it applies the spike-in probe-based method for the second step of
normalization. In case the spike-in probe-based method cannot be applied, a median
normalization is executed instead. Several options allow however to force the execution
of the spike-in probe-based normalization and to fine-tune the resulting correction
functions.

The next step of PM correction is enabled only when numerical values are available
for the MM probes of the input \code{AffyBatch} object. In this case the methods proposed
by \code{NormiR} are provided by the \code{affy} package.

The methods supported by \code{NormiR} for the last step of summarization are also
provided by the \code{affy} package. They do not depend on how the input \code{AffyBatch}
object has been created.
}

\usage{
NormiR(abatch,
       # background correction
       bg.correct=TRUE,
       bgcorrect.method='auto',
       bgcorrect.param=list(),
       # normalize
       normalize=TRUE,
       normalize.method='spikein',
       normalize.param=list(),
       # pm correction (enabled only when MM-values are available)
       pmcorrect.method='pmonly',
       pmcorrect.param=list(),
       # expression values
       summary.method='medianpolish',
       summary.param=list(),
       summary.subset=NULL,
       # misc.
       verbose=FALSE,
       ...)
}

\arguments{
\item{abatch}{
	An \code{AffyBatch} object.
}
\item{bg.correct}{
	Logical. Default is \code{TRUE}: the background correction step will be performed. 
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
\item{normalize}{
	Logical. Default is \code{TRUE}: the normalization step will be performed. 
}
\item{normalize.method}{
	Character vector. It contains the name of normalization method. By default, the
	\code{spikein} method is used. Running \code{NormiR.normalize.methods(abatch)} indicates
	which other methods can be chosen, depending on the raw data contained in the
	\code{abatch} object.
}		
\item{normalize.param}{
	A \code{R} list of the arguments that are used to control the \code{spikein} normalization. Running
        \code{NormiR.spikein.args()} provides a complete list of all the tunable parameters supported by
	\code{NormiR} and explained below.
	\describe{
	\item{figures.output}{
		Character vector. By default, \code{display} is used. Figures are shown to the screen. Using \code{file}
		generates the figures in PDF format in the working directory.
	}
	\item{min.corr}{
		Numeric. Default value is 0.5. Minimal allowed value for the average of the off-diagonal elements of the
		Pearson correlation matrix of the spike-in probeset intensities across the arrays.
	}
	\item{loess.span}{
		Numeric. Default value is -1, which corresponds to a loess smoothing neighbourhood spanning a fraction
		3/(number of spike-in probesets) of the total number of points. Other positive values are allowed, see
		the \code{span} argument of the R \code{loess} function
	}
	\item{extrap.points}{
		Numeric. Default value is 2. The number of spike-in probesets used in the high-intensity extrapolation of
		the normalization correction function.
	}
	\item{extrap.method}{
		Character vector. Default value is \code{mean}. The method used for the high-intensity extrapolation of
		the normalization correction function.
	}
	\item{force.zero}{
		Logical. Default value is \code{FALSE}. If \code{TRUE}, it forces the normalization correction functions to
		have zero values at the lower end of the probe intensity range.
	}
	\item{cover.ext}{
		Numeric. Default value is 1/2. Minimal allowed relative coverage of the spike-in probesets intensities. It
		is computed as the ratio between the intensity range covered by the spike-in probes and the one covered by
		all probes on the array.
	}
	\item{cover.int}{
		Numeric. Default value is 1/3. Maximal allowed relative intensity interval between two consecutive
		spike-in probesets. It is computed as the largest intensity difference between two consecutive spike-in
		probesets divided by the overall probe intensity range.
	}
	\item{verbose}{
		Logical. Default is \code{TRUE}; some details are provided on the console.
	}
	\item{max.log2span}{
		Numeric. Default value is 1. Gives the maximal (log2) intensity interval allowed for the probes belonging
		to one spike-in probeset.
	}
	\item{probeset.list}{
		Vector of probes names that will be used as the "spike-in probes". By default, \code{NormiR} uses the
		probes annotated as "spike-in" by Exiqon or Affymetrix.
	}
	}
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
	A \code{R} list of probe set identifiers. When set to its default \code{NULL} value, the
	summarized expression values are computed for all probe sets available on the array.
}
\item{verbose}{
	Logical. The default value is \code{TRUE}. The details of the function execution are displayed on the console.
}
\item{...}{
	Any additional argument. Used for backward compatibility.
}
}

\details{
See accompanying vignette.  
}

\value{
An \code{ExpressionSet} object containing the normalized expression data.
}

\examples{
data(galenv)
data(GSE20122)
eset.spike <- NormiR(GSE20122,
                     bg.correct=FALSE,
                     normalize.method='spikein',
                     summary.method='medianpolish')
}

\author{
Sylvain Gubian, Alain Sewer, PMP SA
}


\seealso{
\code{bg.correct.miR},
\code{NormiR.bgcorrect.methods},
\code{norm.miR},
\code{NormiR.normalize.methods},
\code{NormiR.spikein.args},
\code{NormiR.pmcorrect.methods},
\code{summarize.miR},
\code{NormiR.summary.methods}.
}

