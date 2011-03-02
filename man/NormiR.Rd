\name{NormiR}
\Rdversion{0.99.0}
\alias{NormiR}
\title{
miRNA raw data normalization function (high level)
}
\description{
	This function converts an AffyBatch object into an ExpressionSet object
	performing both normalization and summarization. By default it uses the spike-in
	probe-based normalization method and the median summarization. In case the
	spike-in probe-based method cannot be applied, a median normalization is
	executed instead. Several options allow to force the execution of the spike-in
	probe-based normalization and to fine-tune the resulting correction functions.	
}
\usage{
	NormiR(	abatch,
		method=c("spikein","mean","median"),
		background.correct=FALSE,
		verbose=TRUE,
		figures.show=TRUE,
		figures.output=c("display","file"),
		out.type=c("ExpressionSet", "data.frame"),
		min.corr=0.5,
		loess.span=-1,
		extrap.points=2,
		extrap.method=c("mean","linear"),
		force.zero=FALSE,
		cover.ext=0.5,
		cover.int=1/3,
		max.log2span=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{abatch}{
		AffyBatch object
	}
	\item{method}{
		Character vector. By default, \code{spikein} method is used. \code{mean} or \code{median} can also be selected and are used in case the 'spike-in' method can not be applied.
	}
  	\item{background.correct}{
		Logical. Default value is \code{FALSE}. If \code{TRUE}, the \code{rma} background correction is applied.
	}
 	\item{verbose}{
		Logical. Default value is \code{TRUE}; some details are provided on the console
	}

	\item{figures.show}{
		Logical. Default value is \code{TRUE}. Control figures are generated for the \code{spikein} method.
	}
	
	\item{figures.output}{
		Character vector. By default, \code{display} is used. Figures are shown to the screen. Using \code{file} generates the figures in PDF format in the working directory.
	}

	\item{out.type}{
		Character vector. Default value is \code{ExpressionSet}. The object type output by \code{NormiR}.
	}
	\item{min.corr}{
		Numeric. Default value is 0.5. Minimal allowed value for the average of the off-diagonal elements of the Pearson correlation matrix of the spike-in probeset intensities across the arrays.
	}
	\item{loess.span}{
		Numeric. Default value is -1, which corresponds to a loess smoothing neighbourhood spanning a fraction 3/(number of spike-in probesets) of the total number of points. Other positive values are allowed, see the \code{span} argument of the R \code{loess} function
	}
	\item{extrap.points}{
		Numeric. Default value is 2. The number of spike-in probesets used in the high-intensity extrapolation of the normalization correction function.
	}
	\item{extrap.method}{
		Character vector. Default value is \code{mean}. The method used for the high-intensity extrapolation of the normalization correction function.
	}
	\item{force.zero}{
		Logical. Default value is \code{FALSE}. If \code{TRUE}, it forces the normalization correction functions to have zero values at the lower end of the probe intensity range.
	}
	\item{cover.ext}{
		Numeric. Default value is 1/2. Minimal allowed relative coverage of the spike-in probesets intensities. It is computed as the ratio between the intensity range covered by the spike-in probes and the one covered by all probes on the array.
	}
	\item{cover.int}{
		Numeric. Default value is 1/3. Maximal allowed relative intensity interval between two consecutive spike-in probesets. It is computed as the largest intensity difference between two consecutive spike-in probesets divided by the overall probe intensity range.      
	}
	\item{max.log2span}{
		Numeric. Default value is 1. Gives the maximal (log2) intensity interval allowed for the probes belonging to one spike-in probeset.
	}

}
\details{
	See accompanying vignette.  
}
\value{
An \code{ExpressionSet} object or a \code{data.frame} object, depending on the \code{out.type} option
}

\examples{
	data(galenv)
	data(GSE20122)
	eset.spike <- NormiR(GSE20122)
	eset.spike	
}

\author{
Sylvain Gubian, Alain Sewer, PMP SA
}
