\name{norm.miR}
\alias{norm.miR}

\title{
ExiMiR low-level function for miRNA raw data normalization.
}

\description{
This function performs low-level normalization on an \code{AffyBatch} object
and returns the result in a new \code{AffyBatch} object.

By default, it applies the spike-in probe-based normalization method. In case the
spike-in probe-based method cannot be applied, a median normalization is executed
instead. Several options allow however to force the execution of the spike-in probe-based
normalization and to fine-tune the resulting correction functions.	
}

\usage{
norm.miR(abatch,
         normalize.method="spikein",
         normalize.param=list(),
         verbose=TRUE,
	 ...)
}

\arguments{
\item{abatch}{
	An \code{AffyBatch} object.
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
	\code{norm.miR} and explained below.
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
		Vector of probes names that will be used as the "spike-in probes". By default, \code{norm.miR} uses the
		probes annotated as "spike-in" by Exiqon or Affymetrix.
	}
	}
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
An \code{AffyBatch} object containing the normalized (but not summarized) expression data.
}

\examples{
data(galenv)
data(GSE20122)
GSE20122.normalized <- norm.miR(GSE20122,
                                normalize.param=list(figures.show=FALSE)) 
# Apply the affy method hist on the generated AffyBatch object GSE20122.normalized
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
hist(GSE20122)
hist(GSE20122.normalized)
layout(1)
}

\author{
Sylvain.Gubian, Alain.Sewer, PMP SA
}

\seealso{
\code{NormiR.normalize.methods},
\code{NormiR.spikein.args},
\code{NormiR}.
}

