\name{ReadExi}
\alias{ReadExi}

\title{
ExiMiR function for reading Exiqon raw data into an AffyBatch object.
}

\description{
This function reads Exiqon raw data in ImaGene file format and creates an \code{AffyBatch}
object.
}

\usage{
ReadExi(txtfile.path=getwd(),
        galname=NULL,
        description=NULL,
	notes='',
	rm.background=FALSE,
	verbose=TRUE)
}

\arguments{
\item{txtfile.path}{
	Character vector. It contains the path to the folder containing the
	\code{samplesinfo.txt} file and the Exiqon raw data files in ImaGene \code{txt} format.
}
\item{galname}{
	Character vector. The default value is \code{NULL}. In this case the GAL annotation
	environment used by \code{ReadExi} for generating the resulting \code{AffyBatch}
	object is lost. Assigning \code{galname} a non-empty value allows to control this GAL
	environment, which is useful in two specific situations. First, it gives a handle to this
	GAL annotation environment for later use. Second, if an adequate GAL annotation
	environment already exists in the memory (e.g. after having been generated by
	\code{ReadExi} or by \code{make.gal.env}), \code{galname} allows to force
	\code{ReadExi} to use it for generating the resulting \code{AffyBatch} object.
}
\item{description}{
	Object of class MIAME, as specified in the documentation of the \code{AffyBatch}
	object provided in the \code{affy} package.
}
\item{notes}{
	Character vector, as specified in the documentation of the \code{AffyBatch} object
	provided in the \code{affy} package.
}
\item{rm.background}{
	Logical. This option is kept for compatibility reasons but it is not used anymore. See
	the \code{NormiR} options for background correction.
}
\item{verbose}{
	Logical. The default value is \code{TRUE}. The details of the function execution are displayed on the console.
}
}

\details{
The Exiqon miRNA raw expression data are normally in ImageGene \code{txt} file format and
accompanied by a \code{samplesinfo.txt} description file. It enumerates the names of the
sample files for each channel. Therefore the \code{txtfile.path} argument of \code{ReadExi}
must be a folder that contains the ImageGene and the \code{samplesinfo.txt} files. If this
is not the case, \code{ReadExi} stops.

The \code{galname} argument of \code{ReadExi} must be the name of a GAL annotation
environment created with the \code{make.gal.env} or the \code{ReadExi} functions. If \code{galname}
is provided a \code{NULL} value, which is the default situation, a minimal GAL annotation
environment is created based on the annotation contained in the ImageGene \code{txt} files.
}

\value{
An \code{AffyBatch} object containing the raw expression data.
}

\examples{
# The folder 'Exiqon' contains the file 'samplesinfo.txt' and  the corresponding raw data files in ImaGene format
\dontrun{ebatch <- ReadExi(txtfile.path='Exiqon')}
# If the GAL environment has already created by the function make.gal.env
\dontrun{ebatch <- ReadExi(galenv='galenv', txtfile.path='Exiqon')}
}


\author{
Sylvain Gubian, Alain Sewer, PMP SA
}

\section{Warning}{
The \code{image} method of the \code{AffyBatch} object might not work properly when the
\code{galname} argument of \code{ReadExi} has not been assigned.
}

\seealso{
\code{make.gal.env},
\code{createAB}.
}
