\name{make.gal.env}
\Rdversion{0.99.0}
\alias{make.gal.env}
\title{
GAL Environment Maker
}
\description{
Reads an Exiqon GAL file and creates an environment used as a hash table for the probeset mapping location
}
\usage{
make.gal.env(filename=NULL, gal.path=getwd(), verbose=FALSE)
}
\arguments{
	\item{filename}{
		Character vector. Filename of the GAL file
	}
	\item{gal.path}{
		Character vector. Path to the GAL file.
	}
	\item{verbose}{
		Logical. If \code{TRUE}, messages are shown.
	}
}
\details{
	This function is designed similarly to \code{make.cdf.env} from the \code{makecdfenv} package.
	If no filename is provided as argument, the function tries to read the first GAL file in the input path.
	The returned environment is a hash table. For every probeset name we have a matrix with 2 columns. The
	first column contains the PM locations and the second column the MM locations. For PM only chips the MM column
	will have NAs.
}
\value{
	An environment, used as a hash table.
}

\examples{
	# The folder 'Exiqon' contains a GAL file
	\dontrun{galenv <- make.gal.env(gal.path='Exiqon')}
}

\author{
Sylvain Gubian, Alain Sewer, PMP SA
}

