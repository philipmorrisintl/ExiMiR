\name{make.gal.env}
\alias{make.gal.env}

\title{
GAL Environment Maker
}
\description{
Reads an Exiqon GAL file and creates an annotation environment used as a hash table for the probeset mapping location.
}
\usage{
make.gal.env(galname=NULL,
             filename=NULL,
             gal.path=getwd(),
             verbose=FALSE)
}

\arguments{
\item{galname}{
	Character vector. Name to be used for the annotation environment.
}
\item{filename}{
	Character vector. Name of the GAL file.
}
\item{gal.path}{
	Character vector. Path to the GAL file.
}
\item{verbose}{
	Logical. The default value is \code{TRUE}. The details of the function execution are displayed on the console.
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
None.
}

\examples{
# The folder 'Exiqon' contains a GAL file
\dontrun{make.gal.env(galname='galenv', gal.path='Exiqon')}
}

\author{
Sylvain Gubian, Alain Sewer, PMP SA
}
