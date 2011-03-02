\name{summarize.miR}
\Rdversion{0.99.0}
\alias{summarize.miR}

\title{
	ExiMiR summarization function.
}
\description{
Apply median summarization on the given \code{AffyBatch} object according to the GAL or CDF environment
}
\usage{
summarize.miR(abatch, out.type=c("ExpressionSet","data.frame"))
}

\arguments{
	\item{abatch}{
		An \code{AffyBatch} Object.
	}
	\item{out.type}{
		Character vector. By default, the output is an \code{ExpressionSet}. \code{data.frame} can be also used.
	}
	
}
\details{
	The GAL or CDF environment hash is used to gather probes median intensity values into the probeset record. 
}
\value{
	An \code{ExpressionSet} or a \code{data.frame} depending on the out.type argument.
}

\examples{
	data(galenv)
	data(GSE20122)
	abatch.spike <- norm.miR(GSE20122, figures.show=FALSE)
	eset.spike <- summarize.miR(abatch.spike)
	data.spike <- summarize.miR(abatch.spike, out.type="data.frame")
}

\author{
Sylvain.Gubian, Alain Sewer, PMP SA
}