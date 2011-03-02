# File: read.exi.R
# 
# Author: sylvain.gubian@pmi.com
# Aim: Read an exiqon file, maps it with the GAL file and return an exi.batch R structure

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#########################################################################################


ReadExi <- function(
		txtfile.path= getwd(),
		galname= NULL,
		description = NULL,
		notes = '',
		rm.background = FALSE,
		verbose=TRUE) {
	
	if (is.null(galname)) {
		is.galfile.ready <- FALSE
		galenv <- NULL
	}
	else if (exists(galname, inherits=FALSE, where=.GlobalEnv)) {
		galenv <- as.environment(get(galname,inherits=FALSE,envir=.GlobalEnv))
		is.galfile.ready <- TRUE
	}
	else {
		is.galfile.ready <- FALSE
	}
	
	
	# Check if the path is a directory
	if (!file.info(txtfile.path)[['isdir']])
		stop('txt.path given as argument is not folder!')
	
		
	# List files in this directory
	filenames <- list.files(txtfile.path)
	
	# Look if there is a samplesinfo.txt file
	sampleInfoIndex = which(regexpr("samplesinfo.txt", filenames) !=-1)
	
	if (length(sampleInfoIndex) == 0) {
		stop("Error: No samplesinfo.txt found!")
	}
	
	# Read the samplesinfo.txt files for having the list of files to read
	path <- file.path(path.expand(txtfile.path),'/samplesinfo.txt')
	con <- file(path, "r", blocking = FALSE)
	samplesinfo <- read.table(con, header=TRUE, sep="\t")
	close(con)
	
	# Create the list of file to read
	filenames <- as.character(samplesinfo[,2])
	filenames <- append(filenames, as.character(samplesinfo[,3]))
	
	if (is.null(description))
	{
		description <- new("MIAME")
		preproc(description)$filenames <- filenames
		preproc(description)$affyversion <- library(help=affy)$info[[2]][[2]][2]
	}
	
	if (length(notes)==0) notes(description) <- notes
	
	pdata <- data.frame(sample=1:length(filenames), row.names=filenames)
	phenoData <- new("AnnotatedDataFrame",
                     data=pdata,
                     varMetadata=data.frame(
                       labelDescription="arbitrary numbering",
                       row.names="sample"))
	
	
	
	# Read the header
	header <- read.exi.header(txtfile.path, filenames[1], verbose)
	
	scanDates <- rep(header$date, length(filenames))
	
	protocol <-
			new("AnnotatedDataFrame",
					data=data.frame("ScanDate"=scanDates, row.names=sampleNames(phenoData),
							stringsAsFactors=FALSE),
					dimLabels=c("sampleNames", "sampleColumns"))
	
	# Read the data and get the affyBatch object
	results <- read.exi.data(txtfile.path, filenames, rm.background, verbose, header, is.galfile.ready, galenv)
	
	if (verbose)
		message(paste("instantiating an AffyBatch (exprs matrix dimensions is", nrow(results$mat), "x", ncol(results$mat), ")"))
	
	if(!is.galfile.ready) {
		ab <- new("AffyBatch",
					exprs  = results$mat,
					se.exprs = results$flags,
					cdfName    = "galenv",
					phenoData  = phenoData,
					nrow       = results$realdim[1],
					ncol       = results$realdim[2],
					annotation = "galenv",
					protocolData  = phenoData[,integer(0)],
					description= description,
					notes      = notes)
	}
	else {
		ab <- new("AffyBatch",
					exprs  = results$mat,
					se.exprs = results$flags,
					cdfName    = galname,
					phenoData  = phenoData,
					nrow       =  results$realdim[1],
					ncol       = results$realdim[2],
					annotation = galname,
					protocolData  = phenoData[,integer(0)],
					description= description,
					notes      = notes)
	}

	if (verbose)
		cat("done.\n")
	
	return(ab)
}