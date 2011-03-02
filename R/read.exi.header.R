# File: read.exi.header.R
# 
# Author: sylvain.gubian@pmi.com
# Aim: read header content of a given Exiqon txt file into an R list object

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

###############################################################################


read.exi.header <- function(txt.path, filename, verbose=TRUE) {
	
	path <- file.path(path.expand(txt.path),filename)
	if (verbose)
		cat("Reading header in file", filename,"...\n")
	con <- file(path, "r", blocking = FALSE)
	source = readLines(con)
	close(con)
	header <- list()
	
	# Get version
	
	header$version <- parseValue(source, attributeName = 'version')
	
	# Get Date
	header$date <- parseValue(source, attributeName = 'Date')
	if (verbose)
		cat("Date of first scan: ", header$date, "\n")
	
	# Get Image File
	header$imageFile <- parseValue(source, attributeName = 'Image File')
	
	# Get Page number
	header$page <- parseValue(source, attributeName = 'Page\t')
	
	# Parse Page Name
	header$pageName <- parseValue(source, attributeName = 'Page Name')
	
	# Parse Inverted
	header$inverted  <- parseValue(source, attributeName = 'Inverted')
	
	# Parse Measurement Parameters
	#header$segmentationMethod <- parseValue(source, attributeName = 'Segmentation Method')
	#header$signalLow <- parseValue(source, attributeName = 'Signal Low')
	#header$signalHigh <- parseValue(source, attributeName = 'Signal High')
	#header$backgroundLow <- parseValue(source, attributeName = 'Background Low')
	#header$backgroundHigh <- parseValue(source, attributeName = 'Background High')
	#header$backgroundBuffer <- parseValue(source, attributeName = 'Background Buffer')
	#header$backgroundWidth <- parseValue(source, attributeName = 'Background Width')
	
	
	# Parse Alert data frame
	beginAlertIndex <- which(regexpr("Begin Alerts", source) != -1) + 1
	endAlertIndex <- which(regexpr("End Alerts", source) != -1) - 1
	alertStr <- ""
	for(i in beginAlertIndex:endAlertIndex) {
		alertStr <- paste(alertStr, sub("\t{2}", "", source[i]), "\n")
	}
	textConn <- textConnection(alertStr)
	header$alerts <- read.table(textConn, header=TRUE, sep="\t")
	close(textConn)
	rm(alertStr)
	
	# Reading Blocks in a data.frame
	beginBlockIndex <- which(regexpr("Begin Field Dimensions", source) != -1) + 1
	endBlockIndex <- which(regexpr("End Field Dimensions", source) != -1) - 1
	blocksStr <- ""
	for(i in beginBlockIndex:endBlockIndex) {
		blocksStr <- paste(blocksStr, sub("\t{2}", "", source[i]), "\n")
	}
	textConn <- textConnection(blocksStr)
	header$blocks <- read.table(textConn, header=TRUE, sep="\t")
	close(textConn)
	rm(blocksStr)

	
	header$beginDataIndex <- which(regexpr("Begin Raw Data", source) != -1) 
	endDataIndex <- which(regexpr("End Raw Data", source) != -1) 
	header$nbRecords <- endDataIndex - header$beginDataIndex - 2
		
	# A Different approach
	#lineIndex <- which(regexpr("End Field Dimensions", allContent) != -1) - 1
	#strRes <- strsplit(allContent[lineIndex], "\t", fixed=TRUE)
	#tmp <- gsub("[^[:digit:]]", " ", strTmp[[1]][3])
	#spl <- strsplit(tmp, " ")[[1]]
	#header$nbBlock <- as.numeric(spl[spl != ""])
	
	# An another approach
	# Then get the block vector
	#header$blocks <- list()
	#count <-1
	#for(i in beginBlockIndex:endBlockIndex) {
	#	header$blocks[[count]] <- list()
	#		# We suppose that blocks are ordered by ID
	#	strRes <-  strsplit(allContent[i], "\t", fixed=TRUE)
	#	header$blocks[[count]]$metaRows <- strRes[[1]][4]
	#	header$blocks[[count]]$metaCols <- strRes[[1]][5]
	#	header$blocks[[count]]$rows <- strRes[[1]][6]
	#	header$blocks[[count]]$cols <- strRes[[1]][7]
	#	count <- count + 1
	#}


	#header$poorSpots <- parseValue(source, attributeName = 'Poor Spots', severity = 2)
	#header$backgroundTestedAgainstSubgridDataOnly <- parseValue(source, attributeName = 'Background tested against subgrid data only')	
	#header$signalContaminationTestConnectedToBackgroundThreshold <- parseValue(source, attributeName = 'Signal contamination test connected to background contamination threshold')
	#header$saturationFlag <- parseValue(source, attributeName = 'Saturation flag', noTab=TRUE)
	#header$negativeSpots <- parseValue(source, attributeName = 'Negative Spots', noTab=TRUE, severity = 2)
	#header$nbEmptySpots <- parseValue(source, attributeName = '# of Empty Spots', noTab=TRUE)
	#header$nbPoorSpots <- parseValue(source, attributeName = '# of Poor Spots', noTab=TRUE)
	#header$nbNegativeSpots <- parseValue(source, attributeName = '# of Negative Spots', noTab=TRUE)
	#header$nbManuallyFlaggedSpots <- parseValue(source, attributeName = '# of Manually Flagged Spots', noTab=TRUE)
	rm(source)
	return(header)

}


parseValue <- function(source, attributeName, severity=1, noTab=FALSE) {
	lineIndex <- which(regexpr(attributeName, source) != -1)
	
	checkAttributeUnit(lineIndex, attributeName, severity)
	strRes <- strsplit(source[lineIndex], "\t", fixed=TRUE)
	if (!noTab) {
		attributePosition <- which(regexpr(attributeName, strRes[[1]])!=-1) + 1
		return(strRes[[1]][attributePosition])
	}
	else {
		attributePosition <- which(regexpr(attributeName, strRes[[1]])!=-1)
		tmp <- gsub("[^[:digit:]]", " ", strRes[[1]][attributePosition])
		spl <- strsplit(tmp, " ")[[1]]
		return(as.numeric(spl[spl != ""]))
	}
}

checkAttributeUnit <- function(index, strToParse, severity=1) {
	strMsg <- ""
	error <- FALSE
	if (length(index) == 0) {
		strMsg <- paste(strToParse, "attribute not found.")
		error <- TRUE
	}
	#else if (length(index) >1) {
	#	strMsg <- paste("Too many", strToParse, "found.")
	#	error <- TRUE
	#}
	if (error) {
		if (severity == 0)
			stop(strMsg)
		else if (severity == 1)
			warning(strMsg)
		else message(strMsg)
	}
}


