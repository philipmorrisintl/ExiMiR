# File: read.exi.data.R
# 
# Author: sylvain.gubian@pmi.com
# Aim: read data content of a given Exiqon txt files and return a matrix of median signal

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

###################################################################################

read.exi.data<- function(txt.path, filenames, rm.background, verbose, header, is.galfile.ready, galref) {
	
	ret <- list()
	
	# Read first part of the files first
	for(i in 1:length(filenames)) {
		if (verbose)
			message("Reading in: ", filenames[i], " ...")
		res <- read.exi.data.single(txt.path, filenames[i], header)
		if (i == 1) {
			nbRows <- nrow(res)
			ret$mat <- matrix(NA, nbRows, length(filenames))
			ret$flags <- matrix(NA, nbRows, length(filenames))
			if (!is.galfile.ready) {
				make.gal.env.fromdf(df=res, type="TXT", verbose)	
			}
			else
			{
				# Check if GAL envi is compatible with TXT file
				if (verbose)
					message("Checking GAL environment/TXT compatibility...")
				l <- mget(ls(galref),galref)
				maxis <- lapply(l, maxi <- function(x) { return(max(x[,1])) })
				
				if ( max(unlist(maxis))> nrow(res)) {
					if (verbose) {
						message("GAL environment is not compatible with assay data.")
					}
					make.gal.env.fromdf(df=res, type="TXT", verbose)
				}
			}
		}
		if (!rm.background) {
			ret$mat[,i] <- as.numeric(res$Signal.Median)
			#ret$mat[,i] <- as.numeric(res$Background.Median)
		}
		else {
			ret$mat[,i] <- as.numeric(res$Signal.Median - res$Background.Median)
		}
		ret$flags[,i] <- as.integer(res$Flag)
	}
	res$b <- as.integer(gsub(" ", "",gsub("[^[:digit:]]", " ", res$Field)))
	res$r <- as.integer(res$Row)
	res$c <- as.integer(res$Column)
	maxRow <- max(res$Row)
	maxCol <- max(res$Column)
	indices <- bcr2i(res$b, res$c, res$r, maxCol, maxRow) 
	row.names(ret$mat) <- indices 
	row.names(ret$flags) <- indices
	colnames(ret$mat) <- filenames
	colnames(ret$flags) <- filenames
	
	# Get the perm.env
	perm.env <- get("perm.env", envir=globalenv())
	
	#if ( maxRow * maxCol * perm.env$block.row * perm.env$block.col != nrow(ret$mat))
	#{
	#	warning("The block structure of the array is currently not supported. The 'image' method of the exibatch object might give incorrect results.")
	#}
	
	#bi <- (res$b-1) %/% block.col  + 1
	#bj <- (res$b-1) %% block.col + 1
	#rindex <- res$r + (bi-1) * ret$rows
	#cindex <- res$c + (bj-1) * ret$cols
	
	m.temp <- matrix(NA, max(perm.env$rindex), max(perm.env$cindex))

	for(i in 1:length(filenames)) {
		for(j in 1:length(ret$mat[,i])) {
			m.temp[perm.env$rindex[j],perm.env$cindex[j]] <- ret$mat[j,i]
			# permute flags
		}
		#m.temp <- m.temp[,rev(1:max(perm.env$cindex))]
		m.temp <- m.temp[,(1:max(perm.env$cindex))]
		ret$mat[,i] <- as.vector(m.temp) 	
	}
	ret$realdim <- dim(m.temp)
	return(ret)
}


read.exi.data.single <- function(txt.path, filename, header) {
	path <- file.path(path.expand(txt.path),filename)
	con <- file(path, "r", blocking = TRUE)
	df <- read.table(con, header=TRUE, sep="\t", skip=header$beginDataIndex, nrows=header$nbRecords)
	close(con)
	return(df)
}



