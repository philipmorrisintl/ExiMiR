# File: make.gal.env.R
# 
# Author: sylvain.gubian@pmi.com

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#########################################################################################

make.gal.env <- function(filename=NULL, gal.path=getwd(), verbose=FALSE) {
	if (verbose)
		message(" Trying to read GAL file:",  filename)
	
	if(is.null(filename)) {
		if (verbose) {
			warning("No filename provided, looking for GAL file...")
		}
		if(is.null(gal.path)) gal.path <- getwd()
		filename <- dir(path=gal.path,pattern="\\.gal$")
		nfiles <- length(filename)
		if(nfiles == 0) stop("Cannot find GAL file, txt file will be used..")
		if(nfiles > 1) {
			filename <- filename[1]
			if (verbose) {
				warning(paste("More than one GAL file found. Reading",filename))
			}
		}
	}
	if(is.null(gal.path)) {
		if (verbose) {
			warning("gal.path is null, using working directory...")
		}
		gal.path <- getwd()
	}
	path <- file.path(path.expand(gal.path),filename)
	if (!file.exists(path)) {
		stop("Can not open GAL file, txt file will be used...")
	}
	con <- file(path, "r", blocking = FALSE)
	source = readLines(con)
	close(con)
	
	# Reading Blocks in a data.frame
	# Look for BlockCount
	block.count.index <- which(regexpr("BlockCount", source) != -1)
	block.count <- as.integer(gsub(" ", "",gsub("[^[:digit:]]", " ", source[block.count.index])))
	
	block.first.index <- which(regexpr("Block1=", source) != -1)
	blocks.list <- source[(block.first.index):(block.first.index+block.count-1)]
	second.col <- vector()
	for(i in 1:(block.count)) {
		second.col[i] <- as.integer(strsplit(blocks.list[i],split=",")[[1]][2])
	}
	block.row <-  length(unique(second.col))
	block.col <- length(which(second.col == second.col[1]))
	

	# Reading Blocks in a data.frame
	beginMapIndex <- which(regexpr("Block\tColumn\tRow\tID\tName", source) != -1)
	dataStr <- tail(source, -beginMapIndex+1)
	textConn <- textConnection(dataStr)
	df <- read.table(textConn, header=TRUE, sep="\t")
	close(textConn)
	rm(dataStr)
	rm(source)

	make.gal.env.fromdf(df, type="GAL", block.row=block.row, block.col=block.col)
}


make.gal.env.fromdf <- function(df, type=c("GAL", "TXT"), block.row=12, block.col=4, verbose=FALSE) {
	type <- match.arg(type)
	
	galenv = new.env(hash=TRUE, parent=emptyenv())
	
	if (type=="TXT") {
		message('Warning: Creating the GAL environment based on TXT file:')
		message("This might cause some functionalities not to work properly (e.g. AffyBatch image).")
		message("If a GAL file is provided, use 'make.gal.env' to create the proper GAL environment.")
		GENEIDS <- df$"Gene.ID"
		BLOCKS <- as.integer(gsub(" ", "",gsub("[^[:digit:]]", " ", df$Field)))	
	}
	else {
		GENEIDS <- df$ID
		BLOCKS <-as.integer(df$Block)
	}
	
	# Create permutation environment
	rows <- as.integer(df$Row)
	cols <- as.integer(df$Column)
	
	# Array of blocks dimension is temporarly hard coded
	perm.list <- list()
	perm.list$block.row <- block.row
	perm.list$block.col <- block.col
	
	bi <- (BLOCKS-1) %/% perm.list$block.col  + 1
	bj <- (BLOCKS-1) %% perm.list$block.col + 1
	
	nrow <- max(df$Row,na.rm=TRUE)
	ncol <- max(df$Column, na.rm=TRUE)
	
	perm.list$rindex <- rows + (bi-1) * nrow
	perm.list$cindex <- cols + (bj-1) * ncol
	perm.list$vect <- vector()
	perm.list$vect <- perm.list$rindex + (perm.list$cindex-1) * max(perm.list$rindex) 
	perm.env = new.env(hash=FALSE, parent=emptyenv())
	multiassign(names(perm.list), perm.list, envir=perm.env)
	assign("perm.env", perm.env, envir=globalenv())
	
	
	map <- list()
	ulist=which(!duplicated(GENEIDS))
	nams = as.character(df$Name[ulist])
	ids = GENEIDS[ulist]
	
	for(i in 1:length(ulist)) {
		values <- list()
		matches <- which(GENEIDS == ids[i])
		if (nams[i] == "") {
			nams[i] <- paste("HIDDEN_", ids[i], sep='')
		}
		else if (nchar(nams[i]) > 255) {
			nams[i] <- substr(nams[i], 1, 255)
		}
		n <- length(matches)
		if (type=="TXT") {
			blocks.match <- as.integer(gsub(" ", "",gsub("[^[:digit:]]", " ", df$Field[matches])))	
		}
		else {
			blocks.match <-as.integer(df$Block[matches])
		}
		columns.match <- as.integer(df$Column[matches])
		rows.match <- as.integer(df$Row[matches])
		
		#for(j in 1:n) {
		#	values[j] <- bcr2i(blocks.match[j], columns.match[j], rows.match[j], ncol, nrow)
		#}
		# User permutation environment to have the appropriate order
		
		for(j in 1:n) {
			values[j] <- perm.list$vect[bcr2i(blocks.match[j], columns.match[j], rows.match[j], ncol, nrow)]
		}
		
		tmp <- as.matrix(as.integer(values))
		mmNA <- matrix(NA, nrow(tmp), 1)
		tmp <- cbind(tmp, mmNA)
		colnames(tmp) <- c("pm", "mm")
		map[[i]] <- tmp
	}
	names(map) <- nams
	multiassign(names(map), map, galenv)
	assign("galenv", galenv, envir=globalenv())
}


i2bcr <- function(index, ncol, nrow) {
	res <- list()
	res$b <- (index-1) %/% (nrow*ncol) + 1 ;
	res$r <- (index-1-((res$b-1)*nrow*ncol)) %/% ncol + 1
	res$c <- (index-1-((res$b-1)*nrow*ncol)) %% ncol + 1
	return(res)
}

bcr2i <- function(block, col, row, ncol, nrow) {
	return((block-1)*nrow*ncol + (row-1)*ncol + (col-1) + 1)
}



