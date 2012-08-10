# File: make.gal.env.R
# 
# Author: Sylvain Gubian, Alain Sewer, PMP SA

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#########################################################################################

make.gal.env <- function(galname=NULL, filename=NULL, gal.path=getwd(), verbose=FALSE) {
	if (is.null(galname)) {
		stop("galname argument has to be provided")
	}
	if (isTRUE(verbose)) {
		cat("Reading GAL file...\n")
	}
	if(is.null(filename)) {
		if (verbose) {
			cat(paste("No filename provided, looking for GAL file in folder:", gal.path, "...\n"))
		}
		if(is.null(gal.path)) gal.path <- getwd()
		filename <- dir(path=gal.path,pattern="\\.gal$")
		nfiles <- length(filename)
		if(nfiles == 0) stop("Cannot find GAL file.")
		if (nfiles == 1) {
			if (verbose) {
				cat(paste('Using GAL file:', filename[1],"\n"))
			}
		}
		if(nfiles > 1) {
			stop("Several GAL files detected, use 'filename' argument to specify which one to use.")
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
		stop("Can not open GAL file")
	}
	
	df <- readGAL(galfile=filename, path=gal.path)
	printer <-  getLayout(df)
	IDS <- df$ID
	GNAMES <- df$Name
	BLOCKS <- df$Block
	COLS <- df$Column
	ROWS <- df$Row
	
	perm.list <- list()
	perm.list$block.row <- printer$ngrid.r
	perm.list$block.col <- printer$ngrid.c
	bi <- (BLOCKS-1) %/% perm.list$block.col  + 1
	bj <- (BLOCKS-1) %% perm.list$block.col + 1
	NROWS <- printer$nspot.r
	NCOLS <- printer$nspot.c
	perm.list$rindex <- ROWS + (bi-1) * NROWS
	perm.list$cindex <- COLS + (bj-1) * NCOLS
	perm.list$vect <- vector()
	perm.list$vect <- perm.list$rindex + (perm.list$cindex-1) * max(perm.list$rindex)
	
	annotation.struct <- list(IDS=IDS, GNAMES=GNAMES, BLOCKS=BLOCKS, COLS=COLS, ROWS=ROWS, NROWS=NROWS, NCOLS=NCOLS)
	if (isTRUE(verbose)) {
		cat(paste("Creating environment using name:", galname,"\n"))
	}	
	create.gal.env(galname = galname, struct=annotation.struct, perm=perm.list)
	
	galfile.env <- new.env(hash=TRUE, parent=emptyenv())
	galfile.env.name <- paste(galname,".envfile", sep="")
	galfile.list <- list(path=gal.path, filename=filename)
	assign("galfile.info", galfile.list, galfile.env)
	assign(galfile.env.name, galfile.env, envir=globalenv())
}

create.gal.env <- function(galname, struct, perm) {
	# Create galenv for the affybatch
	galenv <- new.env(hash=TRUE, parent=emptyenv())
	map <- list()
	#map.bg <- list()
	ulist=which(!duplicated(struct$IDS))
	nams = as.character(struct$GNAMES[ulist])
	ids = struct$IDS[ulist]
	
	for(i in 1:length(ulist)) {
		values <- list()
		#values.bg <- list()
		matches <- which(struct$IDS == ids[i])
		if (nams[i] == "") {
			nams[i] <- paste("HIDDEN_", ids[i], sep='')
		} else if (nchar(nams[i]) > 255) {
			nams[i] <- substr(nams[i], 1, 255)
		}
		n <- length(matches)
		blocks.match <- struct$BLOCKS[matches]
		columns.match <- struct$COLS[matches]
		rows.match <- struct$ROWS[matches]
		
		if (is.null(perm)) {
			for(j in 1:n) {
				values[j] <- bcr2i(blocks.match[j], columns.match[j], rows.match[j], struct$NCOLS, struct$NROWS) 
				#values.bg[j] <- bcr2i(blocks.match[j], columns.match[j], rows.match[j], ncol, nrow) + length(ROWS)
			}
		} else {
			for(j in 1:n) {
				values[j] <- perm$vect[bcr2i(blocks.match[j], columns.match[j], rows.match[j], struct$NCOLS, struct$NROWS)]
				#values.bg[j] <- perm.list$vect[bcr2i(blocks.match[j], columns.match[j], rows.match[j], ncol, nrow)] +length(ROWS)
			}
		}
		tmp <- as.matrix(as.integer(values))
		#tmpbg <- as.matrix(as.integer(values.bg))
		mmNA <- matrix(NA, nrow(tmp), 1)
		tmp <- cbind(tmp, mmNA)
		#tmpbg <- cbind(tmpbg,mmNA)
		colnames(tmp) <- c("pm", "mm")
		#colnames(tmpbg) <- c("pm", "mm")
		map[[i]] <- tmp
		#map.bg[[i]] <- tmpbg
	}
	
	names(map) <- nams
	#names(map.bg) <- nams
	multiassign(names(map), map, galenv)
	#multiassign(names(map.bg), map.bg, galenv.bg)
	assign(galname, galenv, envir=globalenv())
	#assign(paste(galname,".bg",sep=""),  galenv.bg, envir=globalenv())
}

i2bcr <- function(index, ncols, nrows) {
	res <- list()
	res$b <- (index-1) %/% (nrows*ncols) + 1 ;
	res$r <- (index-1-((res$b-1)*nrows*ncols)) %/% ncols + 1
	res$c <- (index-1-((res$b-1)*nrows*ncols)) %% ncols + 1
	return(res)
}

bcr2i <- function(block, col, row, ncols, nrows) {
	return((block - 1)* nrows * ncols + (row-1) * ncols + (col-1) + 1)
}



