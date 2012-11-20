# File: NormiR.R
# 
# Author: Sylvain Gubian, Alain Sewer, PMP SA
# Aim: Set of functions for Exiqon data normalization

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#########################################################################################

NormiR <- function(
		abatch,
		# background correction
		bg.correct = TRUE,
		bgcorrect.method = 'auto',
		bgcorrect.param = list(),
		# normalize
		normalize = TRUE,
		normalize.method = 'spikein',
		normalize.param = list(),
		# pm correction
		pmcorrect.method = 'pmonly',
		pmcorrect.param = list(),
		# expression values
		summary.method = 'medianpolish',
		summary.param = list(),
		summary.subset = NULL,
		# misc.
		verbose = FALSE,
		...
) {
	if (is.null(summary.method)) {
		stop("Please specify a summarization method, for example: 'medianpolish'")
	}
	
	if (is.null(pmcorrect.method)) {
		if (!is.from.createAB(abatch)) {
			stop("Please specify a pmcorrect method, for example: 'pmonly'")
		}
	}
	# For backward compatibility
	#dot.args <- match.call(expand.dots=TRUE)
	dot.args <- list(...)
	if (!is.null(dot.args$background.correct)) {
		bg.correct = dot.args$background.correct
		if (bg.correct) {
			bgcorrect.method = "auto" ;
		}
	}
	if (is.null(bgcorrect.method)) {
		bgcorrect.method = "auto" ;
	}
	
	# End of backward compatibility section
	# Call background correction if asked
	if (isTRUE(bg.correct)) {
		if(isTRUE(verbose)) {
			cat("Background correction...\n")
		}
		ab.bg.corrected <- bg.correct.miR(abatch = abatch, bgcorrect.method=bgcorrect.method, bgcorrect.param=bgcorrect.param, verbose=verbose)
	} else {
		ab.bg.corrected <- abatch
	}
	if (isTRUE(normalize)) {
		if(isTRUE(verbose)) {
			cat("Normalization...\n")
		}
		ab.normalized <- norm.miR(abatch = ab.bg.corrected, normalize.method=normalize.method, normalize.param=normalize.param, verbose = verbose,...)
	} else {
		ab.normalized <- ab.bg.corrected
	}
	
	if(isTRUE(verbose)) {
		cat("Summarization...\n")
	}
	ab.summarized <- summarize.miR(abatch = ab.normalized, pmcorrect.method = pmcorrect.method, pmcorrect.param = pmcorrect.param, summary.method = summary.method, summary.param = summary.param, summary.subset = summary.subset)
	
	return(ab.summarized)
}

bg.correct.miR <- function(abatch, bgcorrect.method="auto", bgcorrect.param = list(), verbose=FALSE) {	
	if (is.from.createAB(abatch)) {
		if (is.null(bgcorrect.param$offset)) {
			bgcorrect.param$offset = 0
		}
		# Creating the limma object for backgroundCorrect function:
		if (is.dual(abatch)) {
			if (bgcorrect.method=="auto")
				bgcorrect.method="normexp"
			if (is.null(bgcorrect.param$normexp.method)) {
				bgcorrect.param$normexp.method = "saddle"
			}
			
			lim.obj <- new("RGList")
			lim.obj$G <- exprs(abatch)[,1:(ncol(exprs(abatch))/2)]
			lim.obj$R <- exprs(abatch)[,((ncol(exprs(abatch))/2)+1):(ncol(exprs(abatch)))]
			if (has.bg(abatch)) {
				lim.obj$Gb <- se.exprs(abatch)[,1:(ncol(exprs(abatch))/2)]
				lim.obj$Rb <- se.exprs(abatch)[,((ncol(exprs(abatch))/2)+1):(ncol(exprs(abatch)))]
			} else {
				lim.obj$Gb <- NULL
				lim.obj$Rb <- NULL
			}
			lim.obj$R <- backgroundCorrect.matrix(lim.obj$R, lim.obj$Rb, method = bgcorrect.method, 
					offset = bgcorrect.param$offset, normexp.method = bgcorrect.param$normexp.method, 
					verbose = verbose)
			lim.obj$G <- backgroundCorrect.matrix(lim.obj$G, lim.obj$Gb, method = bgcorrect.method, 
					offset = bgcorrect.param$offset, normexp.method = bgcorrect.param$normexp.method, 
					verbose = verbose)
			abatch.bg.corrected <- abatch
			exprs(abatch.bg.corrected) <- cbind(lim.obj$G, lim.obj$R)
			if (has.bg(abatch))
				se.exprs(abatch.bg.corrected) <- cbind(lim.obj$Gb, lim.obj$Rb)
		} else {
#			E <- exprs(abatch)[,1:(ncol(exprs(abatch)))]
#			Eb <- NULL
#			if (has.bg(abatch)) {
#				Eb <- se.exprs(abatch)[,1:(ncol(exprs(abatch)))]
#			}
#			E <-  backgroundCorrect.matrix(E=E, Eb=Eb, method=bgcorrect.method,
#					normexp.method= bgcorrect.param$normexp.method,
#					offset=bgcorrect.param$offset,
#					verbose=verbose)
#			abatch.bg.corrected <- abatch
#			exprs(abatch.bg.corrected) <- E
#			if (has.bg(abatch))
#				se.exprs(abatch.bg.corrected) <- Eb
			if (bgcorrect.method=="auto")
				bgcorrect.method="rma"
			if ((bgcorrect.method%in%bgcorrect.methods()))
				abatch.bg.corrected <- bg.correct(abatch, method=bgcorrect.method)
			else {
				stop("Background method not supported for the AffyBatch object provided")
			}
			abatch.bg.corrected <- bg.correct(abatch, method=bgcorrect.method)
		}	
	} else {
		# affy bg.correct
		if (bgcorrect.method=="auto")
			bgcorrect.method="rma"
		if ((bgcorrect.method%in%bgcorrect.methods()))
			abatch.bg.corrected <- bg.correct(abatch, method=bgcorrect.method)
		else {
			stop("Background method not supported for the AffyBatch object provided")
		}
	}
	return(abatch.bg.corrected)
}


norm.miR <- function(abatch, normalize.method ="spikein", normalize.param = list(), verbose=TRUE, ...)
{
	# For backward compatibility
	#dot.args <- match.call(expand.dots=TRUE)
	dot.args <- list(...)
	normalize.method.params <- c(
			"min.corr",
			"loess.span",
			"extrap.points",
			"extrap.method",
			"force.zero",
			"cover.ext",
			"cover.int",
			"max.log2span",
			"figures.show",
			"figures.output")
	if (!is.null(dot.args$method)) {
		normalize.method = dot.args$method;
	}
	for(param in normalize.method.params) {
		if (!is.null(unlist(as.list(dot.args[param])))) {
			normalize.param[param] = dot.args[param];
		}
	}
	if (is.null(normalize.method)) {
		normalize.method = "spikein"
	}
	
	if (normalize.method =="spikein") {
		# Check parameters
		control.params <- list(
				min.corr = 0.5,
				figures.show=TRUE,
				figures.output="display",
				loess.span=-1,
				extrap.points=2,
				extrap.method="mean",
				force.zero=FALSE,
				cover.ext=0.5,
				cover.int=1/3,
				max.log2span=1,
				probeset.list=NULL
		)
		nmsC <- names(control.params)
		control.params[(namc <- names(normalize.param))] <-  normalize.param
		if(length(noNms <- namc[!namc %in% nmsC]))
			warning("Not supported spikein normalization parameters: ", paste(noNms,collapse=", "))
		
		galenv <- getCdfInfo(abatch)
		l <- mget(ls(galenv),galenv)
		
		have.probe.control <- FALSE
		if (!is.null(normalize.param$probeset.list)) {
			if (length(normalize.param$probeset.list)==0) {
				warning("No probe names given as control, using spikes-in controls if any.")	
			}  else {
				# Check if the list matches the probes names we have
				have.probe.control <- TRUE
				#print(normalize.param$probeset.list)
				matches <- match(normalize.param$probeset.list,names(l))
				if (any(is.na(matches))) {
					if (all(is.na(matches))) {
						have.probe.control <- FALSE
						stop("No probe set names given as control match probes names in the AffyBatch object.")
					} else {
						warning("Some probe set names given as control do not match the ones in the Affybatch object.")
					}
				}
				spike.list <- matches[!is.na(matches)]
			}
		}
		
		if (!is.from.createAB(abatch)) {
			# This does not come from createAB, its an original affybatch
			if (!have.probe.control)
				spike.list <-  which(regexpr("^spike_in-control", names(l)) != -1)
			
			if (length(spike.list)==0)
				stop("Can not find any control probes. Please use another normalization method")
			
			# Call the low level spikeinnorm function
			m <- spikeinnorm(abatch, spike.list=spike.list, control.params=control.params, verbose=verbose, channel=0)
		} else {
			if (!have.probe.control)
				spike.list <-  which(regexpr("^spike_control", names(l)) != -1)
			if (length(spike.list)==0)
				stop("Can not find any control probes. Please use another normalization method")
			# If it is dual channel
			if (is.dual(abatch)) {
				if (verbose) {
					cat("Processing first channel...\n")
				}
				m1 <- spikeinnorm(object=abatch, spike.list=spike.list, control.params=control.params, verbose=verbose, channel = 1)
				if (verbose) {
					cat("Processing second channel...\n")
				}
				m2 <- spikeinnorm(object=abatch, spike.list=spike.list, control.params=control.params, verbose=verbose, channel = 2)
				m <- cbind(m1, m2)
			} else {
				# if it is single channel
				m <- spikeinnorm(object=abatch, spike.list=spike.list, control.params=control.params, verbose=verbose, channel = 0)
			}
		}
	} else {
		if (is.from.createAB(abatch)) {
			if (!normalize.method%in%c("quantile","median","mean")) {
				stop("Method not supported for the AffyBatch object provided")
			}
			if (normalize.method=="quantile") {
				if (is.dual(abatch)) {
					lim.obj <- new("RGList")
					lim.obj$G <- exprs(abatch)[,1:(ncol(exprs(abatch))/2)]
					lim.obj$R <- exprs(abatch)[,((ncol(exprs(abatch))/2)+1):(ncol(exprs(abatch)))]
#				if (has.bg(abatch)) {
#					lim.obj$Gb <- se.exprs(abatch)[,1:(ncol(exprs(abatch))/2)]
#					lim.obj$Rb <- se.exprs(abatch)[,((ncol(exprs(abatch))/2)+1):(ncol(exprs(abatch)))]
#				}
					#lim.MA <- normalizeWithinArrays(object = lim.obj, method = "loess")
					#lim.obj$R <- normalizeQuantiles(lim.MA$A + lim.MA$M/2)
					#lim.obj$G <- normalizeQuantiles(lim.MA$A - lim.MA$M/2)
					lim.obj$R <- 2 ^ normalizeQuantiles(log2(lim.obj$R))
					lim.obj$G <- 2 ^ normalizeQuantiles(log2(lim.obj$G))
					m <- cbind(lim.obj$G, lim.obj$R)		
				} else {
					# Single channel
					lim.obj <- new("EListRaw")
					lim.obj$E <- exprs(abatch)[,1:(ncol(exprs(abatch)))]
#				lim.obj$Eb <- NULL
#				if (has.bg(abatch)) {
#					lim.obj$Eb <- se.exprs(abatch)[,1:(ncol(exprs(abatch)))]
#				}
					m <- 2 ^ normalizeQuantiles(log2(lim.obj$E))
				}
			} else {
				if (normalize.method=="median") {
					m <- mediannorm(abatch)
					
				} else if (normalize.method=="mean") {
					m <- meannorm(abatch)
				}
			}
		} else {
#			if (!normalize.method%in%normalize.AffyBatch.methods()) {
#				stop("Method not supported for the AffyBatch object provided")
#			}
			# The abatch comes from Affy pipeline
			m <- exprs(normalize(abatch, method="quantiles"))
		}
	}
	ab <- abatch
	exprs(ab) <- m
	return(ab)
}


spikeinnorm <- function(	object,
		spike.list,
		control.params,
		channel,
		verbose=TRUE) {
	# Spike-in method implementation selects only spikes in the matrix
	galenv <- getCdfInfo(object)
	
	# log2 should not be done because it is done when calling background correction even if option = FALSE
	M <- log2(exprs(object))
	
	if (channel==0) {
		sample.start <- 1;
		sample.end <- length(sampleNames(object))
		channel.name <- "1"
	}
	else if (channel==1) {
		sample.start <- 1;
		sample.end <- length(sampleNames(object)) / 2
		channel.name <- "Green"
	}
	else if (channel==2) {
		sample.start <- length(sampleNames(object)) / 2 + 1;
		sample.end <- length(sampleNames(object))
		channel.name <- "Red"
	}
	
	nbsamples <- (sample.end - sample.start) + 1
	if (nbsamples==1) {
		stop("Can not do the spike calibration with only one sample.")
	}
	l <- mget(ls(galenv),galenv)
	
	# Get the median spike control matrix
	mss <- matrix(0, length(spike.list), nbsamples)
	mss.max <- matrix(0, length(spike.list), nbsamples)
	mss.min <- matrix(0, length(spike.list), nbsamples)
	
	spikeCounts <- vector()
	for(i in 1:length(spike.list)) {
		mat <- l[spike.list[i]]
		mat <- mat[[1]]
		
		for(j in 1:nbsamples) {
			if (j==1) {
				spikeCounts[i] <- length(M[mat[,1],(j+sample.start-1)])
			}
			mss[i,j] <- median(M[mat[,1],(j+sample.start-1)])
			mss.max[i,j] <- max(M[mat[,1],(j+sample.start-1)])
			mss.min[i,j] <- min(M[mat[,1],(j+sample.start-1)])
		}
	}
	# Check mean of coeficients (except diagonal) = Shared variance check
	cmat <- cor(t(mss))
	for(i in 1:nrow(cmat)) {
		cmat[i,i] <- 0
	}
	if (mean(cmat) < control.params$min.corr) {
		if (verbose) {
			mean.display <- floor(mean(cmat) * 100) / 100
			message("Not enough common variance to guarantee good performance of spike-in normalization (", mean.display, " < ", control.params$min.corr, ").\nUsing median normalization...")
		}
		# Call mediannorm
		return(mediannorm(object,channel=channel))
	}
	
	# Checking specificity of spike intensity
	#	vs <- vector()
	#	for(i in 1:nbsamples) {
	#		vs[i] <- (max(M[,i+sample.start-1])-min(M[,i+sample.start-1])) / length(spikeCounts)
	#	}
	if ( median(mss.max-mss.min) > control.params$max.log2span) {
		if (verbose) {
			v.max <- floor(control.params$max.log2span * 100) / 100 
			v.observed <- floor(median(mss.max-mss.min) * 100) / 100
			message("The intensity resolution of the spike-in probe sets is too coarse (",  v.observed, " > ", v.max, ") to guarantee a good performance of spike-in normalization\nUsing median normalization...")
		}
		# Call mediannorm
		return(mediannorm(object, channel=channel))
	}
	
	# Checking external coverage
	ve <- vector()
	for(i in 1:nbsamples) {
		ve[i] <- (max(mss[,i])-min(mss[,i])) / (max(M[,i+sample.start-1])-min(M[,i+sample.start-1]))
	}
	if (median(ve) < control.params$cover.ext) {
		if (verbose) {
			ve.display <- floor(median(ve) * 100) / 100 
			message("The ratio between the intensity range covered by the spike-in probes and the one covered by all probes on the array is too low (",  ve.display, " < ", control.params$cover.ext, ").\nUsing median normalization...")
		}
		# Call mediannorm
		return(mediannorm(object,channel=channel))
	}
	
	# Checking internal coverage
	vi <- vector()
	for(i in 2:nbsamples) {
		vi[i-1] <- max(mss[,i]-mss[,i-1])
	}
	if (median(vi) > control.params$cover.int) {
		if (verbose) {
			vi.display <- floor(median(vi) * 100) / 100 
			message("The size of the largest intensity interval between two consecutive spikes is too large (",  vi.display, " > ", control.params$cover.int, ").\nUsing median normalization...")
		}
		# Call mediannorm
		return(mediannorm(object,channel=channel))
	}
	
	# Use shorts names for spikes control
	if (is.from.createAB(object)) {
		tab <- strsplit(names(l[spike.list]), '_')
	}
	else {
		tab <- strsplit(names(l[spike.list]), '-')
	}
	spikeNames <- list()
	for(i in 1:length(tab)) { spikeNames[i] = tab[[i]][length(tab[[i]])] }
	rownames(mss) <- spikeNames
	colnames(mss) <- (sample.start):(sample.end)
	
	mu <- rowMeans(mss)
	sigma <- apply(mss, 1, sd)
	sigma <- sqrt((ncol(mss)-1)/ncol(mss)) * sigma # We want denominator n not n-1
	
	nss <- (mss - mu) / (sqrt(ncol(mss))*sigma)
	ns <-  colMeans(nss)
	ncss <- matrix(NA, nrow(mss), ncol(mss))
	for(i in 1:ncol(mss)) {ncss[,i] <- nss[,i] - ns[i]}
	mcss <- mu + sqrt(ncol(mss))*sigma * ncss
	
	msp <- vector()
	dmsp <- vector()
	for(s in 1:nbsamples) {
		tmp <- vector()
		dtmp <- vector()
		for(i in 1:length(spike.list)) {
			mat <- l[spike.list[i]]
			mat <- mat[[1]]
			if (i==1) {
				tmp <- as.vector(M[mat[,1],s + sample.start-1])
				dtmp <-as.vector(replicate(length(M[mat[,1],s + sample.start-1]), mcss[i,s]-mss[i,s]))
			}
			else {
				tmp <- append(tmp, as.vector(M[mat[,1],s + sample.start-1]))
				dtmp <- append(dtmp, as.vector(replicate(length(M[mat[,1],s + sample.start-1]), mcss[i,s]-mss[i,s])))
			}
		}
		msp <- cbind(msp, tmp)
		dmsp <- cbind(dmsp, dtmp)
	}
	colnames(msp) <- (sample.start):(sample.end)
	colnames(dmsp) <- (sample.start):(sample.end)
	
	mcsp <- msp + dmsp
	dM <- matrix(NA, nrow(M), nbsamples)
	Mout <- matrix(NA, nrow(M), nbsamples)
	
	if (control.params$figures.show) {
		ordered <- order(mss[,1])
		txt <- paste("Fig.1: Correction of spike-in probe set intensities (", channel.name, " channel)", sep="")
		if (control.params$figures.output=="display") {
			dev.new(title=txt)
		} else {
			filename <- paste("Spike-in_Probe_sets_Intensity_Correction_", channel.name, "_channel.pdf", sep="")
			pdf(file=filename, title=txt)
		}
		
		figure1 <- dev.cur()
		xLegend <- "Array labels"
		yLegend <- "Intensities (log2)"
		#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE),c(3,3,3), c(3,3,1), TRUE)
		layout(matrix(c(1,2,3,4,5,3), 2, 3, byrow = TRUE),c(5,5,3), c(5,5,5), TRUE)
		par(xpd = NA)
		par(oma = c(0, 0, 0, 0))
		
		# Start: Uncorrected
		m <- t(mss[ordered,])
		matplot(m, 
				type="o", col = rainbow(nrow(mss)),
				pch=1:nrow(mss),
				lty=1,
				ylab = yLegend,
				xlab = xLegend,
				main = "Start: Uncorrected",
				axes = FALSE,
				cex = 0.8)
		axis(1, at = seq(1, nrow(m), by = 1))
		axis(2, at = seq(round(min(m)), round(max(m))))
		
		# Step #1: Centered-Normed
		m <- t(nss[ordered,])
		matplot(m, 
				type="o", col = rainbow(nrow(mss)),
				pch=1:nrow(mss),
				lty=1,
				ylab = yLegend,
				xlab = xLegend,
				main = "Step #1: Centered-Normed",
				axes = FALSE,
				cex = 0.8)
		axis(1, at = seq(1, nrow(m), by = 1))
		axis(2, at = seq(-0.75, 0.75, by=0.1))
		
		# Legend
		cols <- rainbow(nrow(mss))
		for(i in 1:nrow(mss)) {
			y <- c(i, i)
			if(i==1) {
				plot(0, i, axes=FALSE, ylim=c(-1.2, nrow(mss)+1), xlim=c(0,1.5), xlab="", ylab="", col=cols[i], pch=i)
			}
			else {
				points(0,i, col=cols[i], pch=i, cex=1.2)
			}
			text(1, y, labels=(spikeNames[ordered])[i]) 
		}
		legend("top", "Spike-in probe sets", bty="n")
		
		# Step #2: Corrected Center-Normed
		m <- t(ncss[ordered,])
		matplot(m, 
				type="o", col = rainbow(nrow(mss)),
				pch=1:nrow(mss),
				lty=1,
				ylab = yLegend,
				xlab = xLegend,
				main = "Step #2: Corrected Center-Normed",
				axes = FALSE,
				cex = 0.8)
		axis(1, at = seq(1, nrow(m), by = 1))
		axis(2, at = seq(round(min(m),1), round(max(m),1), by=0.05))
		
		# Final step: Corrected
		m <- t(mcss[ordered,])
		matplot(m, 
				type="o", col = rainbow(nrow(mss)),
				pch=1:nrow(mss),
				lty=1,	
				ylab = yLegend,
				xlab = xLegend,
				main = "Final step: Corrected",
				axes = FALSE,
				cex = 0.8)
		axis(1, at = seq(1, nrow(m), by = 1))
		axis(2, at = seq(round(min(m)), round(max(m))))
		
		#legend(-3, 3.5, spikeNames[ordered], col=rainbow(nrow(m)), lty=1:nrow(m) ,ncol = 3, cex = 1.0, pch="*")
		layout(1)
		if (!dev.interactive()) dev.off()
		
		txt <- paste("Fig.2: Performance of the spike-in probe set intensity correction (", channel.name, " channel)", sep="")
		if (control.params$figures.output=="display") {
			dev.new(title=txt)
		} else {
			filename <- paste("Performance_of_the_spike-in_probe_set_intensity_correction_", channel.name, "_channel.pdf", sep="")
			pdf(file=filename, title=txt)
		}
		figure2 <- dev.cur()
		nf <- layout(matrix(c(1,2,3,0), 2, 2, byrow=TRUE), c(4,1))
		cmat <- cor(t(mss))
		cmin <- 0
		cmax <- 1
		csize <- ncol(cmat)
		
		ColorRamp <- rgb( seq(0,1,length=256),
				seq(1,0,length=256),
				seq(0,0,length=256))
		ColorLevels <- seq(cmin, cmax, length=length(ColorRamp))
		
		# Data Map
		nb.sources <- 10000
		heat.source <- heat.colors(nb.sources)
		index.map <- vector()
		for(i in 1:nb.sources)
		{
			index.map[i] = 2/pi * asin(i/nb.sources)
		}
		
		heat.panel <-  round(index.map * 128)
		res <- which(!duplicated(heat.panel))
		heat.panel <- heat.source[res]
		
		image(1:csize, 1:csize, t(cmat)[ordered,ordered], col=rev(heat.panel), xlab="Spike-in probe sets", ylab="Spike-in probsets", zlim=c(cmin,cmax), main="Pearson correlation between \nuncorrected spike-in probe sets", axes=FALSE)
		axis(1, at = seq(1, csize, by = 1),(spikeNames[ordered])[1:csize])
		axis(2, at = seq(1, csize, by = 1),(spikeNames[ordered])[1:csize])
		
		# Color Scale
		image(1, ColorLevels,
				matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
				col=rev(heat.panel),
				xlab="",ylab="",
				xaxt="n")
		
		ratio <- (apply(t(mcss),2, var) / apply(t(mss),2, var))[ordered]
		#rownames(ratio) <- spikeNames
		max.ratio <- max(ratio)
		if ( max.ratio <= 1) max.ratio <- 1 
		barplot(ratio,main="Variance ratio after/before correction\nfor spike-in probe set intensities", axes=FALSE, ylim=c(0,max.ratio), col=rainbow(nrow(mcss)), xlab="Spike-in probe sets", ylab="Variance ratio")
		axis(2, at = seq(0,max.ratio,by=0.1))
		# ylab="Corrected variance/Uncorrected variance", xlab="Spike probesets"
		layout(1)
		if (!dev.interactive()) dev.off()
	}
	
	# Set the linit of samples viewed on figure
	max.samplesviewed <- 10
	
	# Set the default loess span value
	if (control.params$loess.span == -1) control.params$loess.span <- 5.0/nrow(mss)
	
	# Compute the number of points to use for extrapolation based of spike control extrap.points needed
	spikeset.order <- order(mss[,1])
	nbpoints <- 0
	
	for(p in 1:(control.params$extrap.points)) {
		nbpoints <- nbpoints + (spikeCounts[spikeset.order])[p] 
	}
	
	points.all <- list()
	for(s in 1:nbsamples) {
		# Sort spike probe signal and correction matrices by intensity
		o <- order(msp[,s])
		x <- msp[o,s]
		y <- dmsp[o,s]
		# Get the number of integer values that need to be extrapoled to cover the all signal
		xa <- seq(from=ceiling(max(x)), to=ceiling(max(M)))
		# Check if number of points needed to be extrapoled are more than 50% of the number of spikes
		#if (length(xa) / length(spikeCounts) > 0.5) {
		#	if (verbose) {
		#		message("The number of points that have to be extrapoled is too high to perform spike-in normalization.\nUsing median normalization...")
		#	}
		#	return(mediannorm(object,channel=channel))
		#}
		
		lenx <- length(x)
		leny <- length(y)
		
		# Compute extrapolation for the queue by doing a regression with the extrap.points last sets
		nbrep <- median(spikeCounts)
		xf <- vector()
		for(i in 1:length(xa)) {
			xf <- c(xf, rep(xa[i], nbrep))
		}
		
		if (control.params$force.zero) {
			# Bug fix, if there are some zeros in the expression matrix, log2 will give -Inf
			# Get the minimum intensity which is not zero
			xb <- rep(min(M[which(M>0)]), (spikeCounts[spikeset.order])[1])
			xe <- c(xb, x, xf)
		}
		else {
			# Bug fix, if there are some zeros in the expression matrix, log2 will give -Inf
			# Get the minimum intensity which is not zero
			xe <- c(min(M[which(M>0)]), x, xf)
		}
		
		if (control.params$extrap.method == "mean") {
			meany <- mean(y[(leny-nbpoints):(leny)])
			prediction <- rep(meany, length(xf))
		}
		else {
			if (control.params$extrap.method == "linear") {
				qdata <- data.frame(x=x[(lenx-nbpoints-1):(lenx)], y=y[(leny-nbpoints-1):(leny)])
				fit <- lm(y ~ x, data=qdata )
				prediction <- predict(fit, data.frame(x=xf))
			}
			else {
				stop("extrap.method not supported")
			}
		}
		
		if (control.params$force.zero) {
			ye <- c(rep(0, (spikeCounts[spikeset.order])[1]), y, prediction)	
		}
		else {
			ye <- c(median(y[1:(spikeCounts[spikeset.order])[1]]), y, prediction)
			#ye <- c(y, prediction)
		}
		
		lo.res <- loess.smooth(xe, ye, span=control.params$loess.span)
		if (control.params$figures.show) {
			points.all[[s]] <- list()
			points.all[[s]]$x <- xe
			points.all[[s]]$y <- ye
			points.all[[s]]$lo <- lo.res
		}
		dM[,s] <- (approx(x=lo.res$x, y=lo.res$y, xout=M[,s]))$y
	}
	
	if (control.params$figures.show) {
		by.samples.figures(points.all, channel.name,figures.output=control.params$figures.output)
	}
	
	Mout <- M[,(sample.start):(sample.end)] + dM
	Mout <- 2 ^ Mout
	return(Mout)
}

by.samples.figures <- function(data, channel, figures.output=c("display","file")) {
	figures.output <- match.arg(figures.output)
	max.samples <- 20
	nbsamples <- length(data)
	nb.figures <- ceiling(nbsamples / max.samples)
	for(i in 1:nb.figures) {
		if ((i*max.samples) > nbsamples) {
			to <- nbsamples
		}
		else {
			to <- i * max.samples
		}
		from <- (i-1) * max.samples + 1		
		txt <- paste("Correction functions (channel ", channel, ") arrays ", from, " to ", to, sep="")
		if (figures.output=="display") {
			dev.new(title=txt)
		} else {
			filename <- paste("Correction_functions_channel_", channel, "_arrays_", from, "_to_", to,".pdf",sep="")
			pdf(file=filename, title=txt)
		}
		def.par <- par(no.readonly = TRUE)
		par(xpd = NA)
		par(oma = c(0, 0, 0, 0))
		nf <- layout(matrix(c(1,2), 1, 2, byrow=TRUE), c(3,1), c(3,3), TRUE)
		by.samples.curve(data, from, to)
		by.samples.legend(from, to)
		layout(1)
		if (!dev.interactive()) dev.off()
	}
}


by.samples.curve <- function(data, sample.from, sample.to) {
	nbsamples <-  sample.to - sample.from + 1
	title <- paste("Correction functions for raw intensity normalization\n(arrays ",sample.from, " to ", sample.to,")", sep="")
	
	dataloy.vect <- vector()
	for(s in 1:nbsamples) {
		index <- s + sample.from - 1
		dataloy.vect <- c(dataloy.vect, data[[index]]$lo$y)
		dataloy.vect <- c(dataloy.vect, data[[index]]$y)
	}
	min.value <- min(dataloy.vect)
	max.value <- max(dataloy.vect)
	
	for(s in 1:nbsamples) {
		index <- s + sample.from - 1
		if (s  == 1) {
			# this is a new figure on the layout
			plot(data[[index]]$x,data[[index]]$y, ylim=c(min.value, max.value), col=s, cex=0.72, pch=s, main=title, xlab="Raw intensities", ylab="Correction functions")
		} else {
			points(data[[index]]$x, data[[index]]$y, col=s, pch=s, cex=0.72)
		}
		lines(data[[index]]$lo,col=s,pch=s)
	}
}

by.samples.legend <- function(sample.from, sample.to) {
	
	nbsamples <-  sample.to - sample.from + 1
	for(s in 1:nbsamples) {
		index <- s + sample.from - 1
		if (s == 1) {
			plot(-0.5, s, axes=FALSE, ylim=c(0, nbsamples+2), xlim=c(0,1.5), xlab="", ylab="", col=s, pch=s)
		}
		else {
			points(-0.5,s, col=s, pch=s, cex=1.2)
		}
		text(1.5, s, labels=(s+sample.from-1), cex=0.6)
	}
	legend(x=-0.5, y=nbsamples+2, "Array labels", bty="n", cex=0.7)
}

meannorm <- function(object, channel=0) {
	if (channel==0) {
		m <- exprs(object)
	} else if (channel==1) {
		m <- exprs(object)[,1:length(sampleNames(object)) / 2]
	} else {
		m <- exprs(object)[,(length(sampleNames(object)) / 2 + 1):(length(sampleNames(object)))]
	}
	m <- exprs(object)
	m <- log2(m)
	cm <- colMeans(m)
	m2 <- m
	for(i in 1:nrow(m)) {
		m2[i,] <- m2[i,] - cm
	}
	m <- m - min(m)
	m <- 2 ^ m
	return(m)
}


mediannorm <- function(object, channel=0) {
	if (channel==0) {
		m <- exprs(object)
	} else if (channel==1) {
		m <- exprs(object)[,1:(length(sampleNames(object)) / 2)]
	} else {
		m <- exprs(object)[,(length(sampleNames(object)) / 2 + 1):(length(sampleNames(object)))]
	}
	m <- log2(m)
	md <- apply(m, 2, median)
	m2 <- m
	for(i in 1:nrow(m)) {
		m2[i,] <- m2[i,] - md
	}
	m <- m - min(m)
	m <- 2 ^ m
	return(m)
}

summarize.miR <- function(
		abatch,
		pmcorrect.method = 'pmonly',
		pmcorrect.param = list(),
		summary.method = 'medianpolish',
		summary.param = list(),
		summary.subset = NULL) {
	if (is.null(summary.method)) {
		stop("Please specify a summarization method, for example: 'medianpolish'")
	}
	
	if (is.null(pmcorrect.method)) {
		if (!is.from.createAB(abatch)) {
			stop("Please specify a pmcorrect method, for example: 'pmonly'")
		} else {
			pmcorrect.method = "pmonly"
		}
	}
	
	
	#pm correct
	if (pmcorrect.method != "pmonly") {
		if (!is.from.createAB(abatch)) {
			if(all(is.na(mm(abatch)))) {
				warning("No mismatch probes specified in the CDF file, using pmcorrect.method='pmonly'")
				pmcorrect.method = "pmonly"
			}
		} else {
			pmcorrect.method = "pmonly"
		}
	}
	return(computeExprSet(abatch, pmcorrect.method=pmcorrect.method,
					summary.method=summary.method, summary.param=summary.param,
					ids=summary.subset, pmcorrect.param=pmcorrect.param ))	
}


is.from.createAB <- function(object) {
	return(any("ExiMiR.channel"%in%names(notes(object))))
}

is.dual <- function(object) {
	if (!is.from.createAB(object))
		FALSE
	else {
		if (notes(object)$ExiMiR.channel == "Dual channel") {
			TRUE
		} else {
			FALSE
		}
	}
}

has.bg <- function(object) {
	if (!is.from.createAB(object))
		FALSE
	else {
		if (notes(object)$ExiMiR.bg == "Contains background") {
			TRUE
		} else {
			FALSE
		}
	}
}
