# File: NormiR.R
# 
# Author: sylvain.gubian@pmi.com
# Aim: Set of functions for Exiqon data normalization

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#########################################################################################

NormiR <- function(	abatch,
					method=c("spikein","mean","median"),
					background.correct=FALSE,
					verbose=TRUE,
					figures.show=TRUE,
					figures.output=c("display","file"),
					out.type=c("ExpressionSet", "data.frame"),
					min.corr=0.5,
					loess.span=-1,
					extrap.points=2,
					extrap.method=c("mean","linear"),
					force.zero=FALSE,
					cover.ext=0.5,
					cover.int=1/3,
					max.log2span=1) {
	method <- match.arg(method)
	out.type <- match.arg(out.type)
	figures.output <- match.arg(figures.output)
	extrap.method <- match.arg(extrap.method)
	
	# Call background correction if asked
	if (background.correct) {
		ab <- backgroud.correction.miR(abatch)
	}
	else {
		ab <- abatch
	}
	
	# Call normalize.exi
	ab <- norm.miR(abatch=ab, method=method, figures.show=figures.show, min.corr=min.corr, loess.span=loess.span, extrap.points=extrap.points, extrap.method=extrap.method, verbose=verbose, force.zero=force.zero,cover.ext=cover.ext,cover.int=cover.int,figures.output=figures.output,max.log2span=max.log2span)
	
	# Call sumarisation
	return(summarize.miR(ab,out.type=out.type))
}

backgroud.correction.miR <- function(object) {
	ab <- bg.correct.rma(object)
	return(ab)
}

norm.miR <- function(	abatch, 
						method=c("spikein","mean","median"),
						figures.show=TRUE,
						figures.output=c("display","file"),
						min.corr=0.5,
						loess.span=-1,
						extrap.points=2,
						extrap.method=c("mean","linear"),
						force.zero=FALSE,
						cover.ext=0.5,
						cover.int=1/3,
						max.log2span=1,
						verbose=TRUE) {
	
	method <- match.arg(method)
	figures.output <- match.arg(figures.output)
	extrap.method <- match.arg(extrap.method)
	
	if (method=="spikein") {
		# Check if it is single or dual channel (should be a dual, we only handle Exiqon for now)
		if (class(abatch)[1]=="AffyBatch") {
			galenv <- getCdfInfo(abatch)
			l <- mget(ls(galenv),galenv)
			# Check if the affybatch is build from Exiqon data
			spikeList <-  which(regexpr("^spike_control", names(l)) != -1)
			if (length(spikeList)==0) {
				# No Exiqon Spike control => Affy data (One channel)
				m <- spikeinnorm(abatch, figures.show=figures.show, channel=0,min.corr=min.corr, loess.span=loess.span, extrap.points=extrap.points, extrap.method=extrap.method, verbose=verbose,force.zero=force.zero,cover.ext=cover.ext,cover.int=cover.int,figures.output=figures.output,max.log2span=max.log2span)
			}
			else {
				# Exiqon spike controls => Exiqon data (Dual channel)
				if (verbose) {
					message("Processing channel Hy3...")
				}
				m1 <- spikeinnorm(object=abatch, figures.show=figures.show, channel=1, min.corr=min.corr, loess.span=loess.span, extrap.points=extrap.points, extrap.method=extrap.method, verbose=verbose,force.zero=force.zero,cover.ext=cover.ext,cover.int=cover.int,figures.output=figures.output,max.log2span=max.log2span)
				if (verbose) {
					message("Processing channel Hy5...")
				}
				m2 <- spikeinnorm(object=abatch, figures.show=figures.show, channel=2,min.corr=min.corr, loess.span=loess.span, extrap.points=extrap.points, extrap.method=extrap.method, verbose=verbose,force.zero=force.zero,cover.ext=cover.ext,cover.int=cover.int,figures.output=figures.output,max.log2span=max.log2span)
				m <- cbind(m1, m2)	
			}
			ab <- abatch
			exprs(ab) <- m
			return(ab)
		}
		else {
			stop("Object of class:", class(abatch)[1], " not supported")
		}
		
	}
	else if (method=="mean") {
		# Mean norm method
		m <- meannorm(abatch)
		ab <- abatch
		exprs(ab) <- m
		return(ab)
		
	}
	else if (method=="median") {
		# Median norm method
		m <- mediannorm(abatch)
		ab <- abatch
		exprs(ab) <- m
		return(ab)
	}
	else {
		stop(paste("Method", method, "is not supported yet."))
	}
}


spikeinnorm <- function(	object,
					figures.show=TRUE,
					figures.output=c("display","file"),
					channel=0,
					min.corr=0.5,
					loess.span=-1,
					extrap.points=2,
					extrap.method=c("mean","linear"),
					force.zero=FALSE,
					cover.ext=0.5,
					cover.int=1/3,
					max.log2span=1,
					verbose=TRUE) {
	# Spike-in method implementation selects only spikes in the matrix
	
	extrap.method <- match.arg(extrap.method)
	figures.output <- match.arg(figures.output)
	
	galenv <- getCdfInfo(object)
	is.exiqon <- TRUE
	
	M <- log2(exprs(object))
	
	if (channel==0) {
		sample.start <- 1;
		sample.end <- length(sampleNames(object))
		channel <- "1"
	}
	else if (channel==1) {
		sample.start <- 1;
		sample.end <- length(sampleNames(object)) / 2
		channel <- "Hy3"
	}
	else if (channel==2) {
		sample.start <- length(sampleNames(object)) / 2 + 1;
		sample.end <- length(sampleNames(object))
		channel <- "Hy5"
	}
	else {
		stop("channel value is incorrect")
	}
	
	nbsamples <- (sample.end - sample.start) + 1
	if (nbsamples==1) {
		stop("Can not do the spike calibration with only one sample.")
	}
	#l <- as.list(galenv)
	l <- mget(ls(galenv),galenv)
	# Try Exiqon Spike control names...
	spikeList <-  which(regexpr("^spike_control", names(l)) != -1)
	if (length(spikeList)==0) {
		if (verbose) {
			message("No Exiqon Spike Control probes, trying with Affy Spike Control probes...")
		}
		is.exiqon <- FALSE
		spikeList <-  which(regexpr("^spike_in-control", names(l)) != -1)
		if (length(spikeList)==0) {
			stop("Can not find any Spike Control probes.")
		}
	}
	
	# Get the median spike control matrix
	mss <- matrix(0, length(spikeList), nbsamples)
	mss.max <- matrix(0, length(spikeList), nbsamples)
	mss.min <- matrix(0, length(spikeList), nbsamples)
	
	spikeCounts <- vector()
	for(i in 1:length(spikeList)) {
		mat <- l[spikeList[i]]
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
	if (mean(cmat) < min.corr) {
		if (verbose) {
			mean.display <- floor(mean(cmat) * 100) / 100
			message("Not enough common variance to guarantee good performance of spike-in normalization (", mean.display, " < ", min.corr, ").\nUsing median normalization...")
		}
		# Call mediannorm
		return(mediannorm(object))
	}
	
	# Checking specificity of spike intensity
	#	vs <- vector()
	#	for(i in 1:nbsamples) {
	#		vs[i] <- (max(M[,i+sample.start-1])-min(M[,i+sample.start-1])) / length(spikeCounts)
	#	}
	if ( median(mss.max-mss.min) > max.log2span) {
		if (verbose) {
			v.max <- floor(max.log2span * 100) / 100 
			v.observed <- floor(median(mss.max-mss.min) * 100) / 100
			message("The intensity resolution of the spike-in probesets is too coarse (",  v.observed, " > ", v.max, ") to guarantee a good performance of spike-in normalization\nUsing median normalization...")
		}
		# Call mediannorm
		return(mediannorm(object))
	}

	# Checking external coverage
	ve <- vector()
	for(i in 1:nbsamples) {
		ve[i] <- (max(mss[,i])-min(mss[,i])) / (max(M[,i+sample.start-1])-min(M[,i+sample.start-1]))
	}
	if (median(ve) < cover.ext) {
		if (verbose) {
			ve.display <- floor(median(ve) * 100) / 100 
			message("The ratio between the intensity range covered by the spike-in probes and the one covered by all probes on the array is too low (",  ve.display, " < ", cover.ext, ").\nUsing median normalization...")
		}
		# Call mediannorm
		return(mediannorm(object))
	}

	# Checking internal coverage
	vi <- vector()
	for(i in 2:nbsamples) {
		vi[i-1] <- max(mss[,i]-mss[,i-1])
	}
	if (median(vi) > cover.int) {
		if (verbose) {
			vi.display <- floor(median(vi) * 100) / 100 
			message("The size of the largest intensity interval between two consecutive spikes is too large (",  vi.display, " > ", cover.int, ").\nUsing median normalization...")
		}
		# Call mediannorm
		return(mediannorm(object))
	}
	
	# Use shorts names for spikes control
	if (is.exiqon) {
		tab <- strsplit(names(l[spikeList]), '_')
	}
	else {
		tab <- strsplit(names(l[spikeList]), '-')
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
		for(i in 1:length(spikeList)) {
			mat <- l[spikeList[i]]
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
	
	if (figures.show) {
		ordered <- order(mss[,1])
		txt <- paste("Fig.1: Correction of spike-in probeset intensities (", channel, " channel)", sep="")
		if (figures.output=="display") {
			dev.new(title=txt)
		} else {
			filename <- paste("Spike-in_Probesets_Intensity_Correction_", channel, "_channel.pdf", sep="")
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
		legend("top", "Spike-in probesets", bty="n")
		
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
		
		txt <- paste("Fig.2: Performance of the spike-in probeset intensity correction (", channel, " channel)", sep="")
		if (figures.output=="display") {
			dev.new(title=txt)
		} else {
			filename <- paste("Performance_of_the_spike-in_probeset_intensity_correction_", channel, "_channel.pdf", sep="")
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
		
		image(1:csize, 1:csize, t(cmat)[ordered,ordered], col=rev(heat.panel), xlab="Spike-in probesets", ylab="Spike-in probsets", zlim=c(cmin,cmax), main="Pearson correlation between \nuncorrected spike-in probesets", axes=FALSE)
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
		barplot(ratio,main="Variance ratio after/before correction\nfor spike-in probeset intensities", axes=FALSE, ylim=c(0,max.ratio), col=rainbow(nrow(mcss)), xlab="Spike-in probesets", ylab="Variance ratio")
		axis(2, at = seq(0,max.ratio,by=0.1))
		# ylab="Corrected variance/Uncorrected variance", xlab="Spike probesets"
		layout(1)
		if (!dev.interactive()) dev.off()
	}
	
	# Set the linit of samples viewed on figure
	max.samplesviewed <- 10
	
	# Set the default loess span value
	if (loess.span == -1) loess.span <- 5.0/nrow(mss)
	
	# Compute the number of points to use for extrapolation based of spike control extrap.points needed
	spikeset.order <- order(mss[,1])
	nbpoints <- 0
	
	for(p in 1:extrap.points) {
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
		#	return(mediannorm(object))
		#}
		
		lenx <- length(x)
		leny <- length(y)
		
		# Compute extrapolation for the queue by doing a regression with the extrap.points last sets
		nbrep <- median(spikeCounts)
		xf <- vector()
		for(i in 1:length(xa)) {
			xf <- c(xf, rep(xa[i], nbrep))
		}
		
		if (force.zero) {
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
	
		if (extrap.method == "mean") {
			meany <- mean(y[(leny-nbpoints):(leny)])
			prediction <- rep(meany, length(xf))
		}
		else {
			if (extrap.method == "linear") {
				qdata <- data.frame(x=x[(lenx-nbpoints-1):(lenx)], y=y[(leny-nbpoints-1):(leny)])
				fit <- lm(y ~ x, data=qdata )
				prediction <- predict(fit, data.frame(x=xf))
			}
			else {
				stop("extrap.method not supported")
			}
		}
		
		if (force.zero) {
			ye <- c(rep(0, (spikeCounts[spikeset.order])[1]), y, prediction)	
		}
		else {
			ye <- c(median(y[1:(spikeCounts[spikeset.order])[1]]), y, prediction)
			#ye <- c(y, prediction)
		}
			
		lo.res <- loess.smooth(xe, ye, span=loess.span)
		if (figures.show) {
			points.all[[s]] <- list()
			points.all[[s]]$x <- xe
			points.all[[s]]$y <- ye
			points.all[[s]]$lo <- lo.res
		}
		dM[,s] <- (approx(x=lo.res$x, y=lo.res$y, xout=M[,s]))$y
	}
	
	if (figures.show) {
		by.samples.figures(points.all, channel,figures.output=figures.output)
	}
	
	Mout <- M[,(sample.start):(sample.end)] + dM
	Mout <- 2 ^ Mout
	return(Mout)
}


#by.samples.figures <- function(data, channel, out.type=c("window","file")) {
#	
#	max.samples <- 10
#	nbsamples <- length(data)
#	nb.figures.total <- (nbsamples-1) %/% max.samples + 1
#	nb.pages <- nb.figures.total %/% 4 + 1
#	
#	for(p in 1:nb.pages) {
#		if (p == nb.pages) {
#			nb.figures <- nb.figures.total %% 4
#		}
#		else nb.figures <- 4
#		txt <- paste("Corrections on channel (", channel, ")", sep="")
#		dev.new(title=txt)
#		def.par <- par(no.readonly = TRUE)
#		par(xpd = NA)
#		par(oma = c(0, 0, 0, 0))
#		#message("Nb figures:", nb.figures)
#		for(i in 1:nb.figures) {
#			if (nb.figures == 1) {
#				nf <- layout(matrix(c(1,2), 1, 2, byrow=TRUE), c(3,1), c(3,3), TRUE)
#			}
#			else if (nb.figures == 2) {
#				nf <- layout(matrix(c(1,2,3,4), 1, 4, byrow=TRUE), c(3,1,3,1), c(4), TRUE)
#			}
#			else if (nb.figures == 3) {
#				nf <- layout(matrix(c(1,2,3,4,5,6,0,0), 2, 4, byrow=TRUE), c(3,1,3,1), c(3,3), TRUE)
#			}
#			else {
#				nf <- layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow=TRUE), c(3,1,3,1), c(3,3), TRUE)
#			}
#			if ((i*max.samples) > nbsamples) {
#			to <- nbsamples
#			}
#			else {
#				to <- i * max.samples
#			}
#			from <- (i-1) * max.samples + 1
#			message("Draw from ", from, " to ", to)
#			by.samples.curve(data, from, to)
#			by.samples.legend(from, to)
#		}
#		par(def.par)
#	}
#}

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
	legend(x=-0.5, y=nbsamples+2, "Spike-in probes", bty="n", cex=0.7)
}

meannorm <- function(object) {
	m <- exprs(object)
	m <- log2(m)
	m <- m - colMeans(m)
	m <- m - min(m)
	m <- 2 ^ m
	return(m)
}

mediannorm <- function(object) {
	m <- exprs(object)
	m <- log2(m)
	m <- m - apply(m, 2, median)
	m <- m - min(m)
	m <- 2 ^ m
	return(m)
}
