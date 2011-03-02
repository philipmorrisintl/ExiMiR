# Author: sylvain.gubian@pmi.com

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

###############################################################################

summarize.miR <- function(abatch,out.type=c("ExpressionSet","data.frame")) {
	
	out.type <- match.arg(out.type)
	
	# Convert galenv into list
	nb.probes = length(indexProbes(abatch, which='pm'))
	nb.samples <- length(sampleNames(abatch))
	m <- matrix(NA, nb.probes, nb.samples)
	
	galenv <- getCdfInfo(abatch)
	l <- mget(ls(galenv),galenv)
	nbsamples <- length(sampleNames(abatch))
	m <- matrix(NA, length(l), nbsamples)
	for(i in 1:length(l)) {
		for(s in 1:nbsamples) {
			mat <- l[i]
			mat <- mat[[1]]
			m[i,s] <-  median(exprs(abatch)[mat[,1],s])
		} 
	}
	m <- log2(m)
	rownames(m) <- names(l)
	colnames(m) <- sampleNames(abatch)
	
	if (out.type=="ExpressionSet") {
		return(new("ExpressionSet",
						phenoData = phenoData(abatch),
						annotation = annotation(abatch),
						protocolData = protocolData(abatch),
						experimentData = experimentData(abatch),
						exprs = m))
	}
	else if(out.type=="data.frame") {
		return(as.data.frame(m))
	}
}

