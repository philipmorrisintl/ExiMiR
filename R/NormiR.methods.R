# File: NormiR.methods.R
# 
# Author: sylvain.gubian@pmi.com

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#########################################################################################

NormiR.bgcorrect.methods <- function(object) {
	if (class(object)!="AffyBatch") {
		stop(paste("Object of class", class(object),"not supported"))
	}
	if(!is.from.createAB(object)) {
		return(bgcorrect.methods())
	} else {
		if (is.dual(object)) {
			return(c("none", "subtract", "half","minimum", "movingmin", "edwards", "normexp"))
		} else {
			return(bgcorrect.methods())
		}
	}
}

NormiR.normalize.methods <- function(object) {
	if (class(object)!="AffyBatch") {
		stop(paste("Object of class", class(object),"not supported"))
	}
	return(c("spikein","quantile","median","mean","none"))
}

NormiR.pmcorrect.methods <- function(object) {
	if (class(object)!="AffyBatch") {
		stop(paste("Object of class", class(object),"not supported"))
	}
	if (!is.from.createAB(object)) {
		if(all(is.na(mm(object)))) {
			return("pmonly")
		} else {
			return(pmcorrect.methods())
		}
	} else {
		warning("pm correction is not applicable for this data type")
		return(NULL)
	}
}

NormiR.summary.methods <- function() {
	return(generateExprSet.methods())
}

NormiR.spikein.args <- function() {
	return(c(	"min.corr",
				"loess.span",
				"extrap.points",
				"extrap.method",
				"force.zero",
				"cover.ext",
				"cover.int",
				"max.log2span",
				"figures.show",
				"figures.output"))
}