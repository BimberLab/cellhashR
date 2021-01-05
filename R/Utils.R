# Reverse the vector x and return the value at the Nth index. If N is larger
# than the length of the vector, return the last value in the reversed vector.
#
# @param x vector of interest
# @param N index in reversed vector
#
# @return returns element at given index
#
MaxN <- function(x, N = 2){
	TopN(x, N = N, useMin = FALSE)
}

MinN <- function(x, N = 2){
	TopN(x, N = N, useMin = TRUE)
}

TopN <- function(x, N = 2, useMin = FALSE){
	len <- length(x)
	if (N > len) {
		warning('N greater than length(x).  Setting N=length(x)')
		N <- length(x)
	}
	sort(x, partial = len - N + 1, decreasing = useMin)[len - N + 1]
}

# Local maxima estimator
#
# Finding local maxima given a numeric vector
#
# @param x A continuous vector
#
# @return Returns a (named) vector showing positions of local maximas
#
# @author Tommy
# @references \url{https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima}
#
# @examples
# x <- c(1, 2, 9, 9, 2, 1, 1, 5, 5, 1)
# LocalMaxima(x = x)
#
LocalMaxima <- function(x) {
	# Use -Inf instead if x is numeric (non-integer)
	y <- diff(x = c(-.Machine$integer.max, x)) > 0L
	y <- cumsum(x = rle(x = y)$lengths)
	y <- y[seq.int(from = 1L, to = length(x = y), by = 2L)]
	if (x[[1]] == x[[2]]) {
		y <- y[-1]
	}
	return(y)
}

GetPlotColors <- function(total, palette = 'Set1') {
 	return(grDevices::colorRampPalette(RColorBrewer::brewer.pal(max(3, min(9, total)), palette))(total))
}

GetTotalPlotPages <- function(totalValues, perPage = 4) {
	return(ceiling(totalValues / perPage))
}

SimplifyHtoNames <- function(x) {
	return(gsub(x, pattern = '-[ATGC]+$', replacement = ''))
}

.InferPerplexityFromSeuratObj <- function(seuratObj, perplexity = 30) {
	return(.InferPerplexity(ncol(seuratObj), perplexity))
}

.InferPerplexity <- function(sampleNumber, perplexity = 30) {
	if (sampleNumber - 1 < 3 * perplexity) {
		print(paste0('Perplexity is too large for the number of samples: ', sampleNumber))
		perplexityNew <- floor((sampleNumber - 1) / 3)
		print(paste0('lowering from ', perplexity, ' to: ', perplexityNew))
		perplexity <- perplexityNew
	}

	return(perplexity)
}

.LogMetric <- function(metricsFile, metricName, metricValue, append = TRUE) {
	if (!is.null(metricsFile)) {
		write(x = paste0(metricName, '\t', metricValue), file = metricsFile, append = append)
	}
}

.AssignCallsToMatrix <- function(object, discrete, suffix, assay = 'HTO') {
	data <- GetAssayData(object = object, assay = assay)

	npositive <- colSums(x = discrete)
	classification.global <- npositive
	classification.global[npositive == 0] <- "Negative"
	classification.global[npositive == 1] <- "Singlet"
	classification.global[npositive > 1] <- "Doublet"
	donor.id = rownames(x = data)
	hash.max <- apply(X = data, MARGIN = 2, FUN = max)
	hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
	hash.second <- apply(X = data, MARGIN = 2, FUN = MaxN, N = 2)
	hash.maxID <- as.character(x = donor.id[sapply(
	X = 1:ncol(x = data),
		FUN = function(x) {
			return(which(x = data[, x] == hash.max[x])[1])
		}
	)])
	hash.secondID <- as.character(x = donor.id[sapply(
	X = 1:ncol(x = data),
		FUN = function(x) {
			return(which(x = data[, x] == hash.second[x])[1])
		}
	)])
	hash.margin <- hash.max - hash.second
	doublet_id <- sapply(
	X = 1:length(x = hash.maxID),
		FUN = function(x) {
			return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), collapse = "_"))
		}
	)

	classification <- classification.global
	classification[classification.global == "Negative"] <- "Negative"
	classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
	classification[classification.global == "Doublet"] <- "Doublet" #doublet_id[which(x = classification.global == "Doublet")]
	classification.metadata <- data.frame(
		hash.maxID,
		hash.secondID,
		hash.margin,
		classification,
		classification.global
	)

	colnames(x = classification.metadata) <- paste(c('maxID', 'secondID', 'margin', 'classification', 'classification.global'), suffix, sep = '.')
	object <- AddMetaData(object = object, metadata = classification.metadata)

	return(object)
}