pkg.env <- new.env(parent=emptyenv())

pkg.env$RANDOM_SEED <- 1234
set.seed(pkg.env$RANDOM_SEED)


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

.LogMetric <- function(metricsFile, metricName, metricValue, category = 'CellHashing', append = TRUE) {
	if (!is.null(metricsFile)) {
		write(x = paste0(category, '\t', metricName, '\t', metricValue), file = metricsFile, append = append)
	}
}

.AssignCallsToMatrix <- function(object, discrete, suffix, assay = 'HTO') {
	npositive <- Matrix::colSums(x = discrete)
	classification.global <- npositive
	classification.global[npositive == 0] <- "Negative"
	classification.global[npositive == 1] <- "Singlet"
	classification.global[npositive > 1] <- "Doublet"

	hash.called <- rownames(object)[apply(X = discrete, MARGIN = 2, FUN = which.max)]

	classification <- classification.global
	classification[classification.global == "Negative"] <- "Negative"
	classification[classification.global == "Singlet"] <- hash.called[which(x = classification.global == "Singlet")]
	classification[classification.global == "Doublet"] <- "Doublet"
	classification.metadata <- data.frame(
		classification,
		classification.global
	)

	# NOTE: leave these as strings to allow easier rbind later:
	#classification.metadata$classification <- naturalsort::naturalfactor(as.character(classification.metadata$classification))
	#classification.metadata$classification.global <- naturalsort::naturalfactor(as.character(classification.metadata$classification.global))

	colnames(x = classification.metadata) <- paste(c('classification', 'classification.global'), suffix, sep = '.')
	object <- AddMetaData(object = object, metadata = classification.metadata)

	return(object)
}

.EnsureNonSparse <- function(mat) {
	if (class(mat)[1] == 'dgCMatrix') {
		mat <- as.matrix(mat)
	}

	return(mat)
}

#' @title Set random seed
#'
#' @description Sets the seed used for R‘s random number generator, which should be used in all internal functions
#' @param seed The random seed
#' @export
SetSeed <- function(seed) {
	pkg.env$RANDOM_SEED <- seed
	set.seed(pkg.env$RANDOM_SEED)
}

#' @title Get random seed
#'
#' @description Sets a random seed, which should be used in all internal functions
#' @export
GetSeed <- function() {
	return(pkg.env$RANDOM_SEED)
}

#' @importFrom magrittr %>%
.RestoreUnderscoreToHtoNames <- function(df, originalBarcodeNames) {
	# Note: Seurat seems to replace underscores with hyphen, so check/replace these:
	for (hto in originalBarcodeNames) {
		changedVersion <- gsub(hto, pattern = '_', replacement = '-')
		if (changedVersion == hto) {
			next
		}

		if (changedVersion %in% unique(df$classification)) {
			isFactor <- is.factor(df$classification)
			df$classification <- as.character(df$classification)
			df$classification[df$classification == changedVersion] <- hto
			if (isFactor) {
				df$classification <- naturalsort::naturalfactor(df$classification)
			}
		}
	}

	return(df)
}

#' @title Estimate Multiplet Rate for 10x Genomics Data
#'
#' @description Estimates multiplet rate by using barcode matrix dimensions and formulae from the Saijita Lab's calculator (https://satijalab.org/costpercell/)
#' @param numCellsRecovered The number of cells recovered (number of columns in the count matrix)
#' @param num10xRuns The number of lanes the cells are loaded into (expected to be 1)
#' @param chemistry The version of 10x chemistry reagent kit used
#' @export
EstimateMultipletRate <- function(numCellsRecovered, num10xRuns = 1, chemistry = "10xV3"){
	if (chemistry == "10xV2"){
		# value extrapolated from inverse cell recovery rate in user guide table (page 6)
		# https://assets.ctfassets.net/an68im79xiti/RT8DYoZzhDJRBMrJCmVxl/6a0ed8015d89bf9602128a4c9f8962c8/CG00052_SingleCell3_ReagentKitv2UserGuide_RevF.pdf
		inverse_recovery_rate <- 1.74
	} else if(chemistry == "10xV3"){
		# value extrapolated from inverse cell recovery rate in user guide table (page 22)
		# https://downloads.ctfassets.net/an68im79xiti/1Y7U6QKKFz7jfQbkErDUoB/fde488a22c6d0a0fcec4597bc1a0338e/CG000416_Chromium_NextGEM_SingleCell3-_HT_v3.1_GeneExp_RevB.pdf
		inverse_recovery_rate <- 1.61
	}

	numCellsLoaded <- inverse_recovery_rate * numCellsRecovered
	m <- 4.597701e-06
	multipletRate <- m*numCellsLoaded/num10xRuns

	#these seem unncessary for now, but just in case
	#numBatches <- nrow(barcodeMatrix)
	#multipletRateIdent <- multipletRate*(numBatches-1)/numBatches
	#multipletRateNonIdent <- multipletRate-multipletRateIdent

	return(multipletRate)
}

.LogProgress <- function(msg) {
	if (Sys.getenv('CELLHASHR_DEBUG', unset = 0) == 1) {
		message(msg)
	}
}

SetAssayData4Or5 <- function(seuratObj, theLayer, new.data, ...) {
	if (!is.matrix(new.data)) {
		warning('Assay data is not a matrix, converting to a sparse matrix!')
		new.data <- Seurat::as.sparse(as.matrix(new.data))
	}

	if (!seuratObj@version < '5.0.0') {
		return(Seurat::SetAssayData(seuratObj, slot = theLayer, new.data = new.data, ...))
	} else {
		return(Seurat::SetAssayData(seuratObj, layer = theLayer, new.data = new.data, ...))
	}
}

GetAssayData4Or5 <- function(seuratObj, theLayer, ...) {
	if (!seuratObj@version < '5.0.0') {
		return(Seurat::GetAssayData(seuratObj, slot = theLayer, ...))
	} else {
		return(Seurat::GetAssayData(seuratObj, layer = theLayer, ...))
	}
}