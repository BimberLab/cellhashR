#' @include Utils.R
#' @include Normalization.R
#' @include Visualization.R



GenerateCellHashCallsThreshold <- function(barcodeMatrix, verbose = TRUE, assay = 'HTO', methodName = 'threshold', label = 'Threshold') {
	if (verbose) {
		print('Starting Threshold Calling')
	}

	tryCatch({
		seuratObj <- Seurat::CreateSeuratObject(barcodeMatrix, assay = assay)
		seuratObj[[assay]]@data <- NormalizeRelative(barcodeMatrix)

		seuratObj <- ThresholdDemux(seuratObj = seuratObj, positivity_threshold = 0.75, assay = assay)

		SummarizeHashingCalls(seuratObj, label = label, columnSuffix = 'threshold', assay = assay)

		df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = methodName, classification = seuratObj$classification.threshold, classification.global = seuratObj$classification.global.threshold, stringsAsFactors = FALSE)
		return(df)
	}, error = function(e){
		print('Error generating threshold calls, aborting')
		print(conditionMessage(e))
		traceback()

		return(NULL)
	})
}

ThresholdDemux <- function(seuratObj, positivity_threshold, assay) {
	barcodeMatrix <- GetAssayData(
		object = seuratObj,
		assay = assay,
		slot = 'data'
	)[, colnames(x = seuratObj)]

	#loop over HTOs in matrix, perform thresholding and store cells that pass the threshold
	#return a discrete matrix, with 1 equal to a call positive for that barcode
	discrete <- GetAssayData(object = seuratObj, assay = assay)
	discrete[discrete > 0] <- 0

	for (hto in rownames(barcodeMatrix)) {
		cells <- barcodeMatrix[hto, colnames(seuratObj), drop = FALSE]

		discrete[hto, colnames(seuratObj)] <- ifelse(cells > positivity_threshold, yes = 1, no = 0)
	}

	# TODO: it should be possible to call doublets too

	seuratObj <- .AssignCallsToMatrix(seuratObj, discrete, suffix = 'threshold', assay = assay)

	return(seuratObj)
}