#' @include Utils.R
#' @include Visualization.R
#' @importFrom dplyr %>%
#'
GenerateCellHashCallsGMMDemux <- function(barcodeMatrix, methodName = 'gmm_demux', label = 'GMM Demux', verbose= TRUE, metricsFile = NULL, doTSNE = TRUE, doHeatmap = TRUE) {
	if (verbose) {
		print('Starting GMM-Demux')
	}

	if (!reticulate::py_available(initialize = TRUE)) {
		stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
	}

	if (!reticulate::py_module_available('GMM_Demux')) {
		stop('GMM_Demux has not been installed!')
	}

	tryCatch({
		#Save to disk:
		inputFile <- tempfile(fileext = '.csv')
		write.table(t(barcodeMatrix), file = inputFile, sep = ',')

		reportPath <- tempfile()
		outPath <- tempfile()
		pyOut <- system2(reticulate::py_exe(), c("-m", "GMM_Demux.GMM_Demux", inputFile, paste0(rownames(barcodeMatrix), collapse=','), '-c', '-f', reportPath, '-o', outPath), stdout = TRUE, stderr = TRUE)
		print(pyOut)

		clusterNames <- read.table(paste0(reportPath, '/GMM_full.config'), header = FALSE, sep = ',')
		names(clusterNames) <- c('Cluster_id', 'classification')
		clusterNames$classification <- gsub(x = clusterNames$classification, pattern = '^ ', replacement = '')

		df <- read.table(paste0(reportPath, '/GMM_full.csv'), header = TRUE, sep = ',', row.names = 1, stringsAsFactors = FALSE)
		df$cellbarcode <- rownames(df)
		df <- merge(df, clusterNames, by = 'Cluster_id')
		df$classification[df$classification == 'negative'] <- 'Negative'
		df$classification[!df$classification %in% c('Negative', rownames(barcodeMatrix))] <- 'Doublet'
		df$classification.global <- df$classification
		df$classification.global[!df$classification.global %in% c('Negative', 'Doublet')] <- 'Singlet'

		unlink(inputFile)
		unlink(reportPath, recursive = TRUE)
		unlink(outPath, recursive = TRUE)

		ret <- data.frame(cellbarcode = df$cellbarcode, method = methodName, classification = df$classification, classification.global = df$classification.global, stringsAsFactors = FALSE)
		assay <- 'HTO'
		seuratObj <- suppressWarnings(Seurat::CreateSeuratObject(barcodeMatrix, assay = assay))

		toMerge <- ret$classification
		names(toMerge) <- ret$cellbarcode
		seuratObj$classification.gmm_demux <- toMerge[colnames(seuratObj)]
		seuratObj$classification.gmm_demux <- naturalsort::naturalfactor(seuratObj$classification.gmm_demux)

		toMerge <- ret$classification.global
		names(toMerge) <- ret$cellbarcode
		seuratObj$classification.global.gmm_demux <- toMerge[colnames(seuratObj)]
		seuratObj$classification.global.gmm_demux <- naturalsort::naturalfactor(seuratObj$classification.global.gmm_demux)
		SummarizeHashingCalls(seuratObj, label = label, columnSuffix = 'gmm_demux', assay = assay, doTSNE = doTSNE, doHeatmap = doHeatmap)

		return(ret)
	}, error = function(e){
		print('Error generating GMMDemux calls, aborting')
		if (!is.null(e)) {
			print(conditionMessage(e))
			traceback()
		}

		return(NULL)
	})
}
