#' @include Utils.R
#' @include Visualization.R

GenerateCellHashCallsGMMDemux <- function(barcodeMatrix, methodName = 'gmm_demux', verbose= TRUE, metricsFile = NULL) {
	if (verbose) {
		print('Starting GMM-Demux')
	}

	invisible(reticulate::py_config())

	if (!reticulate::py_available()) {
		stop('python/reticulate are not installed!')
	}

	if (!reticulate::py_module_available('GMM_Demux')) {
		stop('GMM_Demux has not been installed!')
	}

	tryCatch({
		#Save to disk:
		inputFile <- tempfile(fileext = '.csv')
		write.table(t(barcodeMatrix), file = inputFile, sep = ',')

		reportPath <- tempfile(tmpdir = '/Users/bimber/Downloads/cellhashR')
		pyOut = system2(reticulate::py_exe(), c("-m", "GMM_Demux.GMM_Demux", inputFile, paste0(rownames(barcodeMatrix), collapse=','), '-c', '-f', reportPath), stdout = TRUE, stderr = TRUE)
		print(pyOut)

		clusterNames <- read.table(paste0(reportPath, '/GMM_full.config'), header = FALSE, sep = ',')
		names(clusterNames) <- c('Cluster_id', 'classification')
		clusterNames$classification <- gsub(x = clusterNames$classification, pattern = '^ ', replacement = '')

		df <- read.table(paste0(reportPath, '/GMM_full.csv'), header = TRUE, sep = ',', row.names = 1)
		df$cellbarcode <- rownames(df)
		df <- merge(df, clusterNames, by = 'Cluster_id')
		df$classification[df$classification == 'negative'] <- 'Negative'
		df$classification[!df$classification %in% c('Negative', rownames(barcodeMatrix))] <- 'Doublet'
		df$classification.global <- df$classification
		df$classification.global[!df$classification.global %in% c('Negative', 'Doublet')] <- 'Singlet'

		unlink(inputFile)
		unlink(reportPath, recursive = TRUE)

		return(data.frame(cellbarcode = df$cellbarcode, method = methodName, classification = df$classification, classification.global = df$classification.global, stringsAsFactors = FALSE))
	}, error = function(e){
		print('Error generating GMMDemux calls, aborting')
		if (!is.null(e)) {
			print(conditionMessage(e))
			traceback()
		}

		return(NULL)
	})
}
