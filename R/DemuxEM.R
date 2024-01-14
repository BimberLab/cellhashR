#' @include Utils.R
#' @include Visualization.R

utils::globalVariables(
	names = c('sortOrder'),
	package = 'cellhashR',
	add = TRUE
)

# rawFeatureMatrixH5 is the raw_feature_bc_matrix.h5 file from the gene expression data
GenerateCellHashCallsDemuxEM <- function(barcodeMatrix, rawFeatureMatrixH5, methodName = 'demuxem', label = 'demuxEM', verbose= TRUE, metricsFile = NULL, doTSNE = TRUE) {
	if (verbose) {
		print('Starting demuxEM')
	}

	if (!reticulate::py_available(initialize = TRUE)) {
		stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
	}

	if (!reticulate::py_module_available('demuxEM')) {
		print('Error loading demuxEM')
		tryCatch({
			print(reticulate::import('demuxEM'))
		}, error = function(e){
			print("Error with reticulate::import('demuxEM')")
			print(conditionMessage(e))
			traceback()
		})

		stop('The demuxEM python package has not been installed!')
	}

	tryCatch({
		#Save to disk:
		inputHtoFile <- tempfile(fileext = '.csv')
		df <- as.data.frame(as.matrix(barcodeMatrix))
		df <- cbind(data.frame("HTO"=rownames(barcodeMatrix)), df)

		# demuxEM seems to expect the cellbarcodes in the HTO CSV to lack the suffix, even if the h5 data has them
		newToOldCellbarcode <- NULL
		if (sum(grepl(names(df), pattern = '-[0-9]')) > 0) {
			print('Removing cell barcode suffixes from input HTO matrix')
			newToOldCellbarcode <- data.frame(origCellbarcode = names(df))
			names(df) <- sapply(names(df), function(x){
				return(unlist(strsplit(x, split = '-[0-9]'))[1])
			})

			newToOldCellbarcode$updatedBarcode <- names(df)
		}

		write.table(df, inputHtoFile, row.names=FALSE, sep = ',', quote = FALSE)

		outPath <- tempfile()
		args <- c("-m", "demuxEM", "--random-state", GetSeed(), "--generate-diagnostic-plots", rawFeatureMatrixH5, inputHtoFile, outPath)
		pyOut <- system2(reticulate::py_exe(), args, stdout = TRUE, stderr = TRUE)
		print(pyOut)

		csvOut <- tempfile()
		zip <- paste0(outPath, '_demux.zarr.zip')
		if (!file.exists(zip)) {
			stop(paste0('Unable to find ZIP: ', zip))
		}

		tempScript <- tempfile(fileext = '.py')
		fileConn <- file(tempScript)
		writeLines(c("import pegasusio as io", paste0("data = io.read_input('", zip, "')"), paste0("data.obs.to_csv('", csvOut, "')")), fileConn)
		close(fileConn)

		pyOut2 <- system2(reticulate::py_exe(), tempScript, stdout = TRUE, stderr = TRUE)
		if (!file.exists(csvOut)) {
			stop(paste0('Unable to find CSV: ', csvOut))
		}

		df <- read.table(csvOut, header = TRUE, sep = ',', stringsAsFactors = FALSE)
		names(df) <- c('cellbarcode', 'classification.global', 'classification')
		if (!all(is.null(newToOldCellbarcode))) {
			print('Replacing original cell barcodes in demuxEM output')
			toFix <- data.frame(updatedBarcode = df$cellbarcode, sortOrder = 1:length(df$cellbarcode))
			toFix <- merge(toFix, newToOldCellbarcode, by = 'updatedBarcode', all.x = TRUE)
			toFix <- dplyr::arrange(toFix, sortOrder)

			df$cellbarcode <- toFix$origCellbarcode
		}

		df <- df[df$cellbarcode %in% colnames(barcodeMatrix),]
		df$classification[df$classification == ''] <- 'Negative'
		df$classification[is.na(df$classification) | df$classification == 'unknown'] <- 'Negative'
		df$classification[df$classification == 'singlet'] <- 'Singlet'
		df$classification[df$classification == 'doublet'] <- 'Doublet'

		df$classification[grepl(df$classification, pattern = ',')] <- 'Doublet'
		df$classification.global <- df$classification
		df$classification.global[!df$classification.global %in% c('Negative', 'Doublet')] <- 'Singlet'

		# Ensure order matches input:
		toMerge <- data.frame(cellbarcode = colnames(barcodeMatrix), sortOrder = 1:length(colnames(barcodeMatrix)))
		df <- merge(df, toMerge, by = 'cellbarcode', all.x = FALSE, all.y = FALSE)
		df <- arrange(df, sortOrder)
		df <- df[c('cellbarcode', 'classification.global', 'classification')]
		df$classification[is.na(df$classification)] <- 'Negative'
		df$classification.global[is.na(df$classification.global)] <- 'Negative'

		for (name in c('.ambient_hashtag.hist.pdf', '.real_content.hist.pdf', '.rna_demux.hist.pdf', '.background_probabilities.bar.pdf')) {
			i <- paste0(outPath, name)
			plot(magick::image_read_pdf(i))
			#unlink(i)
		}

		unlink(inputHtoFile)
		unlink(outPath)
		unlink(csvOut)
		unlink(tempScript)
		unlink(zip)

		ret <- data.frame(cellbarcode = df$cellbarcode, method = methodName, classification = df$classification, classification.global = df$classification.global, stringsAsFactors = FALSE)

		assay <- 'HTO'
		seuratObj <- Seurat::CreateSeuratObject(barcodeMatrix, assay = assay)

		toMerge <- ret$classification
		names(toMerge) <- ret$cellbarcode
		seuratObj$classification.demuxem <- toMerge[colnames(seuratObj)]
		seuratObj$classification.demuxem <- naturalsort::naturalfactor(seuratObj$classification.demuxem)

		toMerge <- ret$classification.global
		names(toMerge) <- ret$cellbarcode
		seuratObj$classification.global.demuxem <- toMerge[colnames(seuratObj)]
		seuratObj$classification.global.demuxem <- naturalsort::naturalfactor(seuratObj$classification.global.demuxem)
		SummarizeHashingCalls(seuratObj, label = label, columnSuffix = 'demuxem', assay = assay, doTSNE = F, doHeatmap = F)

		return(ret)
	}, error = function(e){
		print('Error generating demuxEM calls, aborting')
		if (!is.null(e)) {
			print(conditionMessage(e))
			traceback()
		}

		return(NULL)
	})
}
