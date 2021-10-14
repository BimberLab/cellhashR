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
		stop('The demuxEM python package has not been installed!')
	}

	tryCatch({
		#Save to disk:
		inputHtoFile <- tempfile(fileext = '.csv')
		df <- data.frame("HTO"=rownames(barcodeMatrix), barcodeMatrix)

		# demuxEM seems to expect the cellbarcodes in the HTO CSV to lack the suffix, even if the h5 data has them
		newToOldCellbarcode <- NULL
		if (sum(grepl(names(df), pattern = '-[0-9]')) > 0) {
			print('Adding removing cell barcode suffixes from input HTO matrix')
			newToOldCellbarcode <- data.frame(origCellbarcode = names(df))
			names(df) <- unique(sapply(names(df), function(x){
				return(unlist(strsplit(x, split = '-[0-9]'))[1])
			}))

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

		pyOut2 <- system2(reticulate::py_exe(), c("-c", paste0('import pegasusio as io; data = io.read_input("', zip, '"); data.obs.to_csv("', csvOut, '")')), stdout = TRUE, stderr = TRUE)
		print(pyOut2)

		if (!file.exists(csvOut)) {
			stop(paste0('Unable to find CSV: ', csvOut))
		}

		df <- read.table(csvOut, header = TRUE, sep = ',')
		names(df) <- c('cellbarcode', 'classification.global', 'classification')
		if (!all(is.null(newToOldCellbarcode))) {
			toFix <- data.frame(updatedBarcode = df$cellbarcode, sortOrder = 1:length(df$cellbarcode))
			toFix <- merge(toFix, newToOldCellbarcode, by = 'updatedBarcode', all.x = TRUE)
			toFix <- dplyr::arrange(toFix, sortOrder)

			df$cellbarcode <- toFix$origCellbarcode
		}

		df$classification[df$classification == ''] <- 'Negative'
		df$classification[is.na(df$classification)] <- 'Negative'
		df$classification[df$classification == 'singlet'] <- 'Singlet'
		df$classification[df$classification == 'doublet'] <- 'Doublet'

		df$classification[grepl(df$classification, pattern = ',')] <- 'Doublet'
		df$classification.global <- df$classification
		df$classification.global[!df$classification.global %in% c('Negative', 'Doublet')] <- 'Singlet'

		# Ensure order matches input:
		toMerge <- data.frame(cellbarcode = colnames(barcodeMatrix), sortOrder = 1:length(barcodeMatrix))
		df <- merge(df, toMerge, by = 'cellbarcode', all.x = FALSE, all.y = TRUE)
		df <- arrange(df, sortOrder)
		df <- df[c('cellbarcode', 'classification.global', 'classification')]
		df$classification[is.na(df$classification)] <- 'Negative'
		df$classification.global[is.na(df$classification.global)] <- 'Negative'

		for (name in c('.ambient_hashtag.hist.pdf', '.real_content.hist.pdf', '.rna_demux.hist.pdf', '.background_probabilities.bar.pdf')) {
			i <- paste0(outPath, name)
			print(magick::image_read_pdf(i))
			unlink(i)
		}

		unlink(inputHtoFile)
		unlink(outPath)
		unlink(csvOut)

		ret <- data.frame(cellbarcode = df$cellbarcode, method = methodName, classification = df$classification, classification.global = df$classification.global, stringsAsFactors = FALSE)
		assay <- 'HTO'
		seuratObj <- suppressWarnings(Seurat::CreateSeuratObject(barcodeMatrix, assay = assay))

		toMerge <- ret$classification
		names(toMerge) <- ret$cellbarcode
		seuratObj$classification.demuxEM <- toMerge[colnames(seuratObj)]
		seuratObj$classification.demuxEM <- naturalsort::naturalfactor(seuratObj$classification.demuxEM)

		toMerge <- ret$classification.global
		names(toMerge) <- ret$cellbarcode
		seuratObj$classification.global.demuxEM <- toMerge[colnames(seuratObj)]
		seuratObj$classification.global.demuxEM <- naturalsort::naturalfactor(seuratObj$classification.global.demuxEM)
		SummarizeHashingCalls(seuratObj, label = label, columnSuffix = 'demuxEM', assay = assay, doTSNE = F, doHeatmap = F)

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
