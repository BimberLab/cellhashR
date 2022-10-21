#' @include Utils.R

#' @title Process CiteSeq Count Matrix
#'
#' @description The primary entrypoint for parsing and QC of the cell hashing count matrix.
#' @param rawCountData The input barcode file or umi_count folder
#' @param minCountPerCell Cells (columns) will be dropped if their total count is less than this value.
#' @param barcodeWhitelist A vector of barcode names to retain.
#' @param barcodeBlacklist A vector of barcodes names to discard.
#' @param cellbarcodeWhitelist If provided, the raw count matrix will be subset to include only these cells. This allows one to use the cellranger unfiltered matrix as an input, but filter based on target cells, such as those with GEX data. This can either be a character vector of barcodes, or a file with one cell barcode per line.
#' @param doPlot If true, QC plots will be generated
#' @param simplifyBarcodeNames If true, the sequence tag portion will be removed from the barcode names (i.e. HTO-1-ATGTGTGA -> HTO-1)
#' @param saveOriginalCellBarcodeFile An optional file path, where the set of original cell barcodes, prior to filtering, will be written. The primary use-case is if the count matrix was generated using a cell whitelist (like cells with passing gene expression). Preserving this list allows downstream reporting.
#' @param metricsFile If provided, summary metrics will be written to this file.
#' @param minCellsToContinue Demultiplexing generally requires a minimal amount of cells. If the matrix contains fewer than this many cells, it will abort.
#' @import ggplot2
#' @import patchwork
#' @import utils
#' @return The updated count matrix
#' @export
ProcessCountMatrix <- function(rawCountData=NA, minCountPerCell = 5, barcodeWhitelist = NULL, barcodeBlacklist = c('no_match', 'total_reads', 'unmapped'), cellbarcodeWhitelist = NULL, doPlot = TRUE, simplifyBarcodeNames = TRUE, saveOriginalCellBarcodeFile = NULL, metricsFile = NULL, minCellsToContinue = 25) {
	.LogProgress('Processing raw matrix')
	barcodeData <- .LoadCountMatrix(rawCountData = rawCountData, barcodeBlacklist = barcodeBlacklist, simplifyBarcodeNames = simplifyBarcodeNames)

	print(paste0('Initial cell barcodes in hashing data: ', ncol(barcodeData)))

	if (is.null(minCellsToContinue) || is.na(minCellsToContinue)) {
		minCellsToContinue <- 0
	}

	inputBarcodes <- ncol(barcodeData)
	if (!is.null(cellbarcodeWhitelist)) {
		if (is.character(cellbarcodeWhitelist) && length(cellbarcodeWhitelist) == 1) {
			if (file.exists(cellbarcodeWhitelist)) {
				cellbarcodeWhitelist <- read.table(cellbarcodeWhitelist, header = FALSE)[,1]
				print(paste0('cells in cellbarcodeWhitelist file: ', length(cellbarcodeWhitelist)))
			} else {
				warning(paste0('cellbarcodeWhitelist appears to be a filename, but it doesnt exist: ', cellbarcodeWhitelist))
			}
		}

		cellbarcodeWhitelistToUse <- intersect(cellbarcodeWhitelist, colnames(barcodeData))
		print(paste0('intersect between cellbarcodeWhitelist and barcode matrix: ', length(cellbarcodeWhitelistToUse)))
		if (length(cellbarcodeWhitelistToUse) < minCellsToContinue) {
			stop(paste0('cell barcode whitelist length of ', length(cellbarcodeWhitelistToUse), ' is less than minCellsToContinue'))
		}

		# As a sanity check, find the top cells by counts and report intersect:
		topByCounts <- sort(colSums(barcodeData), decreasing = T)
		cellsShared <- length(intersect(names(topByCounts)[1:length(cellbarcodeWhitelist)], cellbarcodeWhitelist))

		barcodeData <- barcodeData[ , cellbarcodeWhitelistToUse, drop = FALSE]
		print(paste0('Subsetting based on whitelist. Cells in whitelist: ', length(cellbarcodeWhitelist), ', cells in matrix after subset: ', ncol(barcodeData)))
		print(paste0('Total cells shared between whitelist and top droplets by count: ', cellsShared, ' (', round(100 * cellsShared / length(cellbarcodeWhitelist), 1),'%)'))
	}
	.LogMetric(metricsFile, 'InitialCellBarcodes', ncol(barcodeData))

	if (ncol(barcodeData) < minCellsToContinue) {
		stop(paste0('Too few cells remain after filtering by allowable cell barcodes, aborting. Cell count: ', ncol(barcodeData), ', original cells: ', inputBarcodes))
	}

	if (!is.null(saveOriginalCellBarcodeFile)) {
		toWrite <- data.frame(cellbarcode = colnames(barcodeData))
		write.table(toWrite, file = saveOriginalCellBarcodeFile, quote = FALSE, row.names = FALSE, col.names = FALSE)
	}

	# Print QC of counts by cell and barcode:
	if (doPlot) {
		PrintRowQc(barcodeData)
	}

	if (!is.null(barcodeWhitelist)) {
		print(paste0('Limiting to barcodes: ', paste0(barcodeWhitelist, collapse = ',')))
		sel <- rownames(barcodeData) %in% barcodeWhitelist
		if (sum(sel) == 0) {
			sel <- SimplifyHtoNames(rownames(barcodeData)) %in% barcodeWhitelist
		}
		
		barcodeData <- barcodeData[sel,]
		if (ncol(barcodeData) < minCellsToContinue) {
			stop(paste0('Too few cells remain after limiting to selected HTOs, aborting. Cell count: ', ncol(barcodeData)))
		}
	}

	if (!is.null(minCountPerCell)) {
		barcodeData <- DoCellFiltering(barcodeData, minCountPerCell = minCountPerCell)
		if (ncol(barcodeData) < minCellsToContinue) {
			stop(paste0('Too few cells remain after dropping low-count cells, aborting. Cell count: ', ncol(barcodeData)))
		}
	}
	.LogMetric(metricsFile, 'PassingCellBarcodes', ncol(barcodeData))

	if (doPlot) {
		PrintColumnQc(barcodeData)
	}

	return(barcodeData)
}

.LoadCountMatrix <- function(rawCountData = NA, barcodeBlacklist = c('no_match', 'total_reads', 'unmapped'), simplifyBarcodeNames = TRUE) {
	if (is.na(rawCountData)){
		stop("Need to provide a directory or file for rawCountData")
	}

	if (all(is.na(barcodeBlacklist)) || all(is.null(barcodeBlacklist))) {
		barcodeBlacklist <- character()
	}

	if (dir.exists(rawCountData)) {
		barcodeData <- Seurat::Read10X(rawCountData, gene.column=1, strip.suffix = TRUE)
		barcodeData <- barcodeData[which(!(rownames(barcodeData) %in% barcodeBlacklist)), , drop = F]
		barcodeData <- as.matrix(barcodeData)
	} else if (!file.exists(rawCountData)){
		stop(paste0("File does not exist: ", rawCountData))
	} else {
		# older CITE-seq-Count versions created a CSV file, so support this:
		barcodeData <- utils::read.table(rawCountData, sep = ',', header = T, row.names = 1)
		barcodeData <- barcodeData[which(!(rownames(barcodeData) %in% barcodeBlacklist)),]
		barcodeData <- as.matrix(barcodeData)
	}

	if (simplifyBarcodeNames) {
		rownames(barcodeData) <- SimplifyHtoNames(rownames(barcodeData))
	}

	return(barcodeData)
}

utils::globalVariables(
	names = c('Barcode', 'Value', 'CellBarcode', 'Freq', 'Str', 'mean_nonzero', 'barcodeMatrix'),
	package = 'cellhashR',
	add = TRUE
)

PrintRowQc <- function(barcodeMatrix) {
	df <- GenerateByRowSummary(barcodeMatrix)
	df$Barcode <- naturalsort::naturalfactor(df$Barcode)
	df$Barcode <- forcats::fct_reorder(df$Barcode, df$rowSums, .desc = TRUE)

	# Total Counts
	P1 <- ggplot(df, aes(x = Barcode, y = rowSums)) +
		geom_bar(stat = 'identity') +
		ggtitle('Raw Counts/Barcode') +
		ylab('Row Sums') +
		egg::theme_presentation(base_size = 14) +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	P2 <- ggplot(df, aes(x = Barcode, y = rowSums)) +
		geom_bar(stat = 'identity') +
		egg::theme_presentation(base_size = 14) +
		ylab('Row Sums (log1p)') +
		scale_y_continuous(trans = scales::log1p_trans()) +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	print(P1 | P2)

	# Mean Counts
	P3 <- ggplot(df, aes(x = Barcode, y = mean)) +
		geom_bar(stat = 'identity') +
		ggtitle('Mean Counts/Cell') +
		ylab('Mean Counts/Cell') +
		egg::theme_presentation(base_size = 14) +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	P4 <- ggplot(df, aes(x = Barcode, y = mean)) +
		geom_bar(stat = 'identity') +
		egg::theme_presentation(base_size = 14) +
		ylab('Mean Counts/Cell (log1p)') +
		scale_y_continuous(trans = scales::log1p_trans()) +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	print(P3 | P4)
}

GenerateByRowSummary <- function(barcodeMatrix) {
	if (nrow(barcodeMatrix) == 0) {
		print('No rows in data, cannot generate summary')
		return(data.frame(Barcode = character()))
	}

	barcodeMatrix <- as.matrix(barcodeMatrix)
	if (sum(is.na(barcodeMatrix)) > 0) {
		print(paste0('The barcodeMatrix has NAs, total: ', sum(is.na(barcodeMatrix))))
	}

	df <- data.frame(Barcode = naturalsort::naturalfactor(SimplifyHtoNames(rownames(barcodeMatrix))), BarcodeFull = naturalsort::naturalfactor(rownames(barcodeMatrix)), min = apply(barcodeMatrix, 1, min), max = apply(barcodeMatrix, 1, max), mean = apply(barcodeMatrix, 1, mean), logmean = log(apply(barcodeMatrix, 1, mean) + 1), nonzero = apply(barcodeMatrix, 1, function(x){
		sum(x > 0)
	}), mean_nonzero = (rowSums(barcodeMatrix) / rowSums(!!barcodeMatrix)), total_gt1 = apply(barcodeMatrix, 1, function(x){
		sum(x > 1)
	}), rowSums = rowSums(barcodeMatrix))

	df$mean_nonzero[is.na(df$mean_nonzero)] <- 0

	df <- df[order(df$Barcode),]

	if (sum(is.na(df)) > 0) {
		print(paste0('The barcodeMatrix had NAs after transform, total: ', sum(is.na(df))))
		for (colName in colnames(df)) {
			if (sum(is.na(df[[colName]])) > 0) {
				print(paste0('Column has NAs: ', colName))
			}
		}
	}

	return(df)
}

DoRowFiltering <- function(barcodeData, minCountPerRow = 1, doLog = TRUE) {
	toDropRows <- rowSums(barcodeMatrix) < minCountPerRow
	if (sum(toDropRows) > 0) {
		if (doLog) {
			print(paste0('Rows dropped due to low counts (<', minCountPerRow, '): ', paste0(rownames(barcodeMatrix)[toDropRows])))
		}

		barcodeData <- barcodeData[!toDropRows, , drop = FALSE]

		if (doLog) {
			print(paste0('Rows after row filter: ', nrow(barcodeData)))
		}
	}

	return(barcodeData)
}

DoCellFiltering <- function(barcodeData, minCountPerCell = 5){
	if (!is.null(minCountPerCell)) {
		toDrop <- sum(colSums(barcodeData) < minCountPerCell)
		if (toDrop > 0){
			print(paste0('Cells dropped due to low total counts per cell (<', minCountPerCell, '): ', toDrop))
			barcodeData <- barcodeData[,which(colSums(barcodeData) >= minCountPerCell), drop = FALSE]
			print(paste0('After cell filter: ', ncol(barcodeData)))
		}
	}

	return(barcodeData)
}

utils::globalVariables(
	names = c('Barcode1', 'pred'),
	package = 'cellhashR',
	add = TRUE
)

PrintColumnQc <- function(barcodeMatrix) {
	df <- data.frame(t(barcodeMatrix))
	colnames(df) <- SimplifyHtoNames(rownames(barcodeMatrix))
	df <- tidyr::gather(df, Barcode, Count)

	values <- sort(unique(df$Count), decreasing = FALSE)
	
	toPlot <- NULL
	
	for (barcode in unique(df$Barcode)) {
		bd <- df[df$Barcode == barcode,]
		dat <- sapply(values, function(x) {
			return(sum(bd$Count > x))
		})
		
		toAdd <- data.frame(Barcode = barcode, Value = values, Count = dat)
		if (is.null(toPlot)) {
			toPlot <- toAdd
		} else {
			toPlot <- rbind(toPlot, toAdd)
		}
	}

	# Counts/cell:
	df2 <- df
	df2$Count <- log10(df2$Count + 1)
	P1 <- ggplot(df2, aes(y = Count, x = Barcode)) +
		geom_violin(position="dodge", alpha=0.5) +
		xlab("") +
		ylab("Log Raw Counts") +
		ggplot2::ggtitle("Raw HTO Counts/Cell") +
		egg::theme_presentation(base_size = 14) +
		theme(axis.text.x = element_text(angle = 90))

	print(P1)

	#normalize columns, print top barcode fraction. note, cells with all zeros will make NAs
	normalizedBarcodes <- sweep(barcodeMatrix, 2, colSums(barcodeMatrix),`/`)
	normalizedBarcodes[is.na(normalizedBarcodes)] <- 0
	topValue <- apply(normalizedBarcodes,2,function(x){
		max(x)
	})

	df <- data.frame(Barcode1 = topValue, TotalPerCell = colSums(barcodeMatrix))
	P3 <- ggplot(df, aes(x = Barcode1)) +
		geom_histogram(binwidth = 0.05) +
		egg::theme_presentation() +
		xlab('Fraction') +
		ylab('# Cells') + ggtitle('Top Barcode Fraction Per Cell') +
		expand_limits(x = c(0, 1))

	P3 <- P3 + plot_annotation(caption = paste0('Total cells where top barcode is >0.75 of counts: ', sum(topValue > 0.75), ' of ', length(topValue))) & theme(plot.caption = element_text(size = 14))

	print(P3)

	P4 <- ggplot(df, aes(x = Barcode1, y = TotalPerCell)) +
		geom_point() +
		egg::theme_presentation() +
		xlab('Top Barcode Fraction') +
		ylab('Total Counts/Cell') +
		ggtitle('Top Barcode Fraction Per Cell')

	P4L <- ggplot(df, aes(x = Barcode1, y = TotalPerCell)) +
		geom_point() +
		egg::theme_presentation() +
		xlab('Top Barcode Fraction') +
		ylab('Total Counts/Cell (log1p)') +
		scale_y_continuous(trans = scales::log1p_trans()) +
		ggtitle('Top Barcode Fraction Per Cell')

	print(P4 + P4L)

	# Top/Second:
	snr <- SNR(t(barcodeMatrix))
	snr$Barcode <- naturalsort::naturalfactor(snr$Barcode)
	snr$Highest <- log10(snr$Highest + 1)
	snr$Second <- log10(snr$Second + 1)
	P1 <- ggplot2::ggplot(snr, aes(x=Highest, y=Second, color=Barcode)) +
		geom_point(cex = 0.5, alpha = 0.5) +
		ggtitle("Ratio of Top Two Barcodes") +
		egg::theme_presentation(base_size = 10) +
		xlab('Highest (log10p)') +
		ylab('Second (log10p)') +
		ylim(0, NA) +
		xlim(0, NA)

	print(P1)

	#MA-plot
	if (nrow(barcodeMatrix) == 2){
		tryCatch({
			df <- data.frame(t(barcodeMatrix))
			M <- log2(df[,1]) - log2(df[,2])
			A <- (log2(df[,1]) + log2(df[,2]))/2

			df$M <- M
			df$A <- A
			df <- df[is.finite(df$M) & is.finite(df$A),]
			o <- order(A)
			a <- A[o]
			m <- M[o]
			ind <- round(seq(1, length(a), len = 5000))
			a <- a[ind]
			m <- m[ind]
			fit <- stats::loess(m ~ a)
			bias <- stats::predict(fit, newdata = data.frame(a = A))
			df$nM <- M - bias
			newdat <- data.frame(a, pred = fit$fitted)

			print(ggplot(df, aes(x=A, y=M)) +
				geom_point() +
				geom_line(data = newdat, aes(x=a, y = pred), size = 1, col=2) +
				ggtitle("MA-plot with Loess Fit of Bias") +
				egg::theme_presentation()
			)
		}, error = function(e){
			print('Error generating MA plot, skipping')
			print(conditionMessage(e))
			traceback()
		})
	}
}

utils::globalVariables(
	names = c('CountsPerCell', 'Saturation'),
	package = 'cellhashR',
	add = TRUE
)

#' @title Plot Library Saturation
#'
#' @description Create a plot of the library saturation per cell
#' @param citeseqCountDir The root of the Cite-seq-Count output folder, which should contain umi_count and read_count folders.
#' @param metricsFile If provided, summary metrics will be written to this file.
#' @return The overall saturation for this library
#' @export
PlotLibrarySaturation <- function(citeseqCountDir, metricsFile = NULL) {
	countData <- Seurat::Read10X(paste0(citeseqCountDir, '/read_count'), gene.column=1, strip.suffix = TRUE)
	countData <- countData[rownames(countData) != 'unmapped',]
	countData <- colSums(countData)

	umiData <- Seurat::Read10X(paste0(citeseqCountDir, '/umi_count'), gene.column=1, strip.suffix = TRUE)
	umiData <- umiData[rownames(umiData) != 'unmapped',]
	umiData <- colSums(umiData)

	saturation <- 1 - (umiData / countData)

	df <- data.frame(CountsPerCell = countData, Saturation = saturation)
	df <- df %>% arrange(CountsPerCell)

	overall <- 1 - round((sum(umiData) / sum(countData)), 2)
	.LogMetric(metricsFile, 'HashingLibrarySaturation', overall)

	print(ggplot(df, aes(x = CountsPerCell, y = Saturation)) +
		labs(x = 'Counts/Cell', y = '% Saturation') + egg::theme_presentation(base_size = 18) +
		geom_point() +
		annotate("text", x = max(df$CountsPerCell), y = min(df$Saturation), hjust = 1, vjust = -1, label = paste0(
			'Total Counts: ', format(sum(countData), big.mark=','), '\n',
			'UMI Counts: ', format(sum(umiData), big.mark=','), '\n',
			'Saturation: ', overall
		)) + ggtitle('Library Saturation')
	)

	return(overall)
}

#' @title Plot Library Saturation By Marker
#'
#' @description Create a plot of the library saturation per cell, separated by marker
#' @param citeseqCountDir The root of the Cite-seq-Count output folder, which should contain umi_count and read_count folders.
#' @export
PlotLibrarySaturationByMarker <- function(citeseqCountDir) {
	countData <- Seurat::Read10X(paste0(citeseqCountDir, '/read_count'), gene.column=1, strip.suffix = TRUE)
	countData <- countData[rownames(countData) != 'unmapped',]

	umiData <- Seurat::Read10X(paste0(citeseqCountDir, '/umi_count'), gene.column=1, strip.suffix = TRUE)
	umiData <- umiData[rownames(umiData) != 'unmapped',]

	saturation <- 1 - (umiData / countData)

	umiData <- reshape2::melt(t(as.matrix(umiData)))
	names(umiData) <- c('CellBarcode', 'Barcode', 'Value')
	umiData$Type <- 'UMI Count'

	saturation <- reshape2::melt(t(as.matrix(saturation)))
	names(saturation) <- c('CellBarcode', 'Barcode', 'Value')
	saturation$Type <- 'Saturation'

	umiData$Barcode <- SimplifyHtoNames(umiData$Barcode)
	saturation$Barcode <- SimplifyHtoNames(saturation$Barcode)
	saturation <- saturation[is.finite(saturation$Value),]

	P1 <- ggplot(umiData[umiData$Value > 0,], aes(x = Barcode, y = Value)) +
		labs(y = 'UMI Counts/Cell', x = '') +
		ggtitle('UMI Counts/Cell By Marker') +
		egg::theme_presentation(base_size = 14) +
		geom_boxplot() +
		scale_y_continuous(trans = scales::log2_trans()) +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	P2 <- ggplot(saturation, aes(x = Barcode, y = Value)) +
		labs(y = 'Saturation', x = '') +
		egg::theme_presentation(base_size = 14) +
		ggtitle('Saturation/Cell By Marker') +
		geom_boxplot() +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	print(P1)
	print(P2)
}


#' @title Calculate Saturation For 10x
#' @description Calculate per-cell saturation for 10x data
#' @param barcodeMatrix The matrix holding HTO count data. The columns should match cells
#' @param molInfoFile The 10x molecule_info.h5 file
#' @param doPlot If true, plots summarizing saturation will be generated.
#' @return A dataframe with the cellbarcode and per-cell saturation
#' @export
CalculateSaturationFor10x <- function(barcodeMatrix, molInfoFile, doPlot = TRUE) {
	.LogProgress('Calculating saturation')
	df <- DropletUtils::get10xMolInfoStats(molInfoFile)
	df$cellbarcode <- df$cell

	if (length(intersect(colnames(barcodeMatrix), df$cellbarcode)) == 0) {
		bc1 <- paste0(head(df$cellbarcode, n = 2), collapse = ';')
		bc2 <- paste0(head(colnames(barcodeMatrix), n = 2), collapse = ';')
		print(paste0('No overlapping barcodes found (example: ', bc1, ' / ', bc2,'), adding gem_group'))
		df$cellbarcode <- paste0(df$cellbarcode, '-', df$gem_group)
	}

	if (length(intersect(colnames(barcodeMatrix), df$cellbarcode)) == 0) {
		bc1 <- paste0(head(df$cellbarcode, n = 2), collapse = ';')
		bc2 <- paste0(head(colnames(barcodeMatrix), n = 2), collapse = ';')
		stop(paste0('No overlapping barcodes found between barcodeMatrix and molecule_info.h5 file, example: ', bc1, ' / ', bc2))
	}

	df <- data.frame(cellbarcode = df$cellbarcode, num.umis = df$num.umis, CountsPerCell = df$num.reads)
	df <- df[df$cellbarcode %in% colnames(barcodeMatrix),]
	df$Saturation <- 1 - (df$num.umis / df$CountsPerCell)

	if (doPlot) {
		overall <- 1 - round((sum(df$num.umis) / sum(df$CountsPerCell)), 2)
		print(ggplot(df, aes(x = CountsPerCell, y = Saturation)) +
				  labs(x = 'Counts/Cell', y = '% Saturation') +
				  egg::theme_presentation(base_size = 18) +
				  geom_point() +
				  annotate("text", x = max(df$CountsPerCell), y = min(df$Saturation), hjust = 1, vjust = -1, label = paste0(
					  'Total Counts: ', format(sum(df$CountsPerCell), big.mark=','), '\n',
					  'UMI Counts: ', format(sum(df$num.umis), big.mark=','), '\n',
					  'Saturation: ', overall
				  )) + ggtitle('Library Saturation')
		)
	}

	print(paste0('Adding saturation for ', nrow(df), ' cells'))

	dat <- data.frame(cellbarcode = colnames(barcodeMatrix), sortOrder = 1:ncol(barcodeMatrix))
	dat <- merge(dat, data.frame(cellbarcode = df$cellbarcode, saturation = df$Saturation), by = 'cellbarcode', all.x = T)
	dat <- dplyr::arrange(dat, sortOrder)

	.LogProgress('Finished calculating saturation')
	return(dat[c('cellbarcode', 'saturation')])
}