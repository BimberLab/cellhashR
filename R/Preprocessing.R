#' @include Utils.R

#' @title Process CiteSeq Count Matrix
#'
#' @description The primary entrypoint for parsing and QC of the cell hashing count matrix.
#' @param rawCountData, The input barcode file or umi_count folder
#' @param minCountPerCell Cells (columns) will be dropped if their total count is less than this value.
#' @param barcodeWhitelist A vector of barcode names to retain.
#' @param barcodeBlacklist A vector of barcodes names to discard.
#' @param doPlot If true, QC plots will be generated
#' @param simplifyBarcodeNames If true, the sequence tag portion will be removed from the barcode names (i.e. HTO-1-ATGTGTGA -> HTO-1)
#' @param saveOriginalCellBarcodeFile An optional file path, where the set of original cell barcodes, prior to filtering, will be written. The primary use-case is if the count matrix was generated using a cell whitelist (like cells with passing gene expression). Preserving this list allows downstream reporting.
#' @param metricsFile If provided, summary metrics will be written to this file.
#' @import ggplot2
#' @import patchwork
#' @import utils
#' @return The updated count matrix
#' @export
ProcessCountMatrix <- function(rawCountData=NA, minCountPerCell = 5, barcodeWhitelist = NULL, barcodeBlacklist = c('no_match', 'total_reads', 'unmapped'), doPlot = TRUE, simplifyBarcodeNames = TRUE, saveOriginalCellBarcodeFile = NULL, metricsFile = NULL) {
	barcodeData <- .LoadCountMatrix(rawCountData = rawCountData, barcodeBlacklist = barcodeBlacklist, simplifyBarcodeNames = simplifyBarcodeNames)

	print(paste0('Initial cell barcodes in hashing data: ', ncol(barcodeData)))
	.LogMetric(metricsFile, 'InitialCellBarcodes', ncol(barcodeData))

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
	}

	if (!is.null(minCountPerCell)) {
		barcodeData <- DoCellFiltering(barcodeData, minCountPerCell = minCountPerCell)
	}
	.LogMetric(metricsFile, 'PassingCellBarcodes', ncol(barcodeData))

	if (doPlot) {
		PrintColumnQc(barcodeData)
	}

	return(barcodeData)
}

.LoadCountMatrix <- function(rawCountData = NA, barcodeBlacklist = c('no_match', 'total_reads', 'unmapped'), simplifyBarcodeNames = TRUE) {
	if (is.na(rawCountData)){
		stop("No file set: change rawCountData")
	}

	if (is.na(barcodeBlacklist) || is.null(barcodeBlacklist)) {
		barcodeBlacklist <- character()
	}

	if (!file.exists(rawCountData)){
		stop(paste0("File does not exist: ", rawCountData))
	}

	if (dir.exists(rawCountData)) {
		barcodeData <- Seurat::Read10X(rawCountData, gene.column=1, strip.suffix = TRUE)
		barcodeData <- barcodeData[which(!(rownames(barcodeData) %in% barcodeBlacklist)), , drop = F]
		barcodeData <- as.matrix(barcodeData)
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
	names = c('Barcode', 'Value', 'CellBarcode', 'Freq', 'Str', 'mean_nonzero'),
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
	P1 <- ggplot(df, aes(x = Barcode, y = mean)) +
		geom_bar(stat = 'identity') +
		ggtitle('Mean Counts/Cell') +
		ylab('Mean Counts/Cell') +
		egg::theme_presentation(base_size = 14) +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	P2 <- ggplot(df, aes(x = Barcode, y = mean)) +
		geom_bar(stat = 'identity') +
		egg::theme_presentation(base_size = 14) +
		ylab('Mean Counts/Cell (log1p)') +
		scale_y_continuous(trans = scales::log1p_trans()) +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	print(P1 | P2)

	# Mean Counts, Non-zero
	P1 <- ggplot(df, aes(x = Barcode, y = mean_nonzero)) +
		geom_bar(stat = 'identity') +
		ggtitle('Mean Counts/Cell') +
		ylab('Mean Counts/Cell') +
		egg::theme_presentation(base_size = 14) +
		theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	P2 <- ggplot(df, aes(x = Barcode, y = mean_nonzero)) +
		geom_bar(stat = 'identity') +
		egg::theme_presentation(base_size = 14) +
		ylab('Mean Counts/Cell (log1p)') +
		scale_y_continuous(trans = scales::log1p_trans()) +
		theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	print(P1 | P2)

	# Density:
	df <- data.frame(t(barcodeMatrix))
	colnames(df) <- SimplifyHtoNames(rownames(barcodeMatrix))
	df <- tidyr::gather(df, Barcode, Count)

	P1 <- ggplot(df[df$Count > 0,], aes(x = Count, color = Barcode)) +
		geom_density() +
		ggtitle('Non-zero Counts/Cell') +
		egg::theme_presentation(base_size = 18) +
		xlab('Counts/Cell') + ylab('Density') +
		scale_x_continuous(trans = scales::log1p_trans()) +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(0.5)),
		)

	# Density, trimmed:
	df <- data.frame(t(barcodeMatrix))
	colnames(df) <- SimplifyHtoNames(rownames(barcodeMatrix))
	df <- tidyr::gather(df, Barcode, Count)
	out <- grDevices::boxplot.stats(df$Count)$out
	out <- out[out > mean(df$Count[df$Count > 0])]
	df <- df[df$Count < min(out),]

	P2 <- ggplot(df[df$Count > 0,], aes(x = Count, color = Barcode)) +
		geom_density() +
		ggtitle('Non-zero Counts/Cell, Outlier Trimmed') +
		egg::theme_presentation(base_size = 18) +
		xlab('Counts/Cell') + ylab('Density') +
		scale_x_continuous(trans = scales::log1p_trans()) +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(0.5))
		)

	print(P1)
	print(P2)
}

GenerateByRowSummary <- function(barcodeMatrix) {
	if (nrow(barcodeMatrix) == 0) {
		print('No rows in data, cannot generate summary')
		return(data.frame(Barcode = character()))
	}

	barcodeMatrix <- as.matrix(barcodeMatrix)
	df <- data.frame(Barcode = naturalsort::naturalfactor(SimplifyHtoNames(rownames(barcodeMatrix))), BarcodeFull = naturalsort::naturalfactor(rownames(barcodeMatrix)), min = apply(barcodeMatrix, 1, min), max = apply(barcodeMatrix, 1, max), mean = apply(barcodeMatrix, 1, mean), logmean = log(apply(barcodeMatrix, 1, mean) + 1), nonzero = apply(barcodeMatrix, 1, function(x){
		sum(x > 0)
	}), mean_nonzero = (rowSums(barcodeMatrix) / rowSums(!!barcodeMatrix)), total_gt1 = apply(barcodeMatrix, 1, function(x){
		sum(x > 1)
	}), mean_gt1 = apply(barcodeMatrix, 1, function(x){
		mean(sapply(x, function(y){if (y > 1) y else NA}), na.rm = T)
	}), rowSums = rowSums(barcodeMatrix))

	df$mean_nonzero[is.na(df$mean_nonzero)] <- 0

	df <- df[order(df$Barcode),]

	return(df)
}

DoCellFiltering <- function(barcodeData, minCountPerCell = 5){
	if (!is.null(minCountPerCell)) {
		toDrop <- sum(colSums(barcodeData) < minCountPerCell)
		if (toDrop > 0){
			print(paste0('Cells dropped due to low total counts per cell (<', minCountPerCell, '): ', toDrop))
			barcodeData <- barcodeData[,which(colSums(barcodeData) >= minCountPerCell), drop = FALSE]
			print(paste0('After filter: ', ncol(barcodeData)))
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

	P1 <- ggplot(toPlot, aes(x = Value, y = Count, color = Barcode)) +
		geom_point() +
		ggtitle('Counts/Cell') +
		xlab('Count') + ylab('Total Cells') +
		egg::theme_presentation(base_size = 18)
	
	out <- grDevices::boxplot.stats(df$Count)$out
	out <- out[out > mean(df$Count[df$Count > 0])]
	toPlot <- toPlot[toPlot$Value < min(out),]
	
	P2 <- ggplot(toPlot, aes(x = Value, y = Count, color = Barcode)) +
		geom_point() +
		ggtitle('Counts/Cell, Outlier Trimmed') +
		xlab('Count') + ylab('Total Cells') +
		egg::theme_presentation(base_size = 18)
	
	print(P1)
	print(P2)


	#normalize columns, print top barcode fraction:
	normalizedBarcodes <- sweep(barcodeMatrix, 2, colSums(barcodeMatrix),`/`)
	topValue <- apply(normalizedBarcodes,2,function(x){
		max(x)
	})

	df <- data.frame(Barcode1 = topValue)
	P1 <- ggplot(df, aes(x = Barcode1)) +
		geom_histogram(binwidth = 0.05) +
		egg::theme_presentation() +
		xlab('Fraction') +
		ylab('# Cells') + ggtitle('Top Barcode Fraction Per Cell') +
		expand_limits(x = c(0, 1))

	P1 <- P1 + plot_annotation(caption = paste0('Total cells where top barcode is >0.75 of counts: ', sum(topValue > 0.75), ' of ', length(topValue))) & theme(plot.caption = element_text(size = 14))

	print(P1)

	#MA-plot
	if (nrow(barcodeMatrix) == 2){
		tryCatch({
			df <- data.frame(t(barcodeMatrix))
			M = log2(df[,1]) - log2(df[,2])
			A = (log2(df[,1]) + log2(df[,2]))/2

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
#' @param citeseqCountDir, The root of the Cite-seq-Count output folder, which should contain umi_count and read_count folders.
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
#' @param citeseqCountDir, The root of the Cite-seq-Count output folder, which should contain umi_count and read_count folders.
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