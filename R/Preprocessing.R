#' @title ProcessCiteSeqCount
#'
#' @description The primary entrypoint for parsing and QC of the cell hashing count matrix.
#' @param rawCountData, The input barcode file.
#' @param minCountPerCell Cells (columns) will be dropped if their total count is less than this value.
#' @import ggplot2
#' @import patchwork
#' @import utils
#' @import stats
#' @return
#' @export
ProcessCountMatrix <- function(rawCountData=NA, minCountPerCell = 5, barcodeWhitelist = NULL, barcodeBlacklist = c('no_match', 'total_reads', 'unmapped'), doPlot = T) {
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
		#CITE-seq-Count 1.4.2 and higher creates a folder
		barcodeData <- Seurat::Read10X(rawCountData, gene.column=1, strip.suffix = TRUE)
		barcodeData <- barcodeData[which(!(rownames(barcodeData) %in% barcodeBlacklist)), , drop = F]
		barcodeData <- as.matrix(barcodeData)
	} else {
		# older CITE-seq-Count versions created a CSV file, so support this:
		barcodeData <- utils::read.table(rawCountData, sep = ',', header = T, row.names = 1)
		barcodeData <- barcodeData[which(!(rownames(barcodeData) %in% barcodeBlacklist)),]
		barcodeData <- as.matrix(barcodeData)
	}

	print(paste0('Initial cell barcodes in hashing data: ', ncol(barcodeData)))

	# Print QC of counts by cell and barcode:
	if (doPlot) {
		PrintRowQc(barcodeData)
	}

	if (!is.null(barcodeWhitelist)) {
		print(paste0('Limiting to barcodes: ', paste0(barcodeWhitelist, collapse = ',')))
		sel <- rownames(barcodeData) %in% barcodeWhitelist
		if (sum(sel) == 0) {
			sel <- simplifyHtoNames(rownames(barcodeData)) %in% barcodeWhitelist
		}
		
		barcodeData <- barcodeData[sel,]
	}

	if (!is.null(minCountPerCell)) {
		barcodeData <- DoCellFiltering(barcodeData, minCountPerCell = minCountPerCell)
	}

	if (doPlot) {
		PrintColumnQc(barcodeData)
	}

	return(barcodeData)
}

PrintRowQc <- function(barcodeMatrix) {
	df <- GenerateByRowSummary(barcodeMatrix)
	df$Barcode <- naturalsort::naturalfactor(df$Barcode)
	df$Barcode <- forcats::fct_reorder(df$Barcode, df$rowSums, .desc = TRUE)

	# Total Counts
	P1 <- ggplot(df, aes(x = Barcode, y = rowSums)) +
		geom_bar(stat = 'identity') +
		ggtitle('Raw Counts/Barcode') +
		ylab('Row Sums') +
		egg::theme_presentation() +
		theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	P2 <- ggplot(df, aes(x = Barcode, y = rowSums)) +
		geom_bar(stat = 'identity') +
		egg::theme_presentation() +
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
		egg::theme_presentation() +
		theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	P2 <- ggplot(df, aes(x = Barcode, y = mean)) +
		geom_bar(stat = 'identity') +
		egg::theme_presentation() +
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
		egg::theme_presentation() +
		theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	P2 <- ggplot(df, aes(x = Barcode, y = mean_nonzero)) +
		geom_bar(stat = 'identity') +
		egg::theme_presentation() +
		ylab('Mean Counts/Cell (log1p)') +
		scale_y_continuous(trans = scales::log1p_trans()) +
		theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	print(P1 | P2)

	# Density:
	df <- data.frame(t(barcodeMatrix))
	colnames(df) <- simplifyHtoNames(rownames(barcodeMatrix))
	df <- tidyr::gather(df, Barcode, Count)

	P1 <- ggplot(df, aes(x = Count, color = Barcode)) +
		geom_density() +
		ggtitle('Counts/Cell') +
		egg::theme_presentation() +
		xlab('Counts/Cell') + ylab('Density') +
		scale_x_continuous(trans = scales::log1p_trans()) +
		theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	print(P1)

	# Density, trimmed:
	df <- data.frame(t(barcodeMatrix))
	colnames(df) <- simplifyHtoNames(rownames(barcodeMatrix))
	df <- tidyr::gather(df, Barcode, Count)
	out <- grDevices::boxplot.stats(df$Count)$out
	out <- out[out > mean(df$Count[df$Count > 0])]
	df <- df[df$Count < min(out),]

	P1 <- ggplot(df, aes(x = Count, color = Barcode)) +
		geom_density() +
		ggtitle('Counts/Cell, Outlier Trimmed') +
		egg::theme_presentation() +
		xlab('Counts/Cell') + ylab('Density') +
		scale_x_continuous(trans = scales::log1p_trans()) +
		theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)

	print(P1)

}

GenerateByRowSummary <- function(barcodeMatrix) {
	if (nrow(barcodeMatrix) == 0) {
		print('No rows in data, cannot generate summary')
		return(data.frame(Barcode = character()))
	}

	barcodeMatrix <- as.matrix(barcodeMatrix)
	df <- data.frame(Barcode = naturalsort::naturalfactor(simplifyHtoNames(rownames(barcodeMatrix))), BarcodeFull = naturalsort::naturalfactor(rownames(barcodeMatrix)), min = apply(barcodeMatrix, 1, min), max = apply(barcodeMatrix, 1, max), mean = apply(barcodeMatrix, 1, mean), logmean = log(apply(barcodeMatrix, 1, mean) + 1), nonzero = apply(barcodeMatrix, 1, function(x){
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

DoCellFiltering <- function(barcodeData, minQuant = 0.05, minCountPerCell = 5){
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
	names = c('Barcode1'),
	package = 'cellhashR',
	add = TRUE
)

PrintColumnQc <- function(barcodeMatrix) {
	df <- data.frame(t(barcodeMatrix))
	colnames(df) <- simplifyHtoNames(rownames(barcodeMatrix))
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
		egg::theme_presentation() 
	
	print(P1)
	
	out <- grDevices::boxplot.stats(df$Count)$out
	out <- out[out > mean(df$Count[df$Count > 0])]
	toPlot <- toPlot[toPlot$Value < min(out),]
	
	P1 <- ggplot(toPlot, aes(x = Value, y = Count, color = Barcode)) +
		geom_point() +
		ggtitle('Counts/Cell') +
		xlab('Count') + ylab('Total Cells') +
		egg::theme_presentation() 
	
	print(P1)


	#normalize columns, print top barcode fraction:
	normalizedBarcodes <- sweep(barcodeMatrix,2,colSums(barcodeMatrix),`/`)
	topValue <- apply(normalizedBarcodes,2,function(x){
		max(x)
	})

	df <- data.frame(Barcode1 = topValue)
	print(ggplot(df, aes(x = Barcode1)) +
		geom_histogram(binwidth = 0.05) +
		egg::theme_presentation() +
		xlab('Top Barcode Fraction') +
		ylab('Count')
	)

	print(paste0('Total cells where top barcode is >0.75 of counts: ', length(topValue > 0.75)))
	
}
