#' @include Utils.R
#' @include Normalization.R
#' @include Visualization.R


GenerateCellHashCallsDropletUtils <- function(barcodeMatrix, verbose = TRUE, assay = 'HTO', methodName = 'dropletutils', label = 'DropletUtils hashedDrops', runEmptyDrops = FALSE) {
	if (verbose) {
		print(paste0('Starting ', label))
	}

	tryCatch({
		seuratObj <- Seurat::CreateSeuratObject(barcodeMatrix, assay = assay)
		seuratObj <- ThresholdHashedDrops(seuratObj = seuratObj, assay = assay, columnSuffix = 'dropletutils', runEmptyDrops = runEmptyDrops)

		SummarizeHashingCalls(seuratObj, label = label, columnSuffix = 'dropletutils', assay = assay, doTSNE = FALSE, doHeatmap = FALSE)

		df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = methodName, classification = seuratObj$classification.dropletutils, classification.global = seuratObj$classification.global.dropletutils, stringsAsFactors = FALSE)
		return(df)
	}, error = function(e){
		print('Error generating hashedDrops calls, aborting')
		print(e)

		return(NULL)
	})
}

ThresholdHashedDrops <- function(seuratObj, assay, columnSuffix, runEmptyDrops = FALSE, seed = 1234) {
	barcodeMatrix <- Seurat::GetAssayData(
		object = seuratObj,
		assay = assay,
		slot = 'data'
	)[, colnames(x = seuratObj)]

	set.seed(seed)

	if (runEmptyDrops) {
		hash.calls <- DropletUtils::emptyDrops(barcodeMatrix, by.rank=min(ncol(barcodeMatrix), 40000))
		is.cell <- which(hash.calls$FDR <= 0.001)

		r <- rank(-hash.calls$Total)
		P1 <- ggplot(data.frame(Rank = r, Total = hash.calls$Total), aes(x = Rank, y = Total)) +
			geom_point() +
			scale_x_continuous(trans = 'log10') +
			scale_y_continuous(trans = 'log10') +
			xlab('Rank') + ylab('Total HTO Count') +
			egg::theme_presentation(base_size = 12)

		P2 <- ggplot(data.frame(Count = hash.calls$Total[is.cell]), aes(x = Count)) +
			geom_histogram(fill = 'grey', color = 'black', bins = 20) +
			scale_x_continuous(trans = 'log10') +
			xlab('Total HTO Count (log10)') + ylab('Count') +
			egg::theme_presentation(base_size = 12)

		print(P1 | P2)

		hash.stats <- DropletUtils::hashedDrops(barcodeMatrix[,is.cell], ambient=S4Vectors::metadata(hash.calls)$ambient)
	} else {

		ambience <- DropletUtils::inferAmbience(barcodeMatrix)
		print('Inferred Ambience:')
		print(ambience)

		hash.stats <- DropletUtils::hashedDrops(barcodeMatrix, ambient = ambience)
	}

	df <- data.frame(LogFC = hash.stats$LogFC, LogFC2 = hash.stats$LogFC2)
	df$Color <- rep("#a9a9a9", nrow(hash.stats))
	df$Color[hash.stats$Doublet] <- "red"
	df$Color[hash.stats$Confident] <- "black"
	df$Color <- as.factor(df$Color)

	df$Category <- rep("Negative", nrow(hash.stats))
	df$Category[hash.stats$Doublet] <- "Doublet"
	df$Category[hash.stats$Confident] <- "Singlet"
	df$Category <- as.factor(df$Category)

	P1 <- ggplot(df, aes(x = LogFC)) +
		geom_histogram(bins = 40) +
		egg::theme_presentation(base_size = 12) +
		xlab("Best to Second HTO (LogFC)") +
		xlab("Log fold-change") +
		ylab('Count')

	P2 <- ggplot(df, aes(x = LogFC, y = LogFC2, color = Category)) +
		geom_point() +
		egg::theme_presentation(base_size = 12) +
		xlab("Best to Second HTO (LogFC)") +
		ylab("Second HTO / Ambient (LogFC)") +
		scale_color_manual(values = levels(df$Color))

	suppressWarnings(print(P1 + P2))

	# This is not especially efficient, but reuses code:
	discrete <- GetAssayData(object = seuratObj, assay = assay)
	discrete[discrete > 0] <- 0

	for (htoIdx in 1:nrow(barcodeMatrix)) {
		# Count singlets or doublets:
		sel <- (!is.na(hash.stats$Best) & hash.stats$Best == htoIdx & hash.stats$Confident) | (!is.na(hash.stats$Doublet) & hash.stats$Doublet & (hash.stats$Best == htoIdx | hash.stats$Second == htoIdx))
		if (sum(sel) > 0) {
			cells <- rownames(hash.stats)[sel]
			discrete[htoIdx, cells] <- 1
		}
	}

	seuratObj <- .AssignCallsToMatrix(seuratObj, discrete, suffix = columnSuffix, assay = assay)

	return(seuratObj)
}