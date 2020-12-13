
NormalizeLog2 <- function(mat, mean.center = TRUE) {
	log2Scaled <- log2(mat)
	for (i in 1:nrow(log2Scaled)) {
		ind <- which(is.finite(log2Scaled[i,]) == FALSE)
		log2Scaled[i,ind] <- 0

		if (mean.center) {
			log2Scaled[i,] <- log2Scaled[i,] - mean(log2Scaled[i,])
		}
	}

	return(as.matrix(log2Scaled))
}

NormalizeCLR <- function(mat) {
	seuratObj <- Seurat::CreateSeuratObject(mat, assay = 'Hashing')
	seuratObj <- Seurat::NormalizeData(seuratObj, assay = 'Hashing', normalization.method = "CLR", verbose = FALSE)

	return(seuratObj@assays$Hashing@data)
}

NormalizeRelative <- function(mat) {
	return(prop.table(mat, 2))
}

#' @export
PlotNormalizationQC <- function(barcodeData) {
	toQC <- list(
		'log2Center' = NormalizeLog2(barcodeData, mean.center = TRUE),
		'CLR' = NormalizeCLR(barcodeData),
		'relative' = NormalizeRelative(barcodeData)
	)

	df <- NULL
	for (norm in names(toQC)) {
		toAdd <- reshape2::melt(t(toQC[[norm]]))
		names(toAdd) <- c('CellBarcode', 'Barcode', 'NormCount')
		toAdd$Normalization <- norm

		if (is.null(df)) {
			df <- toAdd
		} else {
			df <- rbind(toAdd, df)
		}

		df$Barcode <- SimplifyHtoNames(as.character(df$Barcode))
	}

	print(ggplot2::ggplot(df, aes(x = NormCount, color = Barcode)) +
		egg::theme_presentation(base_size = 14) +
		geom_density(size = 1) + labs(y = 'Density', x = 'Value') + ggtitle('Normalized Data') +
		facet_wrap(Barcode ~ Normalization, scales = 'free', ncol = length(unique(df$Normalization)), strip.position = 'top', labeller = labeller(.multi_line = FALSE))
	)
}