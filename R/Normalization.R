#' @include Utils.R

utils::globalVariables(
	names = c('NormCount', 'Saturation', 'Cluster', 'AvgExpression'),
	package = 'cellhashR',
	add = TRUE
)

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

#' @title PlotNormalizationQC
#'
#' @param barcodeData The count matrix
#' @description Generates QC plots related to normalization
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

	for (norm in names(toQC)) {
		PerformHashingClustering(toQC[[norm]], norm = norm)
	}
}

PerformHashingClustering <- function(barcodeMatrix, norm) {
	seuratObj <- CreateSeuratObject(barcodeMatrix, assay = norm)

	# Calculate tSNE embeddings with a distance matrix
	success <- FALSE
	tryCatch({
		perplexity <- .InferPerplexityFromSeuratObj(seuratObj, 100)
		print(paste0('Running tSNE with perplexity: ', perplexity, ' for normalization: ', norm))
		seuratObj[['hto_tsne']] <- RunTSNE(stats::dist(Matrix::t(barcodeMatrix)), assay = norm, perplexity = perplexity)

		.PlotClusters(barcodeMatrix, seuratObj, norm)
	}, error = function(e){
		print(e)
		print('Error generating tSNE, skipping')
	})
}

.PlotClusters <- function(barcodeMatrix, seuratObj, norm) {
	# clara:
	ncenters <- (nrow(x = barcodeMatrix) + 1)
	init.clusters <- clara(
		x = t(x = barcodeMatrix),
		k = ncenters,
		samples = 100
	)
	Idents(object = seuratObj, cells = names(x = init.clusters$clustering), drop = TRUE) <- as.character(init.clusters$clustering)
	seuratObj$cluster.clara <- Idents(seuratObj)
	P <- DimPlot(seuratObj, reduction = 'hto_tsne', group.by = 'cluster.clara', label = TRUE)
	P <- P + ggtitle(paste0('Clusters: ', norm, ' (clara)'))

	P2 <- .CreateClusterTilePlot(seuratObj, assay = norm)
	print(P | P2)

	# kmeans:
	init.clusters <- stats::kmeans(
		x = t(x = barcodeMatrix),
		centers = ncenters,
		nstart = 100
	)
	Idents(object = seuratObj, cells = names(x = init.clusters$cluster), drop = TRUE) <- as.character(init.clusters$cluster)
	seuratObj$cluster.kmeans <- Idents(seuratObj)

	P <- DimPlot(seuratObj, group.by = 'cluster.kmeans', label = TRUE)
	P <- P + ggtitle(paste0('Clusters: ', norm, ' (kmeans)'))
	P2 <- .CreateClusterTilePlot(seuratObj, assay = norm)
	print(P | P2)
}

.CreateClusterTilePlot <- function(seuratObj, assay) {
	average.expression <- AverageExpression(
		object = seuratObj,
		assays = c(assay),
		verbose = FALSE,
		return.seurat = TRUE
	)

	df <- reshape2::melt(GetAssayData(average.expression, assay = assay))
	names(df) <- c('Cluster', 'Barcode', 'AvgExpression')
	df$Cluster <- naturalsort::naturalfactor(df$Cluster)
	df$Barcode <- naturalsort::naturalfactor(df$Barcode)
	df$AvgExpression <- round(df$AvgExpression, 2)

	P2 <- ggplot(df, aes(Cluster, Barcode)) +
		geom_tile(aes(fill = AvgExpression), colour = "white") +
		geom_text(aes(label=AvgExpression)) +
		scale_fill_gradient2(low = "red", mid = "white", high = "green") +
		scale_x_discrete(position = "top") +
		labs(x = "Barcode",y = "Cluster") +
		egg::theme_presentation(base_size = 12) +
		theme(
			legend.position = 'none'
		)

	return(P2)
}