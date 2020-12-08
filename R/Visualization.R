#' @title A Title
#'
#' @description A description
#' @param seuratObj, A Seurat object.
#' @importFrom cluster clara
#' @importFrom Matrix t
#' @return A modified Seurat object.
DebugDemux <- function(seuratObj, assay = 'HTO', reportKmeans = FALSE) {
	print('Debugging information for Seurat HTODemux:')
	data <- GetAssayData(object = seuratObj, assay = assay)
	ncenters <- (nrow(x = data) + 1)

	init.clusters <- clara(
	x = t(x = data),
	k = ncenters,
	samples = 100
	)
	Idents(object = seuratObj, cells = names(x = init.clusters$clustering), drop = TRUE) <- init.clusters$clustering

	# Calculate tSNE embeddings with a distance matrix
	tryCatch({
		perplexity <- .InferPerplexityFromSeuratObj(seuratObj, 100)
		seuratObj[['hto_tsne']] <- RunTSNE(stats::dist(Matrix::t(data)), assay = assay, perplexity = perplexity)
		P <- DimPlot(seuratObj, reduction = 'hto_tsne', label = TRUE)
		P <- P + ggtitle('Clusters: clara')
		print(P)
	}, error = function(e){
		print(e)
		print('Error generating tSNE, skipping')
	})

	average.expression <- AverageExpression(
	object = seuratObj,
	assays = c(assay),
	verbose = FALSE
	)[[assay]]

	print(kableExtra::kbl(average.expression, label = 'clara') %>% kableExtra::kable_styling())

	if (reportKmeans) {
		print('kmeans:')
		init.clusters <- stats::kmeans(
		x = t(x = data),
		centers = ncenters,
		nstart = 100
		)
		Idents(object = seuratObj, cells = names(x = init.clusters$cluster), drop = TRUE) <- init.clusters$cluster

		# Calculate tSNE embeddings with a distance matrix
		P <- DimPlot(seuratObj, label = TRUE)
		P <- P + ggtitle('Clusters: kmeans')
		print(P)

		average.expression <- AverageExpression(
		object = seuratObj,
		assays = c(assay),
		verbose = FALSE
		)[[assay]]

		print(kableExtra::kbl(average.expression, label = 'kmeans') %>% kableExtra::kable_styling())
	}
}

PlotHtoCountData <- function(seuratObj, label, assay = 'HTO') {
	#Plot raw data by HTO:
	data <- GetAssayData(seuratObj, assay = assay, slot = 'counts')
	df2 <- as.data.frame(data)
	df2$HTO <- row.names(data)
	df2 <- tidyr::gather(df2, key = 'CellBarcode', value = 'Count', -HTO)
	df2$HTO <- simplifyHtoNames(as.character(df2$HTO))
	df2$HTO <- naturalsort::naturalfactor(df2$HTO)

	df2$Count <- log10(df2$Count + 0.5)
	print(ggplot(df2, aes(y = Count, x = HTO, fill = HTO)) +
		geom_boxplot() +
		ggtitle(paste0(label, ': Raw counts by HTO (log10)'))
	)

	#Plot normalized data by HTO:
	data <- GetAssayData(seuratObj, assay = assay, slot = 'data')
	df2 <- as.data.frame(data)
	df2$HTO <- row.names(data)
	df2 <- tidyr::gather(df2, key = 'CellBarcode', value = 'Count', -HTO)
	df2$HTO <- simplifyHtoNames(as.character(df2$HTO))
	df2$HTO <- naturalsort::naturalfactor(df2$HTO)

	print(ggplot(df2, aes(y = Count, x = HTO, fill = HTO)) +
		geom_boxplot() +
		ggtitle(paste0(label, ': Normalized data by HTO'))
	)
}

#' @title A Title
#'
#' @description A description
#' @param seuratObj, A Seurat object.
#' @return A modified Seurat object.
HtoSummary <- function(seuratObj, htoClassificationField, globalClassificationField, label, doHeatmap = T, doTSNE = T, assay = 'HTO') {
	#report outcome
	print(table(seuratObj[[htoClassificationField]]))
	print(table(seuratObj[[globalClassificationField]]))

	# Group cells based on the max HTO signal
	seuratObj_hashtag <- seuratObj
	Idents(seuratObj_hashtag) <- htoClassificationField
	htos <- rownames(GetAssayData(seuratObj_hashtag, assay = assay))
	for (hto in naturalsort::naturalsort(htos)){
		print(VlnPlot(seuratObj_hashtag, features = c(hto), assay = assay, ncol = 1, log = F) + ggtitle(paste0(label, ": ", hto)))
	}

	if (doTSNE) {
		perplexity <- .InferPerplexityFromSeuratObj(seuratObj, 100)
		tryCatch({
			seuratObj[['hto_tsne']] <- RunTSNE(stats::dist(Matrix::t(GetAssayData(seuratObj, slot = "data", assay = assay))), assay = assay, perplexity = perplexity)
			print(DimPlot(seuratObj, reduction = 'hto_tsne', group.by = htoClassificationField, label = FALSE) + ggtitle(label))
			print(DimPlot(seuratObj, reduction = 'hto_tsne', group.by = globalClassificationField, label = FALSE) + ggtitle(label))
		}, error = function(e){
			print(e)
			print('Error generating tSNE, skipping')
		})
	}

	if (doHeatmap) {
		print(HTOHeatmap(seuratObj, assay = assay, classification = htoClassificationField, global.classification = globalClassificationField, ncells = min(3000, ncol(seuratObj)), singlet.names = NULL) + ggtitle(label))
	}
}
