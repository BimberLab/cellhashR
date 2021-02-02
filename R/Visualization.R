#' @include Utils.R

utils::globalVariables(
	names = c('Classification'),
	package = 'cellhashR',
	add = TRUE
)

SummarizeHashingCalls <- function(seuratObj, label, columnSuffix, doHeatmap = T, doTSNE = T, assay = 'HTO') {
	htoClassificationField = paste0('classification.', columnSuffix)
	globalClassificationField <- paste0('classification.global.', columnSuffix)

	#report outcome
	df <- data.frame(prop.table(table(Barcode = seuratObj[[htoClassificationField]])))
	df$Barcode <- naturalsort::naturalfactor(df$Barcode)
	P1 <- ggplot(df, aes(x = '', y=Freq, fill=Barcode)) +
		geom_bar(width = 1, stat = "identity", color = "black") +
		coord_polar("y", start=0) +
		scale_fill_manual(values = GetPlotColors(length(unique(unlist(seuratObj[[htoClassificationField]]))))) +
		theme_minimal() +
		theme(
			axis.text.x=element_blank(),
			axis.title = element_blank(),
			axis.ticks = element_blank(),
			panel.grid  = element_blank()
		)

	df <- data.frame(prop.table(table(Classification = seuratObj[[globalClassificationField]])))
	P2 <- ggplot(df, aes(x = '', y=Freq, fill=Classification)) +
		geom_bar(width = 1, stat = "identity", color = "black") +
		coord_polar("y", start=0) +
		scale_fill_manual(values = GetPlotColors(length(unique(unlist(seuratObj[[globalClassificationField]]))))) +
		theme_minimal() +
		theme(
		axis.text.x=element_blank(),
		axis.title = element_blank(),
		axis.ticks = element_blank(),
		panel.grid  = element_blank()
		)

	print(P1 + P2 + plot_annotation(title = label))

	# Group cells based on the max HTO signal
	htos <- rownames(GetAssayData(seuratObj, assay = assay))
	for (hto in naturalsort::naturalsort(htos)){
		print(VlnPlot(seuratObj, group.by = htoClassificationField, features = c(hto), assay = assay, ncol = 1, log = F) + ggtitle(paste0(label, ": ", hto)))
	}

	if (doTSNE) {
		perplexity <- .InferPerplexityFromSeuratObj(seuratObj, 100)
		tryCatch({
			seuratObj[['hto_tsne']] <- RunTSNE(stats::dist(Matrix::t(GetAssayData(seuratObj, slot = "data", assay = assay))), assay = assay, perplexity = perplexity)
			P1 <- DimPlot(seuratObj, reduction = 'hto_tsne', group.by = htoClassificationField, label = FALSE)
			P2 <- DimPlot(seuratObj, reduction = 'hto_tsne', group.by = globalClassificationField, label = FALSE)

			print(P1 + P2  + plot_annotation(label))
		}, error = function(e){
			print('Error generating tSNE, skipping')
			print(conditionMessage(e))
			traceback()
		})
	}

	if (doHeatmap) {
		print(HTOHeatmap(seuratObj, assay = assay, classification = htoClassificationField, global.classification = globalClassificationField, ncells = min(4000, ncol(seuratObj)), singlet.names = NULL) + ggtitle(label))
	}
}
