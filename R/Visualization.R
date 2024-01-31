#' @include Utils.R

utils::globalVariables(
	names = c('Classification'),
	package = 'cellhashR',
	add = TRUE
)

SummarizeHashingCalls <- function(seuratObj, label, columnSuffix, doHeatmap = T, doTSNE = T, assay = 'HTO') {
	.LogProgress('Summarizing Hashing Calls')
	if (is.null(seuratObj)) {
		print('No seurat object provided, cannot run SummarizeHashingCalls')
		return()
	}

	htoClassificationField <- paste0('classification.', columnSuffix)
	globalClassificationField <- paste0('classification.global.', columnSuffix)

	if (!(htoClassificationField %in% names(seuratObj@meta.data))) {
		print(paste0('Unable to find calls in metadata, expected field: ', htoClassificationField))
		return()
	}

	if (!(globalClassificationField %in% names(seuratObj@meta.data))) {
		print(paste0('Unable to find calls in metadata, expected field: ', globalClassificationField))
		return()
	}

	#report outcome
	df <- data.frame(prop.table(table(Barcode = seuratObj@meta.data[[htoClassificationField]])))
	df$Barcode <- naturalsort::naturalfactor(df$Barcode)
	P1 <- ggplot(df, aes(x = '', y=Freq, fill=Barcode)) +
		geom_bar(width = 1, stat = "identity", color = "black") +
		coord_polar("y", start=0) +
		scale_fill_manual(values = GetPlotColors(length(unique(seuratObj@meta.data[[htoClassificationField]])))) +
		theme_minimal() +
		theme(
			axis.text.x=element_blank(),
			axis.title = element_blank(),
			axis.ticks = element_blank(),
			panel.grid  = element_blank()
		)

	df <- data.frame(prop.table(table(Classification = seuratObj@meta.data[[globalClassificationField]])))
	P2 <- ggplot(df, aes(x = '', y=Freq, fill=Classification)) +
		geom_bar(width = 1, stat = "identity", color = "black") +
		coord_polar("y", start=0) +
		scale_fill_manual(values = GetPlotColors(length(unique(seuratObj@meta.data[[globalClassificationField]])))) +
		theme_minimal() +
		theme(
		axis.text.x=element_blank(),
		axis.title = element_blank(),
		axis.ticks = element_blank(),
		panel.grid  = element_blank()
		)

	suppressWarnings(print(P1 + P2 + plot_annotation(title = label)))

	# Group cells based on the max HTO signal
	htos <- rownames(GetAssayData(seuratObj, assay = assay))
	for (hto in naturalsort::naturalsort(htos)){
		suppressWarnings(print(VlnPlot(seuratObj, group.by = htoClassificationField, features = c(hto), assay = assay, ncol = 1, log = F) + ggtitle(paste0(label, ": ", hto))))
	}

	if (doTSNE) {
		.LogProgress('Running tSNE')
		perplexity <- .InferPerplexityFromSeuratObj(seuratObj, 100)
		tryCatch({
			assayObj <- GetAssay(seuratObj, assay = assay)
			slotName <- NULL
			if (class(assayObj)[1] == 'Assay') {
				slotName <- 'data'
			} else if (class(assayObj)[1] == 'Assay5') {
				slotName <- ifelse('data' %in% names(assayObj@layers), yes = 'data', no = 'counts')
			} else {
				stop(paste0('Unknown class: ', class(assayObj)[1]))
			}

			print(paste0('Using slot: ', slotName))
			mat <- GetAssayData(assayObj, slot = slotName)
			seuratObj[['hto_tsne']] <- suppressWarnings(RunTSNE(stats::dist(Matrix::t(mat)), perplexity = perplexity))
			P1 <- DimPlot(seuratObj, reduction = 'hto_tsne', group.by = htoClassificationField, label = FALSE)
			P2 <- DimPlot(seuratObj, reduction = 'hto_tsne', group.by = globalClassificationField, label = FALSE)

			print(P1 + P2  + plot_annotation(label))
		}, error = function(e){
			print('Error generating tSNE, skipping')
			print(conditionMessage(e))
			traceback()
		})
		.LogProgress('Done with tSNE')
	}

	if (doHeatmap) {
		tryCatch({
			P1 <- suppressWarnings(HTOHeatmap(seuratObj, assay = assay, classification = htoClassificationField, global.classification = globalClassificationField, ncells = min(4000, ncol(seuratObj)), singlet.names = NULL))
			if (!is.null(P1)) {
				P1 <- P1 + ggtitle(label)
				suppressWarnings(print(P1))
			}
		}, error = function(e){
			print('Error creating heatmap, skipping')
			print(conditionMessage(e))
			traceback()
		})
	}
}
