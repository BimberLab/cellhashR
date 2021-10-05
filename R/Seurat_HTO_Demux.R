#' @include Utils.R
#' @include Visualization.R

utils::globalVariables(
	names = c('p_val_adj', 'avg_logFC', 'cluster'),
	package = 'cellhashR',
	add = TRUE
)

GenerateCellHashCallsSeurat <- function(barcodeMatrix, positive.quantile = 0.95, methodName = 'htodemux', verbose= TRUE, metricsFile = NULL, doTSNE = TRUE, doHeatmap = TRUE) {
	if (verbose) {
		print('Starting HTODemux')
	}

	seuratObj <- suppressWarnings(CreateSeuratObject(barcodeMatrix, assay = 'HTO'))

	tryCatch({
		seuratObj <- DoHtoDemux(seuratObj, positive.quantile = positive.quantile, verbose = verbose, metricsFile = metricsFile, doTSNE = doTSNE, doHeatmap = doHeatmap)

		df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = methodName, classification = seuratObj$classification.htodemux, classification.global = seuratObj$classification.global.htodemux, stringsAsFactors = FALSE)
		df <- .RestoreUnderscoreToHtoNames(df, rownames(barcodeMatrix))
		return(df)
	}, error = function(e){
		print('Error generating seurat htodemux calls, aborting')
		if (!is.null(e)) {
			print(conditionMessage(e))
			traceback()
		}

		return(NULL)
	})
}


DoHtoDemux <- function(seuratObj, positive.quantile, label = 'Seurat HTODemux', plotDist = FALSE, verbose = TRUE, metricsFile = NULL, doTSNE = TRUE, doHeatmap = TRUE) {
	# Normalize HTO data, here we use centered log-ratio (CLR) transformation
	seuratObj <- NormalizeData(seuratObj, assay = "HTO", normalization.method = "CLR", verbose = FALSE)

	seuratObj <- HTODemux(seuratObj, positive.quantile =  positive.quantile, plotDist = plotDist, verbose = verbose, metricsFile = metricsFile)

	SummarizeHashingCalls(seuratObj, label = label, columnSuffix = 'htodemux', doTSNE = doTSNE, doHeatmap = doHeatmap)

	return(seuratObj)
}


#' @import Seurat
#' @importFrom fitdistrplus fitdist
#' @importFrom cluster clara
#' @importFrom Matrix t
#' @author Seurat
#' url https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/HTODemux
HTODemux <- function(
	object,
	assay = "HTO",
	positive.quantile = 0.95,
	nstarts = 100,
	kfunc = "clara",
	nsamples = 100,
	verbose = TRUE,
	plotDist = FALSE,
	metricsFile = NULL
) {
	#initial clustering
	data <- GetAssayData(object = object, assay = assay)
	counts <- GetAssayData(
		object = object,
		assay = assay,
		slot = 'counts'
	)[, colnames(x = object)]

	ncenters <- (nrow(x = data) + 1)
	switch(
		EXPR = kfunc,
		'kmeans' = {
			init.clusters <- stats::kmeans(
				x = t(x = data),
				centers = ncenters,
				nstart = nstarts
			)
			#identify positive and negative signals for all HTO
			Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
		},
		'clara' = {
			#use fast k-medoid clustering
			init.clusters <- clara(
				x = t(x = data),
				k = ncenters,
				samples = nsamples
			)
			#identify positive and negative signals for all HTO
			Idents(object = object, cells = names(x = init.clusters$clustering), drop = TRUE) <- init.clusters$clustering
		},
		stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'")
	)

	#average hto signals per cluster
	#work around so we don't average all the RNA levels which takes time
	average.expression <- suppressWarnings(Seurat::AverageExpression(
		object = object,
		assays = c(assay),
		slot = 'counts',
		verbose = FALSE
	))[[assay]]

	#create a matrix to store classification result
	discrete <- GetAssayData(object = object, assay = assay)
	discrete[discrete > 0] <- 0
	# for each HTO, we will use the minimum cluster for fitting
	thresholds <- list()
	for (hto in naturalsort::naturalsort(rownames(x = data))) {
		values <- counts[hto, colnames(object)]

		# Take the bottom 2 clusters (top 2 assumed to be HTO and doublet) as background.
		maxPossibleBackgroundCols <- max(nrow(data) - 2, 1)
		numBackgroundCols <- min(2, maxPossibleBackgroundCols)
		backgroundIndices <- order(average.expression[hto, ])[1:numBackgroundCols]

		if (sum(average.expression[hto, backgroundIndices]) == 0) {
			allPossibleBackgroundIndices <- order(average.expression[hto, ])[1:maxPossibleBackgroundCols]
			for (i in 1:maxPossibleBackgroundCols) {
				print('Expanding clusters until non-zero background obtained')
				backgroundIndices <- allPossibleBackgroundIndices[1:i]
				if (sum(average.expression[hto, backgroundIndices]) > 0) {
					break
				}
			}
		}

		if (verbose) {
			print(paste0('Will select bottom ', numBackgroundCols, ' barcodes as background'))
			print(paste0('Background clusters for ', hto, ': ', paste0(backgroundIndices, collapse = ',')))
		}

		if (sum(average.expression[hto, backgroundIndices]) == 0) {
			stop('The background clusters have zero reads, cannot call')
		}

		cutoff <- NULL
		values.use <- values[WhichCells(
			object = object,
			idents = levels(x = Idents(object = object))[backgroundIndices]
		)]

		if (verbose) {
			print(paste0('total cells for background: ', length(values.use)))
		}

		tryCatch(expr = {
			fit <- suppressWarnings(fitdist(data = values.use, distr = "nbinom"))
			if (plotDist) {
				print(plot(fit))
			}

			cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
		}, error = function(e) {
			saveRDS(values.use, file = paste0('./', hto, '.fail.nbinom.rds'))
		})

		if (is.null(cutoff)) {
			print(paste0('Skipping HTO due to failure to fit distribution: ', hto))
			next
		}

		thresholds[[hto]] <- cutoff

		if (verbose) {
			print(paste0("Cutoff for ", hto, " : ", cutoff, " reads"))
		}
		discrete[hto, names(x = which(x = values > cutoff))] <- 1
	}

	# now assign cells to HTO based on discretized values
	object <- .AssignCallsToMatrix(object, discrete, suffix = 'htodemux', assay = assay)

	print("Thresholds:")
	for (hto in names(thresholds)) {
		print(paste0(hto, ': ', thresholds[[hto]]))
		.LogMetric(metricsFile, paste0('cutoff.htodemux.', hto), thresholds[[hto]])
	}

	return(object)
}