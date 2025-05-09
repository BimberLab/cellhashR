#' @include Utils.R
#' @include Visualization.R
#' @include Normalization.R

utils::globalVariables(
	names = c('Proportion', 'Subset'),
	package = 'cellhashR',
	add = TRUE
)

GenerateCellHashCallsMultiSeq <- function(barcodeMatrix, assay = 'HTO', autoThresh = TRUE, quantile = NULL, maxiter = 20, qrange = seq(from = 0.2, to = 0.95, by = 0.05), doRelNorm = FALSE, methodName = 'multiseq', label = 'Multiseq deMULTIplex', verbose = TRUE, metricsFile = NULL, doTSNE = TRUE, doHeatmap = TRUE) {
	if (verbose) {
		print(paste0('Starting ', label))
	}

	tryCatch({
		seuratObj <- CreateSeuratObject(Seurat::as.sparse(barcodeMatrix), assay = assay)
		if (doRelNorm) {
			print('Performing relative normalization instead of log2')
			ad <- NormalizeRelative(barcodeMatrix)
			# NOTE: perform this rename to match the rename Seurat will perform anyway:
			rownames(ad) <- gsub(rownames(ad), pattern = '_', replacement = '-')
			seuratObj <- SetAssayData4Or5(seuratObj, assay = assay, theLayer = 'data', new.data = ad)
		} else {
			ad <- NormalizeLog2(barcodeMatrix)
			rownames(ad) <- gsub(rownames(ad), pattern = '_', replacement = '-')
			seuratObj <- SetAssayData4Or5(seuratObj, assay = assay, theLayer = 'data', new.data = ad)
		}

		seuratObj <- MULTIseqDemux(seuratObj, assay = assay, quantile = quantile, verbose = verbose, autoThresh = autoThresh, maxiter = maxiter, qrange = qrange, metricsFile = metricsFile)

		SummarizeHashingCalls(seuratObj, label = label, columnSuffix = 'multiseq', assay = assay, doTSNE = doTSNE, doHeatmap = doHeatmap)

		df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = methodName, classification = seuratObj$classification.multiseq, classification.global = seuratObj$classification.global.multiseq, stringsAsFactors = FALSE)
		df <- .RestoreUnderscoreToHtoNames(df, rownames(barcodeMatrix))
		return(df)
	}, error = function(e){
		print('Error generating multiseq calls, aborting')
		print(conditionMessage(e))
		traceback()
		return(NULL)
	})
}

#' @rawNamespace import(Matrix, except = c('tail', 'head'))
#' @author Chris McGinnis, Gartner Lab, UCSF
#' @references \url{https://www.biorxiv.org/content/10.1101/387241v1}
MULTIseqDemux <- function(
	object,
	assay = "HTO",
	quantile = 0.7,
	autoThresh = FALSE,
	maxiter = 5,
	qrange = seq(from = 0.1, to = 0.9, by = 0.05),
	verbose = TRUE,
	metricsFile = NULL
) {
	if (is.na(assay) || is.null(assay)) {
		assay <- DefaultAssay(object = object)
	}

	if (!'data' %in% SeuratObject::Layers(Seurat::GetAssay(object, assay = assay))) {
		stop('Missing counts layer from multiseq object!')
	}

	multi_data_norm <- t(x = GetAssayData(
		object = object,
		layer = "data",
		assay = assay
	))

	thresholds <- list()
	if (autoThresh) {
		iter <- 1
		negatives <- c()
		neg.vector <- c()
		while (iter <= maxiter) {
			# Iterate over q values to find ideal barcode thresholding results by maximizing singlet classifications
			bar.table_sweep.list <- list()
			n <- 0
			for (q in qrange) {
				n <- n + 1
				# Generate list of singlet/doublet/negative classifications across q sweep
				cr <- ClassifyCells(data = multi_data_norm, q = q)
				bar.table_sweep.list[[n]] <- cr[['calls']]
				names(x = bar.table_sweep.list)[n] <- paste0("q=" , q)
			}

			# Determine which q values results in the highest pSinglet
			threshold.results1 <- FindThresh(call.list = bar.table_sweep.list)
			res_round <- threshold.results1$res
			res.use <- res_round[res_round$Subset == "pSinglet", ]
			q.use <- res.use[which.max(res.use$Proportion),"q"]

			cr <- ClassifyCells(data = multi_data_norm, q = q.use)
			round.calls <- cr[['calls']]
			thresholds <- cr[['thresholds']]
			#remove negative cells
			neg.cells <- names(x = round.calls)[which(x = round.calls == "Negative")]
			called.cells <- sum(round.calls != "Negative")
			neg.vector <- c(neg.vector, rep(x = "Negative", length(x = neg.cells)))
			negatives <- c(negatives, neg.cells)

			print(ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) +
				geom_line() +
				theme(
				legend.position = "right"
				) +
				geom_vline(xintercept=q.use, lty=2) +
				scale_color_manual(values=c("red","black","blue")) +
				ggtitle(paste0('Multi-Seq Iteration: ', iter, ', quantile: ', q.use, ', input: ', length(round.calls), ', called: ', called.cells))
			)

			if (length(x = neg.cells) == 0) {
				print(paste0('no negatives, breaking after iteration: ', iter))
				break
			}

			multi_data_norm <- multi_data_norm[-which(x = rownames(x = multi_data_norm) %in% neg.cells), ]
			iter <- iter + 1
		}
		names(x = neg.vector) <- negatives
		demux_result <- c(round.calls,neg.vector)
		demux_result <- demux_result[rownames(x = object[[]])]
	} else{
		cr <- ClassifyCells(data = multi_data_norm, q = quantile)
		demux_result <- cr[['calls']]
		thresholds <- cr[['thresholds']]
	}

	demux_result <- demux_result[rownames(x = object[[]])]
	object[['classification.multiseq']] <- factor(x = demux_result)
	Idents(object = object) <- "classification.multiseq"

	doublet.id <- which(x = demux_result == 'Doublet')
	classification.multiseq <- as.character(object$classification.multiseq)
	classification.multiseq[doublet.id] <- 'Doublet'
	object$classification.multiseq <- naturalsort::naturalfactor(x = classification.multiseq)

	object$classification.global.multiseq <- as.character(object$classification.multiseq)
	object$classification.global.multiseq[!(object$classification.multiseq %in% c('Negative', 'Doublet'))] <- 'Singlet'
	object$classification.global.multiseq <- naturalsort::naturalfactor(object$classification.global.multiseq)

	if (!is.null(thresholds)) {
		print("Thresholds:")
		for (hto in names(thresholds)) {
			print(paste0(hto, ': ', thresholds[[hto]]))
			.LogMetric(metricsFile, paste0('cutoff.multiseq.', hto), thresholds[[hto]])
		}
	}
	object@misc[['cutoffs']] <- thresholds

	return(object)
}

# Inter-maxima quantile sweep to find ideal barcode thresholds
# @author Chris McGinnis, Gartner Lab, UCSF
#
FindThresh <- function(call.list) {
	# require(reshape2)
	res <- as.data.frame(x = matrix(
	data = 0L,
	nrow = length(x = call.list),
	ncol = 4
	))
	colnames(x = res) <- c("q","pDoublet","pNegative","pSinglet")
	q.range <- unlist(x = strsplit(x = names(x = call.list), split = "q="))
	res$q <- as.numeric(x = q.range[grep(pattern = "0", x = q.range)])
	nCell <- length(x = call.list[[1]])
	for (i in 1:nrow(x = res)) {
		temp <- table(call.list[[i]])
		if ("Doublet" %in% names(x = temp) == TRUE) {
			res$pDoublet[i] <- temp[which(x = names(x = temp) == "Doublet")]
		}
		if ( "Negative" %in% names(temp) == TRUE ) {
			res$pNegative[i] <- temp[which(x = names(x = temp) == "Negative")]
		}
		res$pSinglet[i] <- sum(temp[which(x = !names(x = temp) %in% c("Doublet", "Negative"))])
	}
	res.q <- res$q
	q.ind <- grep(pattern = 'q', x = colnames(x = res))
	res <- Melt(x = res[, -q.ind])
	res[, 1] <- rep.int(x = res.q, times = length(x = unique(res[, 2])))
	colnames(x = res) <- c('q', 'variable', 'value')
	res[, 4] <- res$value/nCell
	colnames(x = res)[2:4] <- c("Subset", "nCells", "Proportion")
	extrema <- res$q[LocalMaxima(x = res$Proportion[which(x = res$Subset == "pSinglet")])]
	return(list(res = res, extrema = extrema))
}

#' @importFrom KernSmooth bkde
#' @importFrom stats approxfun quantile
ClassifyCells <- function(data, q) {
	## Generate Thresholds: Gaussian KDE with bad barcode detection, outlier trimming
	## local maxima estimation with bad barcode detection, threshold definition and adjustment
	n_cells <- nrow(x = data)
	if (n_cells == 0) {
		print('No cells passed to ClassifyCells')
		return(list(calls = character(length = n_cells), thresholds = NULL))
	}
	bc_calls <- vector(mode = "list", length = n_cells)
	n_bc_calls <- numeric(length = n_cells)
	thresholds <- list()
	for (i in 1:ncol(x = data)) {
		model <- NULL
		tryCatch(expr = {
			model <- approxfun(x = bkde(x = data[, i], kernel = "normal"))
		}, error = function(e) {
			print(paste0("Unable to fit model for ", colnames(x = data)[, i], ", for ", q, "..."))
			saveRDS(data[, i], file = paste0('./', colnames(x = data)[, i], '.fail.approxfun.rds'))
		})

		# This is changed relative to seurat
		if (is.null(x = model)) {
			print(paste0("Unable to fit model for ", rownames(x = data)[i], ", for ", q, ", skipping"))
			next
		}

		x <- seq.int(
			from = quantile(x = data[, i], probs = 0.001),
			to = quantile(x = data[, i], probs = 0.999),
			length.out = 100
		)
		extrema <- LocalMaxima(x = model(x))
		if (length(x = extrema) <= 1) {
			print(paste0("No extrema/threshold found for ", colnames(x = data)[i]))
			next
		}

		low.extremum <- min(extrema)
		high.extremum <- max(extrema)
		thresh <- (x[high.extremum] + x[low.extremum])/2
		## Account for GKDE noise by adjusting low threshold to most prominent peak
		low.extremae <- extrema[which(x = x[extrema] <= thresh)]
		new.low.extremum <- low.extremae[which.max(x = model(x)[low.extremae])]
		thresh <- quantile(x = c(x[high.extremum], x[new.low.extremum]), probs = q)
		thresholds[colnames(x = data)[i]] <- thresh

		## Find which cells are above the ith threshold
		cell_i <- which(x = data[, i] >= thresh)
		n <- length(x = cell_i)
		if (n == 0) { ## Skips to next BC if no cells belong to the ith group
			next
		}
		bc <- colnames(x = data)[i]
		if (n == 1) {
			bc_calls[[cell_i]] <- c(bc_calls[[cell_i]], bc)
			n_bc_calls[cell_i] <- n_bc_calls[cell_i] + 1
		} else {
			# have to iterate, lame
			for (cell in cell_i) {
				bc_calls[[cell]] <- c(bc_calls[[cell]], bc)
				n_bc_calls[cell] <- n_bc_calls[cell] + 1
			}
		}
	}
	calls <- character(length = n_cells)

	for (i in 1:n_cells) {
		if (n_bc_calls[i] == 0) { calls[i] <- "Negative"; next }
		if (n_bc_calls[i] > 1) { calls[i] <- "Doublet"; next }
		if (n_bc_calls[i] == 1) { calls[i] <- bc_calls[[i]] }
	}
	names(x = calls) <- rownames(x = data)
	return(list(calls = calls, thresholds = thresholds))
}

Melt <- function(x) {
	if (!is.data.frame(x = x)) {
		x <- as.data.frame(x = x)
	}
	return(data.frame(
		rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
		cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
		vals = unlist(x = x, use.names = FALSE)
	))
}