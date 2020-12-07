#' @title Perform HTO classification using Seurat's implementation of MULTI-seq classification
#'
#' @description Perform HTO classification using Seurat's implementation of MULTI-seq classification
#' @param seuratObj, A Seurat object.
#' @return A modified Seurat object.
DoMULTIseqDemux <- function(seuratObj, assay = 'HTO', autoThresh = TRUE, quantile = NULL, maxiter = 20, qrange = seq(from = 0.2, to = 0.95, by = 0.05)) {
	## Normalize Data: Log2 Transform, mean-center
	counts <- GetAssayData(seuratObj, assay = assay, slot = 'counts')
	log2Scaled <- as.data.frame(log2(Matrix::t(counts)))
	for (i in 1:ncol(log2Scaled)) {
		ind <- which(is.finite(log2Scaled[,i]) == FALSE)
		log2Scaled[ind,i] <- 0
		log2Scaled[,i] <- log2Scaled[,i]-mean(log2Scaled[,i])
	}
	seuratObjMS <- CreateSeuratObject(counts, assay = 'MultiSeq')
	seuratObjMS[['MultiSeq']]@data <- Matrix::t(as.matrix(log2Scaled))

	PlotHtoCountData(seuratObjMS, label = 'MultiSeq', assay = 'MultiSeq')

	seuratObjMS <- MULTIseqDemux2(seuratObjMS, assay = "MultiSeq", quantile = quantile, verbose = TRUE, autoThresh = autoThresh, maxiter = maxiter, qrange = qrange)

	seuratObjMS$MULTI_classification.global <- as.character(seuratObjMS$MULTI_ID)
	seuratObjMS$MULTI_classification.global[!(seuratObjMS$MULTI_ID %in% c('Negative', 'Doublet'))] <- 'Singlet'
	seuratObjMS$MULTI_classification.global <- as.factor(seuratObjMS$MULTI_classification.global)

	seuratObjMS$MULTI_ID <- naturalsort::naturalfactor(as.character(seuratObjMS$MULTI_ID))

	HtoSummary(seuratObjMS, label = 'MULTI-SEQ', htoClassificationField = 'MULTI_ID', globalClassificationField = 'MULTI_classification.global', assay = 'MultiSeq')

	seuratObj$MULTI_ID <- as.character(seuratObjMS$MULTI_ID)
	seuratObj$MULTI_classification.global <- seuratObjMS$MULTI_classification.global

	return(seuratObj)
}

#' @rawNamespace import(Matrix, except = c('tail', 'head'))
MULTIseqDemux2 <- function(
object,
assay = "HTO",
quantile = 0.7,
autoThresh = FALSE,
maxiter = 5,
qrange = seq(from = 0.1, to = 0.9, by = 0.05),
verbose = TRUE
) {
	if (is.na(assay) || is.null(assay)) {
		assay <- DefaultAssay(object = object)
	}
	multi_data_norm <- t(x = GetAssayData(
	object = object,
	slot = "data",
	assay = assay
	))
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
				bar.table_sweep.list[[n]] <- ClassifyCells(data = multi_data_norm, q = q)
				names(x = bar.table_sweep.list)[n] <- paste0("q=" , q)
			}

			# Determine which q values results in the highest pSinglet
			threshold.results1 <- Seurat:::FindThresh(call.list = bar.table_sweep.list)
			res_round <- threshold.results1$res
			res.use <- res_round[res_round$Subset == "pSinglet", ]
			q.use <- res.use[which.max(res.use$Proportion),"q"]

			round.calls <- ClassifyCells(data = multi_data_norm, q = q.use)
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
		demux_result <- ClassifyCells(data = multi_data_norm, q = quantile)
	}
	demux_result <- demux_result[rownames(x = object[[]])]
	object[['MULTI_ID']] <- factor(x = demux_result)
	Idents(object = object) <- "MULTI_ID"
	bcs <- colnames(x = multi_data_norm)
	bc.max <- bcs[apply(X = multi_data_norm, MARGIN = 1, FUN = which.max)]
	bc.second <- bcs[unlist(x = apply(
	X = multi_data_norm,
	MARGIN = 1,
	FUN = function(x) {
		return(which(x == Seurat:::MaxN(x)))
	}
	))]
	doublet.names <- unlist(x = lapply(
	X = 1:length(x = bc.max),
	FUN = function(x) {
		return(paste(sort(x = c(bc.max[x], bc.second[x])), collapse =  "_"))
	}
	))
	doublet.id <- which(x = demux_result == "Doublet")
	MULTI_classification <- as.character(object$MULTI_ID)
	MULTI_classification[doublet.id] <- doublet.names[doublet.id]
	object$MULTI_classification <- factor(x = MULTI_classification)
	return(object)
}

#' @importFrom KernSmooth bkde
#' @importFrom stats approxfun quantile
ClassifyCells <- function(data, q) {
	## Generate Thresholds: Gaussian KDE with bad barcode detection, outlier trimming
	## local maxima estimation with bad barcode detection, threshold definition and adjustment
	n_cells <- nrow(x = data)
	bc_calls <- vector(mode = "list", length = n_cells)
	n_bc_calls <- numeric(length = n_cells)
	for (i in 1:ncol(x = data)) {
		model <- NULL
		tryCatch(expr = {
			model <- approxfun(x = bkde(x = data[, i], kernel = "normal"))
		}, error = function(e) {
			print(paste0("Unable to fit model for ", rownames(x = data)[i], ", for ", q, "..."))
			saveRDS(data[, i], file = paste0('./', rownames(x = data)[i], '.fail.approxfun.rds'))
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
		extrema <- Seurat:::LocalMaxima(x = model(x))
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
		if (is.na(n_bc_calls[i]) || is.null(n_bc_calls[i]) || !is.finite(n_bc_calls[i])){
			print(head(n_bc_calls))
			stop(paste0('bad value: ', n_bc_calls[i], ' ', i))
		}

		if (n_bc_calls[i] == 0) { calls[i] <- "Negative"; next }
		if (n_bc_calls[i] > 1) { calls[i] <- "Doublet"; next }
		if (n_bc_calls[i] == 1) { calls[i] <- bc_calls[[i]] }
	}
	names(x = calls) <- rownames(x = data)
	return(calls)
}