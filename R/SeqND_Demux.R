#' @include Utils.R
#' @include Normalization.R
#' @include Visualization.R

utils::globalVariables(
  names = c('relative_counts', 'barcode', 'min_average_reads'),
  package = 'cellhashR',
  add = TRUE
)

#Fits a beta distribution to a subset of barcodeMatrix
get_model <- function(barcodeMatrix, breakpoint,  dist_type) {
  barcodeMatrix <- as.vector(barcodeMatrix[barcodeMatrix > breakpoint])
  fit_of_data <- fitdistrplus::fitdist(barcodeMatrix, dist_type)
  return(fit_of_data)
}

#Scans from 1 to .01 and keeps the best fitting beta model
get_right_dist_beta <- function(barcodeMatrix, tolerance=0.1) {
  fits <- lapply(seq(1,0.01, -0.01), FUN= function(x){ 
    tryCatch({
      model <- get_model(barcodeMatrix = barcodeMatrix, breakpoint = x , dist_type = "beta")
      chi_sq <- fitdistrplus::gofstat(model)$chisq
      return(chi_sq)
    }, error = function(e){
      #print('Error fitting beta distribution. Likely not enough strongly hashed cells in this cluster')
			print(conditionMessage(e))
			traceback()
      return(NULL)
    })
  })
  #Issue: Beta distributions with few points will always fit well, so we need to hack around this limitation
  #   by finding the global chi-squared value, then masking near that value to select a more reasonable minimum
  
  #hack: replace NULL values with the max_chisq so that minChisqIndex finds where the global minimum is
  fits[sapply(fits, is.null)] <- max(unlist(fits))
  
  #find the global minimum and determine how much of the vector of fits to mask
  fitsVector <- unlist(fits)
  minChisqIndex <- which(fitsVector %in% min(fitsVector)) + tolerance*100
  
  #hack: mask the chisq values near the global minimum (according to tolerance)
  fitsVector[1:minChisqIndex] <- max(fitsVector)
  #find the index of the new local minimum
  localMin <- which(fitsVector %in% min(fitsVector))
  #transform vector index to (0,1] value for model fitting
  localMin <- (1-localMin/100)
  #Refit the best model and return
  bestModel <- get_model(barcodeMatrix = barcodeMatrix, breakpoint = localMin , dist_type = "beta")
  return(bestModel)
}

GenerateCellHashCallsSeqND <- function(barcodeMatrix, assay = "HTO", min_quantile = 0.01, min_average_reads = 10, verbose = TRUE, methodName = 'seqnd'){
	if (verbose) {
		print('Starting SeqND')
	}

	tryCatch({
		#filter barcodes for low average expression:
		sel <- rowMeans(barcodeMatrix) > min_average_reads
		if (sum(sel) == 0) {
			print(paste0('No passing barcodes after filter using min_average_reads: ', min_average_reads))
			return(NULL)
		}

		barcodeMatrix <- barcodeMatrix[sel,]
		if (verbose) {
			print(paste0('rows dropped for low counts: ', sum(!sel), ' of ', length(sel)))
		}

		seuratObj <- Seurat::CreateSeuratObject(Seurat::as.sparse(barcodeMatrix), assay = assay)
        seuratObj <- SetAssayData4Or5(seuratObj, assay = assay, theLayer = 'data', new.data = NormalizeRelative(barcodeMatrix))

		seuratObj <- SeqNDDemux(seuratObj = seuratObj, min_quantile = min_quantile, assay = assay)

		SummarizeHashingCalls(seuratObj, label = 'SeqND', columnSuffix = 'seqnd', assay = assay)

        df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = methodName, classification = seuratObj$classification.seqnd, classification.global = seuratObj$classification.global.seqnd, stringsAsFactors = FALSE)
        df <- .RestoreUnderscoreToHtoNames(df, rownames(barcodeMatrix))
        return(df)
  }, error = function(e){
		print('Error generating seqnd calls, aborting')
		print(conditionMessage(e))
		traceback()

		return(NULL)
	})
}

SeqNDDemux <- function(seuratObj, assay, min_quantile = 0.01, plotcolor =  "#00BFC4") {
  #Perform thresholding
  if (!'data' %in% SeuratObject::Layers(Seurat::GetAssay(seuratObj, assay = assay))) {
    stop('Missing data layer!')
  }

  barcodeMatrix <- GetAssayData(
		object = seuratObj,
		assay = assay,
		layer = 'data'
  )

  #loop over HTOs in matrix, perform thresholding and store cells that pass the threshold
  #return a discrete matrix, with 1 equal to a call positive for that barcode
  discrete <- GetAssayData(object = seuratObj, assay = assay)
  discrete[discrete > 0] <- 0

  for (hto in rownames(barcodeMatrix)) {
    cells <- barcodeMatrix[hto, colnames(seuratObj), drop = FALSE]

    #finding the best model produces NaNs when scanning, so we suppress them
    bestModel <- suppressWarnings(get_right_dist_beta(cells))
    boot <- fitdistrplus::bootdist(bestModel)
    threshold <- quantile(boot, probs= c(min_quantile, 1))$quantiles[[1]]
    cellNormalizedBarcodes.pivot <- tidyr::pivot_longer(as.data.frame(cells),
                                                        cols = colnames(cells),
                                                        names_to = "barcode", 
                                                        values_to = "relative_counts")
    P1 <- cellNormalizedBarcodes.pivot %>%
      tidyr::gather(barcode, relative_counts) %>% 
      ggplot(aes(x = relative_counts)) +
      geom_histogram(color = plotcolor, fill = plotcolor, binwidth = .01, position="identity", alpha = 0.5)+
      xlab("Fraction of Total Counts") +
      ylab("Density") +
      geom_vline(xintercept=threshold, linetype="dashed") +
      ggtitle("SeqND Threhshold", hto)
    
    print(P1)

		discrete[hto, colnames(seuratObj)] <- ifelse(cells > threshold, yes = 1, no = 0)
  }

	# now assign cells to HTO based on discretized values
  seuratObj <- .AssignCallsToMatrix(seuratObj, discrete, suffix = 'seqnd', assay = assay)

  return(seuratObj)
}
