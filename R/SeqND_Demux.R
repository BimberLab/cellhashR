
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
      #print(e)
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


DoSeqND <- function(barcodeMatrix, assay = "HTO", min_quantile = 0.01, min_average_reads = 10){
  #filter barcodes for low average expression (at least average 1 read per cell)
  barcodeMatrix <- barcodeMatrix[rowMeans(barcodeMatrix) > min_average_reads,]
  seuratObj <- Seurat::CreateSeuratObject(barcodeMatrix, assay = assay)
  seuratObj[[assay]]@data <- prop.table(barcodeMatrix, margin=2)
  
  #TODO remove cellhashR::: implicit call once integrated into codebase
  rownames(seuratObj[[assay]]@data) <- cellhashR:::SimplifyHtoNames(rownames(seuratObj[[assay]]@data))
  
  #TODO: This currently returns a dataframe, not a Seurat Object as suggested
  seuratObj <- SeqNDDemux(seuratObj = seuratObj, min_quantile = min_quantile)
  
  #TODO make calls parseable by SummarizeHashingCalls
  #SummarizeHashingCalls(seuratObj, label = 'SeqND', htoClassificationField = 'classification.seqnd', globalClassificationField = 'classification.global.seqnd', assay = assay)
  return(seuratObj)
}

SeqNDDemux <- function(seuratObj, min_quantile = 0.01, plotcolor =  "#00BFC4") {
  #Perform thresholding
  barcodeMatrix <- seuratObj$HTO
  #loop over HTOs in matrix, perform thresholding and store cells that pass the threshold
  cellList<- lapply(1:nrow(barcodeMatrix), FUN= function(x) {
    cells <- barcodeMatrix[x,] 
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
      ggtitle("Hashing Threhshold for Barcode", rownames(barcodeMatrix)[[x]])
    
    print(P1)
    
    calls <- cells[,cells > threshold]
    return(calls)
  })
  
  allPositiveCalls <- as.data.frame(dplyr::bind_rows(cellList))
  allPositiveCalls[!is.na(allPositiveCalls)]<-1
  allPositiveCalls[is.na(allPositiveCalls)]<-0
  rownames(allPositiveCalls) <- rownames(barcodeMatrix)
  return(allPositiveCalls)
}
