#' @include Utils.R
#' @include Normalization.R
#' @include Visualization.R

utils::globalVariables(
  names = c('relative_counts', 'x', 'y', '..density..'),
  package = 'cellhashR',
  add = TRUE
)

getCountCutoff <- function(data, label) {
  #TODO: See about making this code less dependent on initial conditions
  num_peaks <- 10000
  change <- 10
  j <- 0
  data <- log10(as.vector(data))
  while ((change > 0) | (num_peaks > 2)) {
    j <- j + 0.5
    #TODO?: Potentially add support for different kernels/bandwidths? 
    smooth <- stats::density(data, adjust = j, kernel = 'gaussian',
                      bw = 'SJ', give.Rkern = FALSE)
    deriv <- numeric(length(smooth$x))
    max_list <- c()
    for (i in 2:(length(smooth$x)-1)){
      deriv[i] <- smooth$y[i+1] - smooth$y[i-1]
    }
    for (i in 2:(length(smooth$x)-1)){
      if ((deriv[i+1] < 0 ) && (deriv[i] >= 0) ) {
        max_list <- c(max_list, i)
      }
    }
    change <- num_peaks - length(max_list)
    num_peaks <- length(max_list)
  }
  deriv <- numeric(length(smooth$x))
  max_list <- c()
  for (i in 2:(length(smooth$x)-1)){
    deriv[i] <- smooth$y[i+1] - smooth$y[i-1]
  }
  for (i in 2:(length(smooth$x)-1)){
    if ((deriv[i+1] < 0 ) && (deriv[i] >= 0) ) {
      max_list <- c(max_list, i)
    }
  }
  yvals <- sapply(max_list, FUN = function(x) {
    return(smooth$y[x])
  })
  y1 <- max(yvals)
  index1 <- max_list[which.max(yvals)]
  x1 <- smooth$x[index1]
  y2 <- max(yvals[yvals != y1])
  index2 <- max_list[which(yvals==y2)]
  x2 <- smooth$x[index2]

  nbins <- 100
  P1 <- ggplot(data.frame(Value = data), aes(x = Value)) +
    #geom_histogram(aes(y = stat(count / sum(count))), bins = nbins) +
    geom_histogram(aes(y = ..density..), bins = nbins) +
    ggtitle(paste0("Histogram of ", label)) +
    xlab("log2(HTO Counts)") +
    geom_line(data = data.frame(x = smooth$x, y = smooth$y), mapping = aes(x = x, y = y), color = "blue", size = 1) +
    geom_line(data = data.frame(x = smooth$x, y = 50*deriv), mapping = aes(x = x, y = y), color = "red", size = 1)

  ymax <- max(graphics::hist(data, breaks = nbins, plot=FALSE)$density)
  P1 <- P1 + ylim(c(0, ymax * 1.05))

  if (x2 < x1) {
    #2nd peak must be at least 1/10th the maximum peak
    if (y2 < (y1/10)) {
      print('Second peak was not at least 1/10th the maximum peak, using max value as cutoff')
      P1 <- P1 + geom_vline(xintercept = max(data), size = 1)
      print(P1)

      return(max(data))
    }

    cutoff_indices <- index2:index1
  }
  else {
    cutoff_indices <- index1:index2
  }
  cutoff_index <- which.min(smooth$y[cutoff_indices]) + min(cutoff_indices)
  cutoff <- smooth$x[cutoff_index]

  P1 <- P1 + geom_vline(xintercept = cutoff, size = 1)

  print(P1)

  return(cutoff)
}


GenerateCellHashCallsPeakND <- function(barcodeMatrix, assay = "HTO", verbose = TRUE, methodName = 'peaknd'){
  if (verbose) {
    print('Starting PeakND')
  }
  
  tryCatch({
    seuratObj <- Seurat::CreateSeuratObject(barcodeMatrix, assay = assay)
    seuratObj <- PeakNDDemux(seuratObj = seuratObj, assay = assay)

    SummarizeHashingCalls(seuratObj, label = 'PeakND', columnSuffix = 'peaknd', assay = assay, doHeatmap = FALSE)
    df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = methodName, classification = seuratObj$classification.peaknd, classification.global = seuratObj$classification.global.peaknd, stringsAsFactors = FALSE)
    return(df)
  }, error = function(e){
    warning('Error generating peaknd calls, aborting', e)
    return(NULL)
  })
}

PeakNDDemux <- function(seuratObj, assay, plotcolor =  "#00BFC4") {
  barcodeMatrix <- GetAssayData(
    object = seuratObj,
    assay = assay,
    slot = 'counts'
  )[, colnames(x = seuratObj)]

  #loop over HTOs in matrix, perform thresholding and store cells that pass the threshold
  #return a discrete matrix, with 1 equal to a call positive for that barcode
  discrete <- GetAssayData(object = seuratObj, assay = assay)
  discrete[discrete > 0] <- 0
  for (hto in rownames(barcodeMatrix)) {
    cells <- barcodeMatrix[hto, colnames(seuratObj), drop = FALSE]
    #PeakND uses a log-scale to smooth higher counts, so we transform back once we find the threshold
    threshold <- 10^getCountCutoff(cells, hto)
    discrete[hto, colnames(seuratObj)] <- ifelse(cells > threshold, yes = 1, no = 0)
  }
  # now assign cells to HTO based on discretized values
  seuratObj <- .AssignCallsToMatrix(seuratObj, as.matrix(discrete), suffix = 'peaknd', assay = assay)
  return(seuratObj)
}
