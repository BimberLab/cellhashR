#' @include Utils.R
#' @include Normalization.R
#' @include Visualization.R

utils::globalVariables(
  names = c('relative_counts', 'x', 'y', '..density..', 'Highest', 'Second'),
  package = 'cellhashR',
  add = TRUE
)

SNR <- function(barcodeData) {
  # SNR extracts the minimum data necessary from an input barcode matrix and
  # places this info in a data frame.  Extracted data are (1) the barcode with
  # the highest counts, (2) the highest count value, and (3) the 2nd highest
  # count value.  SNR stands for signal to noise ratio, which the data extracted
  # can be used to calculate.
  df <- data.frame(barcodeData, check.names = FALSE)
  top <- c()
  topbar <- c()
  topval <- c()
  second <- c()
  for (i in 1:nrow(df)) {
    top[i] <- which.max(df[i, ])
    topbar[i] <- colnames(df)[top[i]]
    topval[i] <- df[i, top[i]]
    df[i,top[i]] <- 0
    second[i] <- max(df[i, ])
  }
  outdf <- data.frame("CellID" = rownames(df), check.names = FALSE)
  outdf$Barcode <- topbar
  outdf$Highest <- topval
  outdf$Second <- second
  
  return(outdf)
}

getNegNormedData <- function(discrete, barcodeMatrix) {
  negs <- colnames(discrete[, colSums(discrete)==0])
  negdiscrete <- 1 - discrete
  is.na(negdiscrete) <- negdiscrete==0
  negvals <- negdiscrete * barcodeMatrix
  neg_normed <- NormalizeQuantile(t(as.matrix(negvals)))
  neg_normed_na <- neg_normed
  neg_normed[is.na(neg_normed)] <- 0
  return(neg_normed)
}

getPosNormedData <- function(discrete, barcodeMatrix) {
  all_pos_vals <- data.frame(t(as.matrix(discrete * barcodeMatrix)))
  is.na(all_pos_vals) <- all_pos_vals==0
  all_pos_normed <- NormalizeQuantile(all_pos_vals)
  all_pos_normed[is.na(all_pos_normed)] <- 0
  return(all_pos_normed)
}

getDiscreteFromCutoffs <- function(seuratObj, assay, cutoffs) {
  barcodeMatrix <- GetAssayData(
    object = seuratObj,
    assay = assay,
    slot = 'counts'
  )[, colnames(x = seuratObj)]
  
  #loop over HTOs in matrix, perform thresholding and store cells that pass the threshold
  #return a discrete matrix, with 1 equal to a call positive for that barcode
  discrete <- GetAssayData(object = seuratObj, assay = assay)
  discrete[discrete > 0] <- 0
  for (hto in ls(cutoffs)) {
    cells <- barcodeMatrix[hto, colnames(seuratObj), drop = FALSE]
    threshold <- cutoffs[[hto]]
    discrete[hto, colnames(seuratObj)] <- ifelse(cells > threshold, yes = 1, no = 0)
  }
  # now assign cells to HTO based on discretized values
  return(discrete)
}


PlotCutoff <- function(data, smooth, label) {
  # Function to calculate the threshold between positive and negative peaks and plot the distribution, fit, and
  # derivative peaks.
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
  nbins <- 100
  P1 <- ggplot(data.frame(Value = data), aes(x = Value)) +
    geom_histogram(aes(y = sqrt(..density..)), bins = nbins) +
    ggtitle(paste0("Histogram of ", label)) +
    xlab("Log(HTO Counts)") +
    geom_line(data = data.frame(x = smooth$x, y = sqrt(smooth$y)), mapping = aes(x = x, y = y), color = "blue", size = 1) 
  # + geom_line(data = data.frame(x = smooth$x, y = 50*deriv), mapping = aes(x = x, y = y), color = "red", size = 1)
  
  ymax <- max(graphics::hist(data, breaks = nbins, plot=FALSE)$density)
  # P1 <- P1 + ylim(c(0, ymax * 1.05))
  
  if (length(max_list) > 1) {
    y1 <- max(yvals)
    index1 <- max_list[which.max(yvals)]
    x1 <- smooth$x[index1]
    y2 <- max(yvals[yvals != y1])
    index2 <- max_list[which(yvals==y2)]
    x2 <- smooth$x[index2]
    if (x2 < x1){
      #Negative peak must be at least 1/10th the positive peak
      if (y2 < (y1/10)) {
        print(label)
        print('Negative peak was not at least 1/10th the positive peak, using max value as cutoff')
        P1 <- P1 + geom_vline(xintercept = max(data), size = 1)
        print(P1)
        return(list(max(data), P1))
      }
      cutoff_indices <- index2:index1
    } else {
      cutoff_indices <- index1:index2
    }
    cutoff_index <- which.min(smooth$y[cutoff_indices]) + min(cutoff_indices)
    cutoff <- smooth$x[cutoff_index]
  } else {
    print(label)
    print('Only one peak found, using max value as cutoff')
    P1 <- P1 + geom_vline(xintercept = max(data), size = 1)
    print(P1)
    return(list(max(data), P1))
  }
  P1 <- P1 + geom_vline(xintercept = cutoff, size = 1)
  print(P1)
  return(list(cutoff, P1))
}


getCountCutoff <- function(data, label, num_deriv_peaks, barcodeBlocklist = NULL) {
  # Function to find the threshold between positive and negative peaks of a barcode's distribution
  num_peaks <- 10000
  change <- 10
  j <- 1
  max2_list <- numeric(10)
  data <- log10(as.vector(data+1))
  data <- data[data > 0]
  # Iterate KDE over sequentially larger bandwidths until fake peaks are smoothed away.
  # Stop when there is no change in number of peaks and there are fewer peaks than the number
  #    specified in num_deriv_peaks.
  while ((change > 0) | (length(max2_list) > num_deriv_peaks)) {
    j <- j + 0.5
    smooth <- stats::density(data, adjust = j, kernel = 'gaussian',
                      bw = 'SJ', give.Rkern = FALSE)
    deriv <- numeric(length(smooth$x))
    deriv2 <- numeric(length(smooth$x))
    max_list <- c()
    max2_list <- c()
    smooth_list <- c()
    y_vals <- c()
    for (i in 2:(length(smooth$x)-1)){
      deriv[i] <- smooth$y[i+1] - smooth$y[i-1]
      deriv2[i] <- smooth$y[i+1] - 2*smooth$y[i] + smooth$y[i-1]
    }
    for (i in 2:(length(smooth$x)-1)){
      if ((deriv[i+1] < 0 ) && (deriv[i] >= 0) ) {
        max_list <- c(max_list, i)
        y_vals <- c(y_vals, smooth$y[i])
      }
      if ((deriv2[i+1] < 0 ) && (deriv2[i] >= 0) ) {
        max2_list <- c(max2_list, i)
      }
    }
    
    y_vals <- y_vals[order(y_vals, decreasing = TRUE)][1:2]
    x_vals <- c()
    for (y in y_vals) {
      i <- which(smooth$y==y)
      x_vals <- c(x_vals, smooth$x[i])
    }
    x_vals <- matrix(x_vals[order(x_vals, decreasing=FALSE)], nrow = 1, ncol = 2)
    
    change <- num_peaks - length(max_list)
    num_peaks <- length(max_list)
    # outlabel <- paste(label, change, num_peaks, length(max2_list))
    outlabel <- label
  }
  cutoff_res <- PlotCutoff(data, smooth, outlabel)
  cutoff <- cutoff_res[[1]]
  P1 <- cutoff_res[[2]]
  if ((cutoff == max(data))) {
    barcodeBlocklist <- c(barcodeBlocklist, label)
  }
  return(list(cutoff, barcodeBlocklist, x_vals, P1))
}


getBFFBarcodeBlocklist <- function(barcodeMatrix) {
  # Function used to find barcodes without sufficiently bimodal data to exclude these from further BFF calculation.
  barcodeBlocklist <- NULL
  for (hto in rownames(barcodeMatrix)) {
    cells <- barcodeMatrix[hto, colnames(barcodeMatrix), drop = FALSE]
    #BFF uses a log-scale to smooth higher counts, so we transform back once we find the threshold
    cutoffresults <- getCountCutoff(cells, hto, 4, barcodeBlocklist)
    if (cutoffresults[[1]] < 2) {
      print("Threshold may be placed too low, recalculating with more smoothing.")
      cutoffresults <- getCountCutoff(cells, hto, 2, barcodeBlocklist)
    }
    barcodeBlocklist <- cutoffresults[[2]]
  }
  return(barcodeBlocklist)
}


GenerateCellHashCallsBFF <- function(barcodeMatrix, assay = "HTO", min_average_reads = 10, verbose = TRUE, simple_threshold = FALSE, doublet_thresh = 1/4, neg_thresh = 1, pos_dist = 1, neg_dist = 1/4){
  if (verbose) {
    print('Starting BFF')
  }
  #filter barcodes for low average expression:
  sel <- rowMeans(barcodeMatrix) > min_average_reads
  barcodeMatrix <- barcodeMatrix[sel,]
  
  if (nrow(barcodeMatrix) == 0) {
    print(paste0('No passing barcodes after filter using min_average_reads: ', min_average_reads))
    return(NULL)
  }
  
  if (verbose) {
    print(paste0('rows dropped for low counts: ', sum(!sel), ' of ', length(sel)))
  }

  tryCatch({
    seuratObj <- Seurat::CreateSeuratObject(barcodeMatrix, assay = assay)
    seuratObj <- BFFDemux(seuratObj = seuratObj, assay = assay, simple_threshold = simple_threshold, doublet_thresh = doublet_thresh, neg_thresh = neg_thresh, pos_dist = pos_dist, neg_dist = neg_dist)
    if (as.logical(simple_threshold) == TRUE) {
      SummarizeHashingCalls(seuratObj, label = "bff_threshold", columnSuffix = "bff_threshold", assay = assay, doHeatmap = TRUE)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = "bff_threshold", classification = seuratObj$classification.bff_threshold, classification.global = seuratObj$classification.global.bff_threshold, stringsAsFactors = FALSE)
      return(df)
    } else {
      SummarizeHashingCalls(seuratObj, label = 'bff_quantile', columnSuffix = 'bff_quantile', assay = assay, doHeatmap = TRUE)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = 'bff_quantile', classification = seuratObj$classification.bff_quantile, classification.global = seuratObj$classification.global.bff_quantile, stringsAsFactors = FALSE)
      return(df)
    }
  }, error = function(e){
    warning('Error generating BFF calls, aborting')
    print(conditionMessage(e))
    traceback()
    return(NULL)
  })
}


BFFDemux <- function(seuratObj, assay, simple_threshold=simple_threshold, doublet_thresh=doublet_thresh, neg_thresh=neg_thresh, pos_dist=pos_dist, neg_dist=neg_dist) {
  barcodeMatrix <- GetAssayData(
    object = seuratObj,
    assay = assay,
    slot = 'counts'
  )[, colnames(x = seuratObj)]
  
  print(paste("Simple Threshold: ", simple_threshold))
  print(paste("Doublet thresh: ", doublet_thresh))
  print(paste("Neg thresh: ", neg_thresh))
  print(paste("Pos dist: ", pos_dist))
  print(paste("Neg dist: ", neg_dist))

  #loop over HTOs in matrix, perform thresholding and store cells that pass the threshold
  #return a discrete matrix, with 1 equal to a positive call for that barcode
  discrete <- GetAssayData(object = seuratObj, assay = assay)
  discrete[discrete > 0] <- 0
  cutoffs <- list()
  for (hto in rownames(barcodeMatrix)) {
    cells <- barcodeMatrix[hto, colnames(seuratObj), drop = FALSE]
    #BFF uses a log-scale to smooth higher counts, so we transform back once we find the threshold
    cutoffresults <- getCountCutoff(cells, hto, 4)
    if (cutoffresults[[1]] < 2) {
      print("Threshold may be placed too low, recalculating with more smoothing.")
      cutoffresults <- getCountCutoff(cells, hto, 2)
    }
    kernel_j <- cutoffresults[[3]]
    threshold <- 10^(cutoffresults[[1]])
    discrete[hto, colnames(seuratObj)] <- ifelse(cells > threshold, yes = 1, no = 0)
    cutoffs[[hto]] <- threshold
  }
  
  print("Thresholds:")
  for (cutoff in names(cutoffs)) {
    print(paste0(cutoff, ': ', cutoffs[[cutoff]]))
  }

  if (simple_threshold == TRUE) {
    seuratObj <- .AssignCallsToMatrix(seuratObj, discrete, suffix = 'bff_threshold', assay = assay)
    return(seuratObj)
  } else {
    #   Quantile method compares the distance between the top 2 normalized barcode counts to the
    #   distance between the nearest peak and the threshold.  The bc count for a cell considered for negative 
    #   recovery must be close to the threshold.  Likewise, the 2nd highest bc count for a cell considered for
    #   doublet recovery must be close to the threshold.

    neg_norm <- getNegNormedData(discrete, barcodeMatrix)
    pos_norm <- getPosNormedData(discrete, barcodeMatrix)
    tot_normed <- pos_norm + neg_norm
    
    normed_cutoffs <- list()
    max_list <- c()
    for (hto in colnames(tot_normed)) {
      cells <- tot_normed[rownames(tot_normed), hto, drop = FALSE]
      cutoffresults <- getCountCutoff(cells, hto, 4)
      if (cutoffresults[[1]] < 2) {
        print("Threshold may be placed too low, recalculating with more smoothing.")
        cutoffresults <- getCountCutoff(cells, hto, 2)
      }
      threshold <- (cutoffresults[[1]])
      max_list <- rbind(max_list, cutoffresults[[3]])
      normed_cutoffs[[hto]] <- threshold
    }
    
    norm_cutoff <- mean(unlist(normed_cutoffs))
    maxima <- colMeans(max_list)
    neg_mode <- maxima[1]
    pos_mode <- maxima[2]
    # d is the log-scale distance (after normalization) between the positive peak and threshold.
    d <- pos_mode - norm_cutoff
    # w is the log-scale distance (after normalization) between the negative peak and threshold.
    w <- norm_cutoff - neg_mode
    
    snr <- SNR(log10(tot_normed + 1))

    called <- c()
    
    classification <- c()
    for (i in 1:nrow(snr)) {
      
      if (snr[i, "Highest"] >= norm_cutoff) {
        if (snr[i, "Second"] <= norm_cutoff + d*doublet_thresh) {
          if (snr[i, "Highest"] - snr[i, "Second"] >= d*pos_dist) {
            called <- c(called, snr[i, "CellID"])
            classification[i] <- "Singlet"
          } else {
            classification[i] <- "Doublet"
          }
        } else {
          classification[i] <- "Doublet"
        }
      } else {
        if (snr[i, "Highest"] >= norm_cutoff - w*neg_thresh) {
          if (snr[i, "Highest"] - snr[i, "Second"] >= w*neg_dist) {
            called <- c(called, snr[i, "CellID"])
            classification[i] <- "Singlet"
          } else {
            classification[i] <- "Negative"
          }
        } else {
          classification[i] <- "Negative"
        }
      }
    }

    for (cell in called) {
      discrete[, cell] <- 0
      row_max <- which.max(tot_normed[cell,])
      discrete[row_max, cell] <- 1
    }
    
    joined <- cbind(snr, classification)
    
    boundaries1 <- data.frame(x = c(norm_cutoff, norm_cutoff+d*(doublet_thresh+pos_dist), max(joined$Highest) + 0.1),
                              y = c(norm_cutoff-d*pos_dist,  norm_cutoff+d*doublet_thresh, norm_cutoff+d*doublet_thresh))
    boundaries2 <- data.frame(x = c(norm_cutoff- w*neg_thresh, norm_cutoff- w*neg_thresh, norm_cutoff),
                              y = c(min(joined$Second) - 0.1, norm_cutoff- w*(neg_thresh+neg_dist), norm_cutoff-w*neg_dist))
    
    print(ggplot2::ggplot(joined, aes(x=Highest, y=Second, color=classification)) + 
            geom_point(cex=0.25) + geom_hline(yintercept = norm_cutoff) +
            geom_vline(xintercept = norm_cutoff) + 
            geom_line(aes(x=x, y=y), data = boundaries1, color="black", linetype="dashed") +
            geom_line(aes(x=x, y=y), data = boundaries2, color="black", linetype="dashed"))
    
    seuratObj <- .AssignCallsToMatrix(seuratObj, discrete, suffix = 'bff_quantile', assay = assay)
    return(seuratObj)
  }
}
