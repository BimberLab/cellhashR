#' @include Utils.R
#' @include Normalization.R
#' @include Visualization.R

utils::globalVariables(
  names = c('relative_counts', 'x', 'y', '..density..'),
  package = 'cellhashR',
  add = TRUE
)


get_Top2 <- function(mat) {
  outmat <- mat
  for (barcode in rownames(outmat)) {
    outmat[barcode,][which(outmat[barcode,] < TopN(outmat[barcode,], N = 2)[[1]])] <- 0
  }
  return(outmat)
}


rel_Threshold <- function(top2, discrete, threshold) {
  outdiscrete <- discrete
  for (barcode in rownames(top2)) {
    if (max(top2[barcode,]) >= threshold) {
      outdiscrete[, barcode] <- 0
      outdiscrete[which.max(top2[barcode,])[[1]], barcode] <- 1
    }
  }
  return(outdiscrete)
}


filter_lowmax <- function(mat, minval) {
  barcodes <- c()
  for (barcode in rownames(mat)) {
    if (max(mat[barcode,]) > minval) {
      barcodes <- c(barcodes, barcode)
    }
  }
  return(barcodes)
}


top2_negcounts <- function(discrete, barcodeMatrix) {
  negs <- colnames(discrete[, colSums(discrete)==0])
  negdiscrete <- 1 - discrete
  is.na(negdiscrete) <- negdiscrete==0
  negvals <- negdiscrete * barcodeMatrix
  neg_normed <- NormalizeQuantile(t(as.matrix(negvals)))
  neg_normed_na <- neg_normed
  neg_normed[is.na(neg_normed)] <- 0
  filtered_negs <- filter_lowmax(t(barcodeMatrix)[negs,], 10)
  top2_negs <- get_Top2(neg_normed[filtered_negs,])
  top2_negs_relnormed <- t(NormalizeRelative(as.matrix(t(top2_negs))))
  top2_negs <- get_Top2(neg_normed[negs,])
  return(list(neg_normed, top2_negs_relnormed, top2_negs))
}


top2_doublets <- function(discrete, barcodeMatrix) {
  called <- colnames(discrete[, colSums(discrete)>=1])
  called_vals <- t(as.matrix(discrete[,called]* barcodeMatrix[,called]))
  called_vals_df <- as.data.frame(called_vals)
  is.na(called_vals_df) <- called_vals_df==0
  pos_normed <- NormalizeQuantile(called_vals_df)
  
  all_pos_vals <- data.frame(t(as.matrix(discrete * barcodeMatrix)))
  is.na(all_pos_vals) <- all_pos_vals==0
  all_pos_normed <- NormalizeQuantile(all_pos_vals)
  all_pos_normed[is.na(all_pos_normed)] <- 0
  
  multi <- colnames(discrete[, colSums(discrete)>=1])
  multi_pos_normed <- pos_normed[multi,]
  multi_pos_normed[is.na(multi_pos_normed)] <- 0
  
  top2_multi_pos_normed <- get_Top2(multi_pos_normed)
  top2_pos_relnormed <- t(NormalizeRelative(as.matrix(t(top2_multi_pos_normed))))
  return(list(pos_normed, top2_pos_relnormed, top2_multi_pos_normed, all_pos_normed))
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
        return(max(data))
      }
      cutoff_indices <- index2:index1
    }
    else {
      cutoff_indices <- index1:index2
    }
    cutoff_index <- which.min(smooth$y[cutoff_indices]) + min(cutoff_indices)
    cutoff <- smooth$x[cutoff_index]
  }
  else {
    print(label)
    print('Only one peak found, using max value as cutoff')
    P1 <- P1 + geom_vline(xintercept = max(data), size = 1)
    print(P1)
    return(max(data))
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
    doublet_res <- top2_doublets(discrete, barcodeMatrix)
    top2_pos_relnormed <- doublet_res[[2]]
    neg_res <- top2_negcounts(discrete, barcodeMatrix)
    top2_negs <- neg_res[[3]]
    top2_negs_log <- log10(top2_negs + 1)
    top2_multi <- doublet_res[[3]]
    top2_multi_log <- log10(top2_multi + 1)
    neg_norm <- neg_res[[1]]
    pos_norm <- doublet_res[[4]]
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
    
    norm_cutoff <- mean(as.numeric(normed_cutoffs))
    maxima <- colMeans(max_list)
    neg_mode <- maxima[1]
    pos_mode <- maxima[2]
    # d is the log-scale distance (after normalization) between the positive peak and threshold.
    d <- pos_mode - norm_cutoff
    # w is the log-scale distance (after normalization) between the negative peak and threshold.
    w <- norm_cutoff - neg_mode
    
    called <- c()
    # this loop classifies negatives vs. singlets (because the highest count is less than the threshold found by KDE)
    for (cell in rownames(top2_negs_log)) {
      max_val <- max(top2_negs_log[cell,])
      min_val <- min(top2_negs_log[cell,])
      # if max value >= threshold - alpha (w*neg_thresh)
      if (max_val >= norm_cutoff - w*neg_thresh) {
        # and if diff between top 2 counts >= phi
        if (max_val - min_val >= w*neg_dist) {
          # add this droplet to singlet group
          called <- c(called, cell)
        }
      }
    }
    
    called_multi <- c()
    # this loop classifies doublets vs. singlets (because the highest count is greater than the threshold found by KDE)
    for (cell in colnames(discrete[, colSums(discrete)>=1])) {
      # we are comparing counts on the log-scale
      top2 <- get_Top2(log10(tot_normed[cell,]+1))
      max_val <- max(top2)
      min_val <- MaxN(top2)[[1]]
      # if 2nd highest count is <= the threshold + beta (d*doublet_thresh)
      if (min_val <= norm_cutoff + d*doublet_thresh) {
        # and if the diff between top 2 counts is >= theta
        if ((max_val - min_val) >= d*pos_dist) {
          # add this droplet to the singlet group
          called_multi <- c(called_multi, cell)
        }
      }
    }
    
    for (cell in c(called, called_multi)) {
      discrete[, cell] <- 0
      row_max <- which.max(tot_normed[cell,])
      discrete[row_max, cell] <- 1
    }
    
    seuratObj <- .AssignCallsToMatrix(seuratObj, discrete, suffix = 'bff_quantile', assay = assay)
    return(seuratObj)
  }
}
