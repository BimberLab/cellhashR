#' @include Utils.R
#' @include Normalization.R
#' @include Visualization.R

utils::globalVariables(
  names = c('relative_counts', 'x', 'y', '..density..'),
  # package = 'cellhashR',
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
  
  multi <- colnames(discrete[, colSums(discrete)>1])
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
    #BFF uses a log-scale to smooth higher counts, so we transform back once we find the threshold
    threshold <- cutoffs[[hto]]
    discrete[hto, colnames(seuratObj)] <- ifelse(cells > threshold, yes = 1, no = 0)
  }
  # now assign cells to HTO based on discretized values
  return(discrete)
}


PlotCutoff <- function(data, smooth, label) {
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
    #geom_histogram(aes(y = stat(count / sum(count))), bins = nbins) +
    geom_histogram(aes(y = ..density..), bins = nbins) +
    ggtitle(paste0("Histogram of ", label)) +
    xlab("log2(HTO Counts)") +
    geom_line(data = data.frame(x = smooth$x, y = smooth$y), mapping = aes(x = x, y = y), color = "blue", size = 1) +
    geom_line(data = data.frame(x = smooth$x, y = 50*deriv), mapping = aes(x = x, y = y), color = "red", size = 1)
  
  ymax <- max(graphics::hist(data, breaks = nbins, plot=FALSE)$density)
  P1 <- P1 + ylim(c(0, ymax * 1.05))
  
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
  return(cutoff)
}


getCountCutoff <- function(data, label, num_deriv_peaks, barcodeBlocklist = NULL) {
  num_peaks <- 10000
  change <- 10
  j <- 1
  max2_list <- numeric(10)
  data <- log10(as.vector(data+1))
  data <- data[data > 0]
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
    outlabel <- paste(label, change, num_peaks, length(max2_list))
  }
  cutoff <- PlotCutoff(data, smooth, outlabel)
  if ((cutoff == max(data))) {
    barcodeBlocklist <- c(barcodeBlocklist, label)
  }
  return(list(cutoff, barcodeBlocklist, x_vals))
}


getBFFBarcodeBlocklist <- function(barcodeMatrix) {
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


GenerateCellHashCallsBFF <- function(barcodeMatrix, assay = "HTO", min_average_reads = 10, verbose = TRUE, methodName = 'bff', recover, doublet_thresh, neg_thresh, rec_meth){
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
    seuratObj <- BFFDemux(seuratObj = seuratObj, assay = assay, recover=recover, doublet_thresh=doublet_thresh, neg_thresh=neg_thresh, rec_meth=rec_meth)
    print(paste("Recover: ", recover))
    print(paste("Doublet thresh: ", doublet_thresh))
    print(paste("Neg thresh: ", neg_thresh))
    print(paste("Recovery method: ", rec_meth))
    
    if (recover==FALSE){
      SummarizeHashingCalls(seuratObj, label = 'bff', columnSuffix = 'bff', assay = assay, doHeatmap = FALSE)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = 'bff', classification = seuratObj$classification.bff, classification.global = seuratObj$classification.global.bff, stringsAsFactors = FALSE)
      return(df)
    } else if (recover==TRUE && rec_meth == 2){
      SummarizeHashingCalls(seuratObj, label = 'bff_rec2', columnSuffix = 'bff_rec2', assay = assay, doHeatmap = FALSE)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = 'bff_rec2', classification = seuratObj$classification.bff_rec2, classification.global = seuratObj$classification.global.bff_rec2, stringsAsFactors = FALSE)
      return(df)
    } else if (recover==TRUE && rec_meth == 3){
      SummarizeHashingCalls(seuratObj, label = 'bff_dist', columnSuffix = 'bff_dist', assay = assay, doHeatmap = FALSE)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = 'bff_dist', classification = seuratObj$classification.bff_dist, classification.global = seuratObj$classification.global.bff_dist, stringsAsFactors = FALSE)
      return(df)
    }
  }, error = function(e){
    warning('Error generating BFF calls, aborting')
    print(conditionMessage(e))
    traceback()
    return(NULL)
  })
}


BFFDemux <- function(seuratObj, assay, recover, doublet_thresh, neg_thresh, rec_meth) {
  barcodeMatrix <- GetAssayData(
    object = seuratObj,
    assay = assay,
    slot = 'counts'
  )[, colnames(x = seuratObj)]

  #loop over HTOs in matrix, perform thresholding and store cells that pass the threshold
  #return a discrete matrix, with 1 equal to a call positive for that barcode
  discrete <- GetAssayData(object = seuratObj, assay = assay)
  discrete[discrete > 0] <- 0
  cutoffs <- hash::hash()
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
  print(cutoffs)
  if (recover == FALSE) {
    seuratObj <- .AssignCallsToMatrix(seuratObj, as.matrix(discrete), suffix = 'bff', assay = assay)
    return(seuratObj)
  } else if (recover == TRUE) {
    doublet_res <- top2_doublets(discrete, barcodeMatrix)
    top2_pos_relnormed <- doublet_res[[2]]
    if (rec_meth==2){
      top2_negs_relnormed <- top2_negcounts(discrete, barcodeMatrix)[[2]]
      discrete <- rel_Threshold(top2_negs_relnormed, discrete, doublet_thresh)
      discrete <- rel_Threshold(top2_pos_relnormed, discrete, doublet_thresh)
      seuratObj <- .AssignCallsToMatrix(seuratObj, as.matrix(discrete), suffix = 'bff_rec2', assay = assay)
      return(seuratObj)
    } else if (rec_meth==3) {
      neg_res <- top2_negcounts(discrete, barcodeMatrix)
      top2_negs <- neg_res[[3]]
      top2_negs_log <- log10(top2_negs + 1)
      top2_multi <- doublet_res[[3]]
      top2_multi_log <- log10(top2_multi + 1)
      neg_norm <- neg_res[[1]]
      pos_norm <- doublet_res[[4]]
      tot_normed <- pos_norm + neg_norm
      
      normed_cutoffs <- hash::hash()
      max_list <- c()
      for (hto in colnames(tot_normed)) {
        cells <- tot_normed[rownames(tot_normed), hto, drop = FALSE]
        #BFF uses a log-scale to smooth higher counts, so we transform back once we find the threshold
        cutoffresults <- getCountCutoff(cells, hto, 4)
        if (cutoffresults[[1]] < 2) {
          print("Threshold may be placed too low, recalculating with more smoothing.")
          cutoffresults <- getCountCutoff(cells, hto, 2)
        }
        threshold <- (cutoffresults[[1]])
        max_list <- rbind(max_list, cutoffresults[[3]])
        normed_cutoffs[[hto]] <- threshold
      }
      
      norm_cutoff <- mean(as.numeric(hash::values(normed_cutoffs)))
      maxima <- colMeans(max_list)
      neg_mode <- maxima[1]
      pos_mode <- maxima[2]
      d <- pos_mode - norm_cutoff
      w <- norm_cutoff - neg_mode
      
      called <- c()
      
      for (cell in rownames(top2_negs_log)) {
        max_val <- max(top2_negs_log[cell,])
        min_val <- min(top2_negs_log[cell,])
        if (max_val >= norm_cutoff - w/2) {
          if (max_val - min_val >= w/4) {
            called <- c(called, cell)
          }
        }
      }
      
      called_multi <- c()
      
      for (cell in rownames(top2_multi_log)) {
        max_val <- max(top2_multi_log[cell,])
        min_val <- MaxN(top2_multi_log[cell,])[[1]]
        if (min_val <= norm_cutoff + 2*d/3) {
          if ((max_val - min_val) >= d/3) {
            called_multi <- c(called_multi, cell)
          }
        }
      }
      for (cell in c(called, called_multi)) {
        discrete[, cell] <- 0
        row_max <- which.max(tot_normed[cell,])
        discrete[row_max, cell] <- 1
      }
      seuratObj <- .AssignCallsToMatrix(seuratObj, as.matrix(discrete), suffix = 'bff_dist', assay = assay)
      return(seuratObj)
    }
  }
}
