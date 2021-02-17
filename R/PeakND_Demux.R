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

callFile2 <- "/Users/boggy/bimberlab/dist_normed.csv"
callFile3 <- "/Users/boggy/bimberlab/neg_normed.csv"
neg_disc_outfile <- "/Users/boggy/bimberlab/neg_discrete.csv"
neg_vals_outfile <- "/Users/boggy/bimberlab/neg_vals.csv"
disc_outfile <- "/Users/boggy/bimberlab/disc.csv"

top2_1 <- function(discrete, barcodeMatrix) {
  negs <- colnames(discrete[, colSums(discrete)==0])
  filtered_negs <- filter_lowmax(t(barcodeMatrix)[negs,], 10)
  dist_normed <- NormalizeQuantile(t(barcodeMatrix))
  write.table(as.data.frame(dist_normed), file = callFile2, sep = '\t', row.names = TRUE, quote = FALSE, col.names = NA)
  top2_negs <- get_Top2(dist_normed[filtered_negs,])
  top2_negs_relnormed <- t(NormalizeRelative(as.matrix(t(top2_negs))))
  return(list(dist_normed, top2_negs_relnormed))
}

top2_2 <- function(discrete, barcodeMatrix) {
  negs <- colnames(discrete[, colSums(discrete)==0])
  negdiscrete <- 1 - discrete
  is.na(negdiscrete) <- negdiscrete==0
  write.table(as.data.frame(negdiscrete), file = neg_disc_outfile, sep = '\t', row.names = TRUE, quote = FALSE, col.names = NA)
  negvals <- negdiscrete * barcodeMatrix
  write.table(as.data.frame(negvals), file = neg_vals_outfile, sep = '\t', row.names = TRUE, quote = FALSE, col.names = NA)
  neg_normed <- NormalizeQuantile(t(negvals))
  neg_normed[is.na(neg_normed)] <- 0
  write.table(as.data.frame(neg_normed), file = callFile3, sep = '\t', row.names = TRUE, quote = FALSE, col.names = NA)
  filtered_negs <- filter_lowmax(t(barcodeMatrix)[negs,], 10)
  top2_negs <- get_Top2(neg_normed[filtered_negs,])
  top2_negs_relnormed <- t(NormalizeRelative(as.matrix(t(top2_negs))))
  return(list(neg_normed, top2_negs_relnormed))
}

top2_doublets <- function(discrete, barcodeMatrix) {
  called <- colnames(discrete[, colSums(discrete)>=1])
  called_vals <- t(discrete[,called]* barcodeMatrix[,called])
  called_vals_df <- as.data.frame(called_vals)
  is.na(called_vals_df) <- called_vals_df==0
  pos_normed <- NormalizeQuantile(called_vals_df)
  multi <- colnames(discrete[, colSums(discrete)>1])
  multi_pos_normed <- pos_normed[multi,]
  multi_pos_normed[is.na(multi_pos_normed)] <- 0
  
  top2_multi_pos_normed <- get_Top2(multi_pos_normed)
  top2_pos_relnormed <- t(NormalizeRelative(as.matrix(t(top2_multi_pos_normed))))
  return(list(pos_normed, top2_pos_relnormed))
}

#New PeakND

evaluate_objective <- function(seuratObj, cutoffs, penalty="log"){
  if (penalty=='square'){
    obj_value <- sum((colSums(as.matrix(getDiscreteFromCutoffs(seuratObj, 'HTO', cutoffs)))-1)^2)
  }
  if (penalty == "log"){
    obj_value <- sum(10^(abs(colSums(as.matrix(getDiscreteFromCutoffs(seuratObj, 'HTO', cutoffs)))-1)))
  }
  if (penalty == "magnitude"){
    obj_value <- sum(abs(colSums(as.matrix(getDiscreteFromCutoffs(seuratObj, 'HTO', cutoffs)))-1))
  }
  return(obj_value)
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
    #PeakND uses a log-scale to smooth higher counts, so we transform back once we find the threshold
    threshold <- cutoffs[[hto]]
    discrete[hto, colnames(seuratObj)] <- ifelse(cells > threshold, yes = 1, no = 0)
  }
  # now assign cells to HTO based on discretized values
  
  return(discrete)
}

gridSearch <- function(seuratObj, barcodeMatrix, hto, cutoffs_in, penalty, num = 50) {
  # returns vector of objective function values vs. threshold for an HTO
  cutoffs <- hash::copy(cutoffs_in)
  cutoff <- cutoffs[[hto]]
  data <- t(barcodeMatrix)[,hto]
  sorted <- data[order(data)]
  higher <- sorted[sorted > cutoff]
  lower <- sorted[sorted < cutoff]
  #TODO: This has edge cases that break the algorithm
  num <- min(num, length(lower)-1, length(higher)-1)
  lower_eval <- tail(lower, num + 1)
  higher_eval <- head(higher, num + 1)
  if (num < 3) {
    obj_vals <- data.frame("cutoff" = cutoff, "obj" = evaluate_objective(seuratObj, cutoffs))
    names(obj_vals) <- c("cutoff", "obj")
    return(obj_vals)
  }
  cutoff_temp <- mean(c(lower_eval[[1]], lower_eval[[2]]))
  cutoffs[[hto]] <- cutoff_temp
  obj_vals <- data.frame("cutoff" = cutoff_temp, "obj" = evaluate_objective(seuratObj, cutoffs))
  names(obj_vals) <- c("cutoff", "obj")
  #Scan cutoff values below original PeakND threshold
  
  for (i in 2:num) {
    cutoff_temp <- mean(c(lower_eval[[i]], lower_eval[[i+1]]))
    cutoffs[[hto]] <- cutoff_temp
    new_row <- data.frame("cutoff"= cutoff_temp,"obj"= evaluate_objective(seuratObj, cutoffs))
    obj_vals <- rbind(obj_vals, new_row)
  }
  cutoffs[[hto]] <- cutoff
  new_row <- data.frame("cutoff"= cutoff,"obj"= evaluate_objective(seuratObj, cutoffs))
  obj_vals <- rbind(obj_vals, new_row)
  #Scan cutoff values above original PeakND threshold
  for (i in 1:num) {
    cutoff_temp <- mean(c(higher_eval[[i]], higher_eval[[i+1]]))
    cutoffs[[hto]] <- cutoff_temp
    new_row <- data.frame("cutoff"= cutoff_temp,"obj"= evaluate_objective(seuratObj, cutoffs))
    obj_vals <- rbind(obj_vals, new_row)
  }
  return(obj_vals)
}

getCutoffsFromGridSearch <- function(seuratObj, barcodeMatrix, cutoffs, penalty ) {
  # loops over HTos to find optimal threshold values
  for (hto in rownames(barcodeMatrix)) {
    gs <- gridSearch(seuratObj, barcodeMatrix, hto, cutoffs, penalty = penalty)
    vals <- which(gs[['obj']]==min(gs[['obj']]))
    cutoffs[[hto]] <- gs[max(vals),'cutoff']
    plot(gs, main=hto)
    abline(v=gs[vals,'cutoff'])
  }
  return(cutoffs)
}

#loops over getCutoffsFromGridSearch() until the thresholds converge
loopsearch <- function(seuratObj, barcodeMatrix, cutoffs, penalty = 'log') {
  cutoffs2 <- hash::copy(cutoffs)
  cutoffs2 <- getCutoffsFromGridSearch(seuratObj, barcodeMatrix, cutoffs2)
  i <- 1
  #TODO: find best convergence criterion
  while (sum(as.numeric(hash::values(cutoffs2)) - as.numeric(hash::values(cutoffs))) > 10 ) {
    i <- i + 1
    cutoffs <- hash::copy(cutoffs2)
    cutoffs2 <- getCutoffsFromGridSearch(seuratObj, barcodeMatrix, cutoffs, penalty)
  }
  return(cutoffs2)
}

# Old PeakND


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
  }
  else {
    print('Second peak was not at least 1/10th the maximum peak, using max value as cutoff')
    P1 <- P1 + geom_vline(xintercept = max(data), size = 1)
    print(P1)
    
    return(max(data))
  }
  
  P1 <- P1 + geom_vline(xintercept = cutoff, size = 1)
  print(P1)
  return(cutoff)

}

getCountCutoff <- function(data, label, barcodeBlocklist = NULL) {
  #TODO: See about making this code less dependent on initial conditions
  num_peaks <- 10000
  change <- 10
  j <- 1
  max2_list <- numeric(10)
  data <- log10(as.vector(data))
  data <- data[data > 0]
  while ((change > 0) | (length(max2_list) > 4)) {
    j <- j + 0.5
    #TODO?: Potentially add support for different kernels/bandwidths? 
    smooth <- stats::density(data, adjust = j, kernel = 'gaussian',
                      bw = 'SJ', give.Rkern = FALSE)
    deriv <- numeric(length(smooth$x))
    deriv2 <- numeric(length(smooth$x))
    max_list <- c()
    max2_list <- c()
    for (i in 2:(length(smooth$x)-1)){
      deriv[i] <- smooth$y[i+1] - smooth$y[i-1]
      deriv2[i] <- smooth$y[i+1] - 2*smooth$y[i] + smooth$y[i-1]
    }
    for (i in 2:(length(smooth$x)-1)){
      if ((deriv[i+1] < 0 ) && (deriv[i] >= 0) ) {
        max_list <- c(max_list, i)
      }
      if ((deriv2[i+1] < 0 ) && (deriv2[i] >= 0) ) {
        max2_list <- c(max2_list, i)
      }
    }
    change <- num_peaks - length(max_list)
    num_peaks <- length(max_list)
    outlabel <- paste(label, change, num_peaks, length(max2_list))
    
  }
  cutoff <- PlotCutoff(data, smooth, outlabel)
  if (cutoff == max(data)) {
    barcodeBlocklist <- c(barcodeBlocklist, label)
  }

  return(list(cutoff, barcodeBlocklist))
}

getPeakNDBarcodeWhitelist <- function(barcodeMatrix) {
  # barcodeMatrix <- GetAssayData(
  #   object = seuratObj,
  #   assay = 'HTO',
  #   slot = 'counts'
  # )[, colnames(x = seuratObj)]
  barcodeBlocklist <- NULL
  
  for (hto in rownames(barcodeMatrix)) {
    cells <- barcodeMatrix[hto, colnames(barcodeMatrix), drop = FALSE]
    #PeakND uses a log-scale to smooth higher counts, so we transform back once we find the threshold
    cutoffresults <- getCountCutoff(cells, hto, barcodeBlocklist)
    barcodeBlocklist <- cutoffresults[[2]]
  }
  PeakNDBarcodeWhitelist <- row.names(barcodeMatrix[!row.names(barcodeMatrix) %in% barcodeBlocklist,])
  return(PeakNDBarcodeWhitelist)
}


GenerateCellHashCallsPeakND <- function(barcodeMatrix, assay = "HTO", min_average_reads = 10, verbose = TRUE, methodName = 'peaknd', optimize_cutoffs, recover, doublet_thresh, neg_thresh, rec_meth){
  if (verbose) {
    print('Starting PeakND')
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
    seuratObj <- PeakNDDemux(seuratObj = seuratObj, assay = assay, optimize_cutoffs = optimize_cutoffs, recover=recover, doublet_thresh=doublet_thresh, neg_thresh=neg_thresh, rec_meth=rec_meth)
    print(paste("Optimize: ", optimize_cutoffs))
    print(paste("Recover: ", recover))
    print(paste("Doublet thresh: ", doublet_thresh))
    print(paste("Neg thresh: ", neg_thresh))
    print(paste("Recovery method: ", rec_meth))
    if (optimize_cutoffs==TRUE && recover==TRUE && rec_meth == 1){
      SummarizeHashingCalls(seuratObj, label = 'PeakND_opt_rec', columnSuffix = 'peaknd_opt_rec', assay = assay, doHeatmap = FALSE)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = 'peaknd_opt_rec', classification = seuratObj$classification.peaknd_opt_rec, classification.global = seuratObj$classification.global.peaknd_opt_rec, stringsAsFactors = FALSE)
      return(df)
    } else if (optimize_cutoffs==TRUE && recover==TRUE && rec_meth == 2){
      SummarizeHashingCalls(seuratObj, label = 'PeakND_opt_rec2', columnSuffix = 'peaknd_opt_rec2', assay = assay, doHeatmap = FALSE)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = 'peaknd_opt_rec2', classification = seuratObj$classification.peaknd_opt_rec2, classification.global = seuratObj$classification.global.peaknd_opt_rec2, stringsAsFactors = FALSE)
      return(df)
    } else if (optimize_cutoffs==FALSE && recover==FALSE){
      SummarizeHashingCalls(seuratObj, label = 'PeakND', columnSuffix = 'peaknd', assay = assay, doHeatmap = FALSE)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = 'peaknd', classification = seuratObj$classification.peaknd, classification.global = seuratObj$classification.global.peaknd, stringsAsFactors = FALSE)
      return(df)
    } else if (optimize_cutoffs==FALSE && recover==TRUE && rec_meth == 1){
      SummarizeHashingCalls(seuratObj, label = 'PeakND_rec', columnSuffix = 'peaknd_rec', assay = assay, doHeatmap = FALSE)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = 'peaknd_rec', classification = seuratObj$classification.peaknd_rec, classification.global = seuratObj$classification.global.peaknd_rec, stringsAsFactors = FALSE)
      return(df)
    } else if (optimize_cutoffs==FALSE && recover==TRUE && rec_meth == 2){
      SummarizeHashingCalls(seuratObj, label = 'PeakND_rec2', columnSuffix = 'peaknd_rec2', assay = assay, doHeatmap = FALSE)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = 'peaknd_rec2', classification = seuratObj$classification.peaknd_rec2, classification.global = seuratObj$classification.global.peaknd_rec2, stringsAsFactors = FALSE)
      return(df)
    } else if (optimize_cutoffs==TRUE && recover==FALSE){
      SummarizeHashingCalls(seuratObj, label = 'PeakND_opt', columnSuffix = 'peaknd_opt', assay = assay, doHeatmap = FALSE)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = 'peaknd_opt', classification = seuratObj$classification.peaknd_opt, classification.global = seuratObj$classification.global.peaknd_opt, stringsAsFactors = FALSE)
      return(df)
    }
    
    
  }, error = function(e){
    warning('Error generating peaknd calls, aborting')
    print(conditionMessage(e))
    traceback()
    return(NULL)
  })
}

PeakNDDemux <- function(seuratObj, assay, optimize_cutoffs, recover, doublet_thresh, neg_thresh, rec_meth) {
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
    #PeakND uses a log-scale to smooth higher counts, so we transform back once we find the threshold
    cutoffresults <- getCountCutoff(cells, hto)
    threshold <- 10^(cutoffresults[[1]])
    discrete[hto, colnames(seuratObj)] <- ifelse(cells > threshold, yes = 1, no = 0)
    cutoffs[[hto]] <- threshold
  }
  # barcodeMatrix <- barcodeMatrix[!row.names(barcodeMatrix) %in% barcodeBlocklist,]
  write.table(as.data.frame(discrete), file = disc_outfile, sep = '\t', row.names = TRUE, quote = FALSE, col.names = NA)
  print("Unoptimized Thresholds:")
  print(cutoffs)
  #New PeakND code
  if (optimize_cutoffs == FALSE) {
    if (recover == FALSE) {
      seuratObj <- .AssignCallsToMatrix(seuratObj, as.matrix(discrete), suffix = 'peaknd', assay = assay)
      return(seuratObj)
    } else if (recover == TRUE) {
      top2_pos_relnormed <- top2_doublets(discrete, barcodeMatrix)[[2]]
      if (rec_meth==1) {
        top2_negs_relnormed <- top2_1(discrete, barcodeMatrix)[[2]]
        discrete <- rel_Threshold(top2_negs_relnormed, discrete, doublet_thresh)
        discrete <- rel_Threshold(top2_pos_relnormed, discrete, doublet_thresh)
        seuratObj <- .AssignCallsToMatrix(seuratObj, as.matrix(discrete), suffix = 'peaknd_rec', assay = assay)
        return(seuratObj)
      } else if (rec_meth==2){
        top2_negs_relnormed <- top2_2(discrete, barcodeMatrix)[[2]]
        discrete <- rel_Threshold(top2_negs_relnormed, discrete, doublet_thresh)
        discrete <- rel_Threshold(top2_pos_relnormed, discrete, doublet_thresh)
        seuratObj <- .AssignCallsToMatrix(seuratObj, as.matrix(discrete), suffix = 'peaknd_rec2', assay = assay)
        return(seuratObj)
      }
    }
  } else if (optimize_cutoffs == TRUE){
    new_thresholds <- loopsearch(seuratObj, barcodeMatrix, cutoffs, penalty = "log")
    discrete <- getDiscreteFromCutoffs(seuratObj, assay, new_thresholds)
    print("Optimized Thresholds:")
    print(new_thresholds)
    if (recover == FALSE) {
      seuratObj <- .AssignCallsToMatrix(seuratObj, as.matrix(discrete), suffix = 'peaknd_opt', assay = assay)
      return(seuratObj)
    } else if (recover == TRUE) {
      top2_pos_relnormed <- top2_doublets(discrete, barcodeMatrix)[[2]]
      if (rec_meth==1) {
        top2_negs_relnormed <- top2_1(discrete, barcodeMatrix)[[2]]
        discrete <- rel_Threshold(top2_negs_relnormed, discrete, doublet_thresh)
        discrete <- rel_Threshold(top2_pos_relnormed, discrete, doublet_thresh)
        seuratObj <- .AssignCallsToMatrix(seuratObj, as.matrix(discrete), suffix = 'peaknd_opt_rec', assay = assay)
        return(seuratObj)
      } else if (rec_meth==2){
          top2_negs_relnormed <- top2_2(discrete, barcodeMatrix)[[2]]
          discrete <- rel_Threshold(top2_negs_relnormed, discrete, doublet_thresh)
          discrete <- rel_Threshold(top2_pos_relnormed, discrete, doublet_thresh)
          seuratObj <- .AssignCallsToMatrix(seuratObj, as.matrix(discrete), suffix = 'peaknd_opt_rec2', assay = assay)
          return(seuratObj)
      }
    }
  }
}
