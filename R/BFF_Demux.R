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
  secondbar <- c()
  topval <- c()
  second <- c()
  for (i in 1:nrow(df)) {
    top[i] <- which.max(df[i, ])
    topbar[i] <- colnames(df)[top[i]]
    topval[i] <- df[i, top[i]]
    df[i,top[i]] <- -1000
    second[i] <- max(df[i, ])
    secondbar[i] <- colnames(df)[which.max(df[i, ])]
  }
  outdf <- data.frame("CellID" = rownames(df), check.names = FALSE)
  outdf$Barcode <- topbar
  outdf$Highest <- topval
  outdf$Second <- second
  outdf$Barcode2 <- secondbar
  
  return(outdf)
}

getNegNormedData <- function(discrete, barcodeMatrix) {
  if (!identical(dim(discrete), dim(barcodeMatrix))) {
    warning('getNegNormedData being passed inputs with different dimensions')
  }

  negdiscrete <- 1 - discrete
  is.na(negdiscrete) <- negdiscrete==0
  negvals <- negdiscrete * barcodeMatrix
  neg_normed <- NormalizeQuantile(t(as.matrix(negvals)))
  neg_normed_na <- neg_normed
  neg_normed[is.na(neg_normed)] <- 0

  return(neg_normed)
}

getPosNormedData <- function(discrete, barcodeMatrix) {
  if (!identical(dim(discrete), dim(barcodeMatrix))) {
    warning('getPosNormedData being passed inputs with different dimensions')
  }

  all_pos_vals <- data.frame(t(as.matrix(discrete * barcodeMatrix)), check.names = FALSE)
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
  
  plotdata <- data.frame(Value = data)
  linedata <- data.frame(x = smooth$x, y = sqrt(smooth$y))
  
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
        print('Negative peak was not at least 1/10th the positive peak, using max value as cutoff')
        cutoff <- max(data)
        return(list(cutoff, plotdata, linedata))
      }
      cutoff_indices <- index2:index1
    } else {
      cutoff_indices <- index1:index2
    }
    cutoff_index <- which.min(smooth$y[cutoff_indices]) + min(cutoff_indices)
    cutoff <- smooth$x[cutoff_index]
  } else {
    print('Only one peak found, using max value as cutoff')
    cutoff <- max(data)
    return(list(cutoff, plotdata, linedata))
  }
  return(list(
  	cutoff = cutoff,
  	plotdata = plotdata,
  	linedata = linedata
  ))
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
      # x_vals is a vector containing the count values where maxima occur in the density curve
      x_vals <- c(x_vals, smooth$x[i])
    }
    x_vals <- matrix(x_vals[order(x_vals, decreasing=FALSE)], nrow = 1, ncol = 2)
    
    change <- num_peaks - length(max_list)
    num_peaks <- length(max_list)
    outlabel <- label
  }
  cutoff_res <- PlotCutoff(data, smooth, outlabel)
  cutoff <- cutoff_res[['cutoff']]
  plotdata <- cutoff_res[['plotdata']]
  linedata <- cutoff_res[['linedata']]
  if ((cutoff == max(data))) {
    barcodeBlocklist <- c(barcodeBlocklist, label)
  }
  return(list(
		cutoff = cutoff,
		barcodeBlocklist = barcodeBlocklist,
  	x_vals = x_vals,
  	plotdata = plotdata,
  	linedata = linedata
  ))
}

generateBFFGridPlot <- function(barcodeMatrix, xlab, maintitle, universal_cutoff = NULL, smoothingThreshold = 2) {
  barcodeBlocklist = NULL
  plotdata <- NULL
  linedata <- NULL
  cutoffs <- NULL
  discrete <- barcodeMatrix
  discrete[discrete > 0] <- 0

  for (hto in rownames(barcodeMatrix)) {
    cells <- barcodeMatrix[hto, colnames(barcodeMatrix), drop = FALSE]
    cutoffresults <- getCountCutoff(cells, hto, 4, barcodeBlocklist)
    cutoffval <- cutoffresults[['cutoff']]
    x_vals <- cutoffresults[['x_vals']]

    if (!is.null(universal_cutoff)) {
      cutoffval <- universal_cutoff
    }

    #TODO: should this be tunable??
    if ((is.null(universal_cutoff)) & (cutoffval < smoothingThreshold)) {
      print(paste0("Threshold for ", hto, " may be placed too low, recalculating with more smoothing."))
      cutoffresults <- getCountCutoff(cells, paste0(hto, "*"), 2, barcodeBlocklist)
      cutoffval <- cutoffresults[['cutoff']]
      
      toAdd <- cutoffresults[['plotdata']]
      toAdd$Barcode <- hto

      if (is.null(plotdata)) {
        plotdata <- toAdd
      } else {
        plotdata <- rbind(toAdd, plotdata)
      }
      
      toAdd <- cutoffresults[['linedata']]
      toAdd$Barcode <- hto

      if (is.null(linedata)) {
        linedata <- toAdd
      } else {
        linedata <- rbind(toAdd, linedata)
      }
      
      toAdd <- data.frame(cutoff = cutoffval)
      toAdd$Barcode <- hto
      toAdd$y <- max(linedata$y) * 1.1
      
      if (is.null(cutoffs)) {
        cutoffs <- toAdd
      } else {
        cutoffs <- rbind(toAdd, cutoffs)
      }
      toAdd <- data.frame(cutoff = cutoffval)
      toAdd$Barcode <- hto
      toAdd$y <- -0.1
      cutoffs <- rbind(toAdd, cutoffs)
      
    } else {
      toAdd <- cutoffresults[['plotdata']]
      toAdd$Barcode <- hto

      if (is.null(plotdata)) {
        plotdata <- toAdd
      } else {
        plotdata <- rbind(toAdd, plotdata)
      }
      
      toAdd <- cutoffresults[['linedata']]
      toAdd$Barcode <- hto

      if (is.null(linedata)) {
        linedata <- toAdd
      } else {
        linedata <- rbind(toAdd, linedata)
      }
      
      toAdd <- data.frame(cutoff = cutoffval)
      toAdd$Barcode <- hto
      toAdd$y <- max(linedata$y) * 1.1
      
      if (is.null(cutoffs)) {
        cutoffs <- toAdd
      } else {
        cutoffs <- rbind(toAdd, cutoffs)
      }
      toAdd <- data.frame(cutoff = cutoffval)
      toAdd$Barcode <- hto
      toAdd$y <- -0.1
      cutoffs <- rbind(toAdd, cutoffs)
    }
    barcodeBlocklist <- cutoffresults[['barcodeBlocklist']]
    discrete[hto, colnames(barcodeMatrix)] <- ifelse(cells > 10^cutoffval, yes = 1, no = 0)
  }
  plotdata$Barcode <- naturalsort::naturalfactor(plotdata$Barcode)
  linedata$Barcode <- naturalsort::naturalfactor(linedata$Barcode)
  cutoffs$Barcode <- naturalsort::naturalfactor(cutoffs$Barcode)
  
  cutoffsout <- unique(cutoffs[, c("cutoff", "Barcode")])

  cutoffslist <- list()
  for (i in 1:length(cutoffsout$Barcode)) {
    cutoffslist[[as.character(cutoffsout[[i, "Barcode"]])]] <- cutoffsout[i, "cutoff"]
  }

  nbins <- 100
  maxPerPlot <- 12
  totalPages <- GetTotalPlotPages(totalValues = length(unique(plotdata$Barcode)), perPage = maxPerPlot)
  for (i in 1:totalPages) {
    print(ggplot2::ggplot(plotdata, aes(x = Value)) + geom_line(data = linedata, mapping = aes(x = x, y = y), color = "blue", size = 1) +
            egg::theme_presentation(base_size = 12) + geom_line(data=cutoffs, aes(x=cutoff, y = y), size = 1) + 
            geom_histogram(aes(y = sqrt(..density..)), size = 1, bins = nbins) + scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) + 
            labs(y = 'sqrt(Density)', x = xlab) + ggtitle(maintitle)  +
            ggforce::facet_wrap_paginate(~Barcode, scales = 'free', strip.position = 'top', nrow = min(3, length(unique(plotdata$Barcode))), labeller = labeller(.multi_line = FALSE), page = i)
    )
  }

  return(list(
		barcodeBlocklist = barcodeBlocklist,
  	cutoffslist = cutoffslist,
  	discrete = discrete,
  	x_vals = x_vals
  ))
}

GenerateCellHashCallsBFF <- function(barcodeMatrix, assay = "HTO", min_average_reads = 10, verbose = TRUE, simple_threshold = FALSE, doublet_thresh = 0.05, neg_thresh = 0.05, dist_frac = 0.1){
  if (verbose) {
    print('Starting BFF')
  }
  #filter barcodes for low average expression:
  sel <- rowMeans(barcodeMatrix) > min_average_reads
  barcodeMatrix <- barcodeMatrix[sel,,drop = FALSE]
  
  if (nrow(barcodeMatrix) == 0) {
    print(paste0('No passing barcodes after filter using min_average_reads: ', min_average_reads))
    return(NULL)
  }
  
  if (verbose) {
    print(paste0('rows dropped for low counts: ', sum(!sel), ' of ', length(sel)))
  }

  tryCatch({
    seuratObj <- Seurat::CreateSeuratObject(barcodeMatrix, assay = assay)
    seuratObj <- BFFDemux(seuratObj = seuratObj, assay = assay, simple_threshold = simple_threshold, doublet_thresh = doublet_thresh, neg_thresh = neg_thresh, dist_frac=dist_frac)
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

BFFDemux <- function(seuratObj, assay, simple_threshold=simple_threshold, doublet_thresh=doublet_thresh, neg_thresh=neg_thresh, dist_frac=dist_frac) {
  barcodeMatrix <- GetAssayData(
    object = seuratObj,
    assay = assay,
    slot = 'counts'
  )[, colnames(x = seuratObj)]
  
  print(paste("Simple Threshold: ", simple_threshold))
  print(paste("Doublet thresh: ", doublet_thresh))
  print(paste("Neg thresh: ", neg_thresh))
  print(paste("Min distance as fraction of distance between peaks: ", dist_frac))

  thresholdres <- generateBFFGridPlot(barcodeMatrix, "Log(Counts + 1)", "Raw Count Distributions with BQN Thresholds")
  
  cutofflist <- thresholdres[['cutoffslist']]
  cutoffs <- list()
  for (bar in names(cutofflist)) {
    cutoffs[[bar]] <- 10^cutofflist[[bar]]
  }
  
  discrete <- thresholdres[['discrete']]
  x_vals <- thresholdres[['x_vals']]

  
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

    #outputs from NormalizeBimodalQuantile: discrete, tot_normed, log10(tot_normed+1), barcodeBlocklist
    normedres <- NULL
    tryCatch({
      normedres <- NormalizeBimodalQuantile(barcodeMatrix)
    }, error = function(e){
      print("No valid barcodes, skipping BFF")
      print(conditionMessage(e))
      traceback()
    })
    
    if (is.null(normedres)){
      return()
    }
    
    discrete <- normedres[['discrete']]
    tot_normed <- normedres[['tot_normed']]
    lognormedcounts <- normedres[['lognormedcounts']]
    barcodeBlocklist <- normedres[['barcodeBlocklist']]
    
    #generateBFFGridPlot outputs barcodeBlocklist, cutoffslist, discrete, x_vals
    normedplotres <- generateBFFGridPlot(t(tot_normed), "Log(Counts + 1)", "Normalized Count Distributions with Fitted Threshold")
    
    x_vals <- normedplotres[['x_vals']]
    normed_cutoffs <- normedplotres[['cutoffslist']]
    max_list <- x_vals
    norm_cutoff <- mean(unlist(normed_cutoffs))

    #TODO: why result ignored?
    foo <- generateBFFGridPlot(t(tot_normed), "Log(Counts + 1)", "Normalized Count Distributions with Final Threshold", universal_cutoff = norm_cutoff)
    
    maxima <- colMeans(max_list)
    neg_mode <- maxima[1]
    pos_mode <- maxima[2]

    snr <- SNR(lognormedcounts)
    called <- c()
    
    highest_dist <- stats::density(snr$Highest, adjust = 1, kernel = 'gaussian',
                                   bw = 'SJ', give.Rkern = FALSE)
    second_dist <- stats::density(snr$Second, adjust = 1, kernel = 'gaussian',
                                  bw = 'SJ', give.Rkern = FALSE)
    neg_cutoff <- highest_dist$x[[min(which((abs(highest_dist$y - doublet_thresh*max(highest_dist$y))) < 0.01))]]
    doublet_cutoff <- second_dist$x[[max(which((abs(second_dist$y - neg_thresh*max(second_dist$y))) < 0.01))]]

    classification <- c()
    for (i in 1:nrow(snr)) {
      if (snr[i, "Highest"] >= neg_cutoff) {
        if (snr[i, "Second"] <= doublet_cutoff) {
          if (snr[i, "Highest"] - snr[i, "Second"] >= dist_frac * (pos_mode - neg_mode)) {
            called <- c(called, snr[i, "CellID"])
            classification[i] <- "Singlet"
          } else {
            classification[i] <- "Doublet"
          }
          
        } else {
          classification[i] <- "Doublet"
        }
      } else {
        classification[i] <- "Negative"
      }
    }

    for (cell in called) {
      discrete[, cell] <- 0
      row_max_name <- colnames(tot_normed)[which.max(tot_normed[cell,])]
      discrete[row_max_name, cell] <- 1
    }

    joined <- cbind(snr, classification)
    
    P1 <- ggplot2::ggplot(joined, aes(x=Highest, y=Second, color=classification)) + 
            geom_point(cex=0.25) + geom_hline(yintercept = doublet_cutoff) +
            geom_vline(xintercept = neg_cutoff) + ggtitle("BFF Droplet Classifications") +
            geom_segment(x=neg_cutoff, y=neg_cutoff - dist_frac * (pos_mode - neg_mode),
                   xend= doublet_cutoff+dist_frac * (pos_mode - neg_mode), yend=doublet_cutoff,
                   linetype="dashed", color="black") +
            egg::theme_presentation(base_size = 10) +
            theme(legend.position = c(0.1, 0.65), legend.text=element_text(size=10)) +
            guides(colour = guide_legend(override.aes = list(size=3)))
    P2 <- ggExtra::ggMarginal(P1, size=4)
    grid::grid.newpage()
    grid::grid.draw(P2)

    seuratObj <- .AssignCallsToMatrix(seuratObj, discrete, suffix = 'bff_quantile', assay = assay)
    return(seuratObj)
  }
}
