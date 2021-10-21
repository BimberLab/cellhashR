#' @include Utils.R
#' @include Normalization.R
#' @include Visualization.R

utils::globalVariables(
  names = c('relative_counts', 'x', 'y', '..density..', 'Highest', 'Second'),
  package = 'cellhashR',
  add = TRUE
)


#' @title Perform a parameter scan with BFF parameters
#'
#' @param lognormedcounts Count matrix after BQN and log transformation
#' @description Prints plots that demonstrate the impact of BFF parameters.
#' Parameters alpha_c, beta_c, and delta_c can take values between 0 and 1.
#' Reasonable parameter values are less than 0.25.  The scan is over parameter
#' values (0.05, 0.1, 0.15, 0.2, 0.25, 0.5) for each of the parameters.
#' @export
ParameterScan <- function(lognormedcounts) {
  snr <- SNR(lognormedcounts)
  tot_normed <- normedres[['tot_normed']]

  normedplotres <- generateBFFGridPlot(t(tot_normed), "Log(Counts + 1)", "Normalized Count Distributions with Fitted Threshold")
  normed_cutoffs <- normedplotres[['cutoffslist']]
  univ_thresh <- mean(unlist(normed_cutoffs))
  neg_mode <- mean(normedplotres[['neglist']]$neg_mode)
  pos_mode <- mean(normedplotres[['poslist']]$pos_mode)
  
  highest_dist <- stats::density(snr$Highest, adjust = 1, kernel = 'gaussian',
                                 bw = 'SJ', give.Rkern = FALSE)
  second_dist <- stats::density(snr$Second, adjust = 1, kernel = 'gaussian',
                                bw = 'SJ', give.Rkern = FALSE)
  #loop over parameter values
  neg_thresh_list <- c(NA, 0.5, 0.25,0.2,0.15,0.1,0.05)
  doublet_thresh_list <- c(NA, 0.5, 0.25,0.2,0.15,0.1,0.05)
  dist_frac_list <- c(0.5, 0.25,0.2,0.15,0.1,0.05)
  neg_counts <- c()
  neg_cutoffs <- c()
  for (neg_thresh in neg_thresh_list){
    vals <- c()
    if (!is.na(neg_thresh)) {
      for (i in 2:length(highest_dist$y)) {
        if (highest_dist$y[[i-1]] <= neg_thresh*max(highest_dist$y) & highest_dist$y[[i]] > neg_thresh*max(highest_dist$y)) {
          vals <- c(vals, i)
        }
      }
      if (length(vals) == 0) {
        print('Cannot find negative threshold.  Exiting BFF')
        return(NULL)
      }
      val <- min(vals)
      neg_cutoff <- highest_dist$x[[val]]
    } else {
      neg_cutoff <- univ_thresh
    }
    
    neg_cutoffs <- c(neg_cutoffs, neg_cutoff)
    neg_mat <- snr[snr$Highest < neg_cutoff,]
    neg_counts <- c(neg_counts, dim(neg_mat)[1])
  }
  neg_df <- data.frame(alpha_c = c(neg_thresh_list), alpha = neg_cutoffs, Negatives = neg_counts)

  Pneg <- ggplot2::ggplot(neg_df[-1,], aes(alpha_c, alpha)) +
    geom_tile(aes(fill = Negatives), colour = "black") +
    geom_hline(yintercept = neg_df[1,"alpha"]) +
    geom_text(aes(0.25,neg_df[1,"alpha"],label = neg_df[1,"Negatives"], vjust = -1)) +
    geom_text(aes(label=Negatives)) +
    ggthemes::theme_clean(base_size = 12) +
    scale_x_continuous(breaks = unique(neg_df$alpha_c)) +
    scale_y_continuous(breaks = unique(neg_df$alpha), labels = scales::number_format(accuracy = 0.001)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(neg_df$Negatives)) +
    ggtitle("Alpha vs. Alpha_c and resulting negative classifications", subtitle= "Horizontal line shows BQN threshold") +
    theme(
      plot.background = element_blank(),
      panel.grid.major.x = element_line(colour = "gray", linetype = "dotted")
    )
  
  print(Pneg)
  
  doublet_counts <- c()
  doublet_cutoffs <- c()
  
  for (doublet_thresh in doublet_thresh_list) {
    vals <- c()
    if (!is.na(doublet_thresh)) {
      for (i in 2:length(second_dist$y)) {
        if (second_dist$y[[i-1]] > doublet_thresh*max(second_dist$y) & second_dist$y[[i]] <= doublet_thresh*max(second_dist$y)) {
          vals <- c(vals, i)
        }
      }
      if (length(vals) == 0) {
        print('Cannot find doublet threshold.  Exiting BFF')
        return(NULL)
      }
      val <- max(vals)
      doublet_cutoff <- second_dist$x[[val]]
    } else {
      doublet_cutoff <- univ_thresh
    }
    
    doublet_cutoffs <- c(doublet_cutoffs, doublet_cutoff)
    doublet_mat <- snr[snr$Second > doublet_cutoff,]
    doublet_counts <- c(doublet_counts, dim(doublet_mat)[1])
  }
  doublet_df <- data.frame(beta_c = doublet_thresh_list, beta = doublet_cutoffs, Doublets = doublet_counts)

  Pdub <- ggplot2::ggplot(doublet_df[-1,], aes(beta_c, beta)) +
    geom_tile(aes(fill = Doublets), colour = "black") +
    geom_hline(yintercept = doublet_df[1,"beta"]) +
    geom_text(aes(0.25,doublet_df[1,"beta"],label = doublet_df[1,"Doublets"], vjust = 2)) +
    geom_text(aes(label=Doublets)) +
    ggthemes::theme_clean(base_size = 12) +
    scale_x_continuous(breaks = unique(doublet_df$beta_c)) +
    scale_y_continuous(breaks = unique(doublet_df$beta), labels = scales::number_format(accuracy = 0.001)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(doublet_df$Doublets)) +
    ggtitle("Beta vs. Beta_c and resulting doublet classifications", subtitle= "Horizontal line shows BQN threshold") +
    theme(
      plot.background = element_blank(),
      panel.grid.major.x = element_line(colour = "gray", linetype = "dotted")
    )
  
  print(Pdub)
  
  unidentified_counts <- c()
  distances <- c()
  for (dist_frac in dist_frac_list) {
    distance <- dist_frac * (pos_mode - neg_mode)
    distances <- c(distances, distance)
    unidentified_mat <- snr[(snr$Highest - snr$Second) < distance,]
    unidentified_counts <- c(unidentified_counts, dim(unidentified_mat)[1])
  }
  unid_df <- data.frame(delta_c = dist_frac_list, delta = distances, Unidentifieds = unidentified_counts)

  Pun <- ggplot2::ggplot(unid_df, aes(delta_c, delta)) +
    geom_tile(aes(fill = Unidentifieds), colour = "black") +
    geom_text(aes(label=Unidentifieds)) +
    ggthemes::theme_clean(base_size = 12) +
    scale_x_continuous(breaks = unique(unid_df$delta_c)) +
    scale_y_continuous(breaks = unique(unid_df$delta), labels = scales::number_format(accuracy = 0.001)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(unid_df$Unidentifieds)) +
    ggtitle("Delta vs. Delta_c and resulting unidentified classifications") +
    theme(
      plot.background = element_blank(),
      panel.grid.major.x = element_line(colour = "gray", linetype = "dotted")
    )
  
  print(Pun)
  
  P1 <- (ggplot2::ggplot(snr, aes(x=Highest, y=Second)) +
           geom_point(cex=0.25) +
           egg::theme_presentation(base_size = 12)) + theme(legend.position = c(0.1, 0.65), legend.text=element_text(size=10)) +
    geom_vline(xintercept = neg_df$alpha[-1], size=0.25) + geom_vline(xintercept = neg_df$alpha[1], color="red", size=0.25) + 
    geom_hline(yintercept = doublet_df$beta[-1], size=0.25) + geom_hline(yintercept = doublet_df$beta[1], color="red", size=0.25) + 
    geom_abline(slope=1, intercept = -unid_df$delta, size=0.25) +
    geom_label(aes(x=Inf, y=-Inf, hjust=1, vjust=-0.1, label="Red lines: BQN thresholds 
    Black lines: Parameter effects 
    Values: (0.05, 0.1, 0.15, 0.2, 0.25, 0.5) 
    Beta_c (horizontal): Top: 0.05, Bottom: 0.5 
    Alpha_c (vertical): Left: 0.05, Right: 0.5 
    Delta_c (diag): Top: 0.05, Bottom: 0.5 "), size=2.5)
  
  
  P2 <- suppressWarnings(ggExtra::ggMarginal(P1, size=6, groupColour = FALSE))
  print(P2)
}

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

getCountCutoff <- function(data, label, num_deriv_peaks, barcodeBlocklist = NULL, random.seed = GetSeed()) {
  if (!is.null(random.seed)) {
    set.seed(random.seed)
  }

  j <- 1
  max_list <- numeric(10)

  # adding a small amount of noise (0-1) to count data removes distortions from the log-scale count histogram
  data <- log10(as.vector(data+1) + matrix(stats::runif(dim(data)[[1]]*dim(data)[[2]]), nrow=dim(data)[[1]]))
  data <- data[data > 0]
  # Iterate KDE over sequentially larger bandwidths until fake peaks are smoothed away.
  # Stop when there are fewer peaks than the number specified in num_deriv_peaks.
  while (length(max_list) > num_deriv_peaks) {
    j <- j + 1
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
    
    plotdata <- data.frame(Value = data)
    linedata <- data.frame(x = smooth$x, y = sqrt(smooth$y))
  }

  peak_df <- data.frame(index = max_list, dens = smooth$y[max_list], count = smooth$x[max_list])
  peak_df <- peak_df[order(-peak_df$dens),]
  if ((length(peak_df$index)<2)|(peak_df$dens[1] > 100*peak_df$dens[2])) {
    print(paste0('Only one peak found, using max value as cutoff: ', label))
    ind_min <- length(smooth$x)
    ind1 <- peak_df$index[1]
    ind2 <- ind_min
  } else if (length(peak_df$index) == 2) {
    top2_df <- peak_df
    top2_df <- top2_df[order(top2_df$index),]
    ind1 <- top2_df$index[1]
    ind2 <- top2_df$index[2]
    minval <- min(smooth$y[ind1:ind2])
    ind_range <- which(smooth$y[ind1:ind2] == minval)
    ind_min <- ind_range + ind1 - 1
  } else if (length(peak_df$index) > 2){
    if (peak_df$dens[2] > 10*peak_df$dens[3]) {
      top2_df <- peak_df[1:2,]
      top2_df <- top2_df[order(top2_df$index),]
      ind1 <- top2_df$index[1]
      ind2 <- top2_df$index[2]
      minval <- min(smooth$y[ind1:ind2])
      ind_range <- which(smooth$y[ind1:ind2] == minval)
      ind_min <- ind_range + ind1 - 1
    } else {
      top2_df_a <- peak_df[1:2,]
      top2_df_a <- top2_df_a[order(top2_df_a$index),]
      ind1_a <- top2_df_a$index[1]
      ind2_a <- top2_df_a$index[2]
      minval_a <- min(smooth$y[ind1_a:ind2_a])
      diff_a <- min(top2_df_a$dens) - minval_a
      
      top2_df_b <- peak_df[c(1,3),]
      top2_df_b <- top2_df_b[order(top2_df_b$index),]
      ind1_b <- top2_df_b$index[1]
      ind2_b <- top2_df_b$index[2]
      minval_b <- min(smooth$y[ind1_b:ind2_b])
      diff_b <- min(top2_df_b$dens) - minval_b
      
      if (diff_a > diff_b) {
        ind_range <- which(smooth$y[ind1_a:ind2_a] == minval_a)
        ind_min <- ind_range + ind1_a - 1
        ind1 <- ind1_a
        ind2 <- ind2_a
      } else {
        ind_range <- which(smooth$y[ind1_b:ind2_b] == minval_b)
        ind_min <- ind_range + ind1_b - 1
        ind1 <- ind1_b
        ind2 <- ind2_b
      }
    }
  }
  
  cutoff <- smooth$x[ind_min]
  neg_mode <- smooth$x[ind1]
  pos_mode <- smooth$x[ind2]
 
  if (cutoff == max(smooth$x)) {
    barcodeBlocklist <- c(barcodeBlocklist, label)
  }

  return(list(
    cutoff = cutoff,
    neg_mode = neg_mode,
    pos_mode = pos_mode,
    barcodeBlocklist = barcodeBlocklist,
    plotdata = plotdata,
    linedata = linedata
  ))
}

generateBFFGridPlot <- function(barcodeMatrix, xlab, maintitle, universal_cutoff = NULL) {
  barcodeBlocklist <- NULL
  plotdata <- NULL
  linedata <- NULL
  cutoffs <- NULL
  poslist <- NULL
  neglist <- NULL
  discrete <- barcodeMatrix
  discrete[discrete > 0] <- 0

  for (hto in rownames(barcodeMatrix)) {
    cells <- barcodeMatrix[hto, colnames(barcodeMatrix), drop = FALSE]
    cutoffresults <- getCountCutoff(cells, hto, 3, barcodeBlocklist)
    cutoffval <- cutoffresults[['cutoff']]
    pos_mode <- cutoffresults[['pos_mode']]
    neg_mode <- cutoffresults[['neg_mode']]

    if (!is.null(universal_cutoff)) {
      cutoffval <- universal_cutoff
    }

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
    
    toAdd <- data.frame(pos_mode = pos_mode)
    toAdd$Barcode <- hto
    
    if (is.null(poslist)) {
      poslist <- toAdd
    } else {
      poslist <- rbind(toAdd, poslist)
    }

    toAdd <- data.frame(neg_mode = neg_mode)
    toAdd$Barcode <- hto
    
    if (is.null(neglist)) {
      neglist <- toAdd
    } else {
      neglist <- rbind(toAdd, neglist)
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
    print(ggplot2::ggplot(plotdata, aes(x = Value)) +
            geom_line(data = linedata, mapping = aes(x = x, y = y), color = "blue", size = 1) +
            egg::theme_presentation(base_size = 12) +
            geom_line(data=cutoffs, aes(x=cutoff, y = y), size = 1, color = "red") +
            geom_histogram(aes(y = sqrt(..density..)), size = 1, bins = nbins) +
            scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
            labs(y = 'sqrt(Density)', x = xlab) + ggtitle(maintitle)  +
            ggforce::facet_wrap_paginate(~Barcode, scales = 'free', strip.position = 'top', nrow = min(3, length(unique(plotdata$Barcode))), labeller = labeller(.multi_line = FALSE), page = i)
    )
  }

  return(list(
    barcodeBlocklist = barcodeBlocklist,
    cutoffslist = cutoffslist,
    discrete = discrete,
    poslist = poslist,
    neglist = neglist
  ))
}


GenerateCellHashCallsBFF <- function(barcodeMatrix, assay = "HTO", min_average_reads = 10, verbose = TRUE, simple_threshold = FALSE, doublet_thresh = 0.05, neg_thresh = 0.05, dist_frac = 0.1, metricsFile = NULL, doTSNE = TRUE){
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
    seuratObj <- suppressWarnings(Seurat::CreateSeuratObject(barcodeMatrix, assay = assay))
    seuratObj <- BFFDemux(seuratObj = seuratObj, assay = assay, simple_threshold = simple_threshold, doublet_thresh = doublet_thresh, neg_thresh = neg_thresh, dist_frac=dist_frac, metricsFile = metricsFile)
    if (as.logical(simple_threshold) == TRUE) {
      SummarizeHashingCalls(seuratObj, label = "bff_raw", columnSuffix = "bff_raw", assay = assay, doTSNE = doTSNE, doHeatmap = F)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = "bff_raw", classification = seuratObj$classification.bff_raw, classification.global = seuratObj$classification.global.bff_raw, stringsAsFactors = FALSE)
      df <- .RestoreUnderscoreToHtoNames(df, rownames(barcodeMatrix))
      return(df)
    } else {
      SummarizeHashingCalls(seuratObj, label = 'bff_cluster', columnSuffix = 'bff_cluster', assay = assay, doTSNE = doTSNE, doHeatmap = F)
      df <- data.frame(cellbarcode = as.factor(colnames(seuratObj)), method = 'bff_cluster', classification = seuratObj$classification.bff_cluster, classification.global = seuratObj$classification.global.bff_cluster, stringsAsFactors = FALSE)
      df <- .RestoreUnderscoreToHtoNames(df, rownames(barcodeMatrix))
      return(df)
    }
  }, error = function(e){
    warning('Error generating BFF calls, aborting')
    print(conditionMessage(e))
    traceback()
    return(NULL)
  })
}

BFFDemux <- function(seuratObj, assay, simple_threshold=simple_threshold, doublet_thresh=doublet_thresh, neg_thresh=neg_thresh, dist_frac=dist_frac, metricsFile = NULL) {
  barcodeMatrix <- GetAssayData(
    object = seuratObj,
    assay = assay,
    slot = 'counts'
  )[, colnames(x = seuratObj)]
  
  if (!simple_threshold) {
    print('Running BFF_cluster')
    print(paste("Doublet threshold: ", doublet_thresh))
    print(paste("Neg threshold: ", neg_thresh))
    print(paste("Min distance as fraction of distance between peaks: ", dist_frac))
  } else {
    print('Running BFF_raw')
  }

  thresholdres <- generateBFFGridPlot(barcodeMatrix, "Log(Counts + 1)", "Raw Count Distributions with BQN Thresholds")
  
  cutofflist <- thresholdres[['cutoffslist']]
  cutoffs <- list()
  for (bar in names(cutofflist)) {
    cutoffs[[bar]] <- 10^cutofflist[[bar]]
  }
  
  print("Thresholds:")
  for (hto in names(cutoffs)) {
    print(paste0(hto, ': ', cutoffs[[hto]]))
    .LogMetric(metricsFile, paste0('cutoff.bff.', hto), cutoffs[[hto]])
  }
  
  discrete <- thresholdres[['discrete']]

  if (simple_threshold == TRUE) {
    seuratObj <- .AssignCallsToMatrix(seuratObj, discrete, suffix = 'bff_raw', assay = assay)
    return(seuratObj)
  } else {
    #   Cluster method compares the distance between the top 2 normalized barcode counts to the
    #   distance between the nearest peak and the threshold.  The bc count for a cell considered for negative 
    #   recovery must be close to the threshold.  Likewise, the 2nd highest bc count for a cell considered for
    #   doublet recovery must be close to the threshold.

    # outputs from NormalizeBimodalQuantile: discrete, tot_normed, log10(tot_normed+1), barcodeBlocklist
    normedres <- NULL
    tryCatch({
      normedres <- NormalizeBimodalQuantile(barcodeMatrix)
    }, error = function(e){
      print("No valid barcodes, skipping BFF")
      print(conditionMessage(e))
      traceback()
    })

    if (is.null(normedres)){
      return(NULL)
    }
    
    discrete <- normedres[['discrete']]
    tot_normed <- normedres[['tot_normed']]
    lognormedcounts <- normedres[['lognormedcounts']]

    #generateBFFGridPlot outputs barcodeBlocklist, cutoffslist, discrete, x_vals
    normedplotres <- generateBFFGridPlot(t(tot_normed), "Log(Counts + 1)", "Normalized Count Distributions with Fitted Threshold")
    normed_cutoffs <- normedplotres[['cutoffslist']]
    neg_mode <- mean(normedplotres[['neglist']]$neg_mode)
    pos_mode <- mean(normedplotres[['poslist']]$pos_mode)
    
    snr <- SNR(lognormedcounts)
    singlets <- c()
    negs <- c()
    doublets <- c()
    singlethi <- c()
    doublethi <- c()
    doubletsecond <- c()
    
    highest_dist <- stats::density(snr$Highest, adjust = 1, kernel = 'gaussian',
                                   bw = 'SJ', give.Rkern = FALSE)
    second_dist <- stats::density(snr$Second, adjust = 1, kernel = 'gaussian',
                                  bw = 'SJ', give.Rkern = FALSE)
    
    vals <- c()
    for (i in 2:length(highest_dist$y)) {
      if (highest_dist$y[[i-1]] <= neg_thresh*max(highest_dist$y) & highest_dist$y[[i]] > neg_thresh*max(highest_dist$y)) {
        vals <- c(vals, i)
      }
    }
    if (length(vals) == 0) {
      print('Cannot find negative threshold.  Exiting BFF')
      return(NULL)
    }
    val <- min(vals)
    neg_cutoff <- highest_dist$x[[val]]
    
    
    vals <- c()
    for (i in 2:length(second_dist$y)) {
      if (second_dist$y[[i-1]] > doublet_thresh*max(second_dist$y) & second_dist$y[[i]] <= doublet_thresh*max(second_dist$y)) {
        vals <- c(vals, i)
      }
    }
    if (length(vals) == 0) {
      print('Cannot find doublet threshold.  Exiting BFF')
      return(NULL)
    }
    val <- max(vals)
    doublet_cutoff <- second_dist$x[[val]]

    classification <- c()
    for (i in 1:nrow(snr)) {
      if (snr[i, "Highest"] >= neg_cutoff) {
        if (snr[i, "Second"] <= doublet_cutoff) {
          if (snr[i, "Highest"] - snr[i, "Second"] >= dist_frac * (pos_mode - neg_mode)) {
            singlets <- c(singlets, snr[i, "CellID"])
            classification[i] <- "Singlet"
            singlethi <- c(singlethi, snr[i, "Barcode"])
          } else {
            classification[i] <- "Doublet"
            doublets <- c(doublets, snr[i, "CellID"])
            doublethi <- c(doublethi, snr[i, "Barcode"])
            doubletsecond <- c(doubletsecond, snr[i, "Barcode2"])
          }
          
        } else {
          classification[i] <- "Doublet"
          doublets <- c(doublets, snr[i, "CellID"])
          doublethi <- c(doublethi, snr[i, "Barcode"])
          doubletsecond <- c(doubletsecond, snr[i, "Barcode2"])
        }
      } else {
        classification[i] <- "Negative"
        negs <- c(negs, snr[i, "CellID"])
      }
    }
    
    i <- 0
    for (cell in singlets) {
      i <- i+1
      discrete[, cell] <- 0
      discrete[singlethi[i], cell] <- 1
    }
    
    i <- 0
    for (cell in doublets) {
      i <- i+1
      discrete[, cell] <- 0
      discrete[c(doublethi[i], doubletsecond[i]), cell] <- 1
    }
    
    for (cell in negs) {
      discrete[, cell] <- 0
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

    seuratObj <- .AssignCallsToMatrix(seuratObj, discrete, suffix = 'bff_cluster', assay = assay)
    seuratObj@misc[['cutoffs']] <- cutoffs
    return(seuratObj)
  }
}
