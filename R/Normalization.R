#' @include Utils.R

utils::globalVariables(
  names = c('NormCount', 'Saturation', 'Cluster', 'AvgExpression', 'cutoff', 'count', 'Normalization'),
  package = 'cellhashR',
  add = TRUE
)

NormalizeQuantile <- function(mat) {
  dat <- preprocessCore::normalize.quantiles(as.matrix(mat), copy=TRUE)
  row.names(dat) <- row.names(mat)
  colnames(dat) <- colnames(mat)
  return(as.data.frame(dat))
}

NormalizeBimodalQuantile <- function(barcodeMatrix) {
  barcodeMatrix <- .EnsureNonSparse(barcodeMatrix)

  cutoffs <- list()
  threshold <- list()
  barcodeBlocklist <- NULL
  for (hto in rownames(barcodeMatrix)) {
    cells <- barcodeMatrix[hto, colnames(barcodeMatrix), drop = FALSE]
    cutoffresults <- getCountCutoff(cells, hto, 4, barcodeBlocklist)

    barcodeBlocklist <- cutoffresults[['barcodeBlocklist']]
    threshold[[hto]] <- 10^(cutoffresults[['cutoff']])
  }
  # barcodeBlocklist, defined in the for loop above, is used to subset
  # barcodeMatrix to exclude htos w/o bimodal distributions
  mat <- barcodeMatrix[!rownames(barcodeMatrix) %in% barcodeBlocklist,]

  if (length(rownames(mat)) == 0) {
    return()
  }

  discrete <- mat
  discrete[discrete > 0] <- 0

  if (!all(rownames(discrete) == rownames(mat))) {
    stop('rownames dont match')
  }

  if (!all(colnames(discrete) == colnames(mat))) {
    stop('colnames dont match')
  }

  if (!identical(dim(discrete), dim(mat))) {
    stop(paste0('discrete/mat have different dimensions, were: ', paste0(dim(discrete), collapse = ','), ' and ', paste0(dim(mat), collapse = ',')))
  }

  for (hto in rownames(mat)) {
    cells <- mat[hto,]
    discrete[hto,] <- ifelse(cells > threshold[[hto]], yes = 1, no = 0)
    cutoffs[[hto]] <- threshold[[hto]]
  }

  if (!identical(dim(discrete), dim(mat))) {
    stop(paste0('discrete/mat have different dimensions after update, were: ', paste0(dim(discrete), collapse = ','), ' and ', paste0(dim(mat), collapse = ',')))
  }

  #NOTE: these could have different dimensions if there is a barcodeBlocklist
  neg_norm <- getNegNormedData(discrete, mat)
  pos_norm <- getPosNormedData(discrete, mat)
  tot_normed <- pos_norm + neg_norm

  return(list(
    discrete = discrete,
    tot_normed = tot_normed,
    lognormedcounts = log10(tot_normed+1),
    barcodeBlocklist = barcodeBlocklist)
  )
}

TransposeDF <- function(df) {
  t_df <- data.table::transpose(df)
  rownames(t_df) <- colnames(df)
  colnames(t_df) <- rownames(df)
  return(t_df)
}

NormalizeLog2 <- function(mat, mean.center = TRUE) {
  log2Scaled <- log2(mat)
  for (i in 1:nrow(log2Scaled)) {
    ind <- which(is.finite(log2Scaled[i,]) == FALSE)
    log2Scaled[i,ind] <- 0

    if (mean.center) {
      log2Scaled[i,] <- log2Scaled[i,] - mean(log2Scaled[i,])
    }
  }

  return(as.matrix(log2Scaled))
}

NormalizeCLR <- function(mat) {
  seuratObj <- Seurat::CreateSeuratObject(mat, assay = 'Hashing')
  seuratObj <- Seurat::NormalizeData(seuratObj, assay = 'Hashing', normalization.method = "CLR", verbose = FALSE)

  return(seuratObj@assays$Hashing@data)
}

NormalizeRelative <- function(mat) {
  return(prop.table(mat, 2))
}

#' @title Plot Normalization QC
#'
#' @param barcodeData The count matrix
#' @param methods The normalizations to use. Allowable options are: bimodalQuantile, Quantile, log2Center, CLR and Raw
#' @description Generates QC plots related to normalization
#' @export
PlotNormalizationQC <- function(barcodeData, methods = c('bimodalQuantile', 'Quantile', 'log2Center', 'CLR')) {
  .LogProgress('Plotting QC')
  toQC <- list()
  for (method in methods) {
    .LogProgress(paste0('Plotting: ', method))
    if (method == 'bimodalQuantile') {
      bqn <- NULL
      
      tryCatch({
        normedres <- NormalizeBimodalQuantile(barcodeData)
        if (!all(is.null(normedres))) {
          temp <- normedres[['lognormedcounts']]
          ParameterScan(temp)
          toQC[['bimodalQuantile']] <- TransposeDF(data.frame(temp, check.names=FALSE))
        } else {
          print('Error running BQN, skipping')
        }
      }, error = function(e){
        print("No valid barcodes, skipping Bimodal quantile normalization")
        print(conditionMessage(e))
        traceback()
      })
    } else if (method  == 'Raw') {
      toQC[['Raw']] <- data.frame(log10(barcodeData + 1), check.names = FALSE)
    } else if (method  == 'Quantile') {
      toQC[['Quantile']] <- TransposeDF(data.frame(log10(NormalizeQuantile(t(barcodeData))+1), check.names = FALSE))
    } else if (method  == 'CLR') {
      toQC[['CLR']] <- NormalizeCLR(barcodeData)
    } else if (method  == 'log2Center') {
      toQC[['log2Center']] <- NormalizeLog2(barcodeData, mean.center = TRUE)
    } else {
      stop(paste0('Unknown method: ', method))
    }

    .LogProgress(paste0('Done normalizing: ', method))
  }
  uniqueRows <- c()
  uniqueCols <- c()

  for (norm in names(toQC)) {
    uniqueRows <- unique(c(uniqueRows, nrow(toQC[[norm]])))
    uniqueCols <- unique(c(uniqueCols, ncol(toQC[[norm]])))

    if (length(uniqueRows) != 1) {
      print(paste0('inconsistent row size in normalized data: ', paste0(uniqueRows, collapse = ','), ', encountered for: ', norm))
    }

    if (length(uniqueCols) != 1) {
      print(paste0('inconsistent column size in normalized data: ', paste0(uniqueCols, collapse = ','), ', encountered for: ', norm))
    }
  }

  df <- NULL
  for (norm in names(toQC)) {
    toAdd <- reshape2::melt(t(toQC[[norm]]))
    names(toAdd) <- c('CellBarcode', 'Barcode', 'NormCount')
    toAdd$Normalization <- norm

    if (is.null(df)) {
      df <- toAdd
    } else {
      df <- rbind(toAdd, df)
    }

    df$Barcode <- SimplifyHtoNames(as.character(df$Barcode))
    df$Barcode <- naturalsort::naturalfactor(df$Barcode)
  }

  maxPerPlot <- 3
  totalPages <- GetTotalPlotPages(totalValues = length(unique(df$Barcode)), perPage = maxPerPlot)
  if ('Raw' %in% levels(df$Normalization)) {
    df$Normalization <- forcats::fct_relevel(Normalization, 'Raw')
  }

  for (i in 1:totalPages) {
    print(ggplot2::ggplot(df, aes(x = NormCount, color = Barcode)) +
      egg::theme_presentation(base_size = 14) +
      geom_density(aes(y = sqrt(after_stat(density))), size = 1) +
      labs(y = 'sqrt(Density)', x = 'Value') + ggtitle('Normalized Data') +
      ggforce::facet_wrap_paginate(Barcode ~ Normalization, scales = 'free', ncol = length(unique(df$Normalization)), nrow = min(3, length(unique(df$Barcode))), strip.position = 'top', drop = FALSE, labeller = labeller(.multi_line = FALSE), page = i)
    )
  }

  for (norm in names(toQC)) {
    if (norm == "Raw") {
      p1title <- "Raw Counts"
      p2title <- "Raw Count Density"
      announcement <- paste0("Printing Raw Count QC Plots")
    } else {
      p1title <- paste0(norm, ' Normalized Counts')
      p2title <- paste0(norm, ' Normalized Count Density')
      announcement <- paste0("Printing ", norm, " Normalization QC Plots")
    }
    print(announcement)
    .PlotViolin(t(toQC[[norm]]), norm = norm)
    PerformHashingClustering(toQC[[norm]], norm = norm)
    snr <- SNR(t(toQC[[norm]]))
    snr$Barcode <- naturalsort::naturalfactor(snr$Barcode)

    P1 <- (ggplot2::ggplot(snr, aes(x=Highest, y=Second, color=Barcode)) +
                  geom_point(cex=0.25) + ggtitle(p1title) +
      egg::theme_presentation(base_size = 14)) + theme(legend.position = c(0.1, 0.65), legend.text=element_text(size=10)) +
      guides(colour = guide_legend(override.aes = list(size=3)))

    P2 <- ggplot(snr, aes(x=Highest, y=Second) ) +
      stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
      scale_fill_distiller(palette=16, direction=-1) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      egg::theme_presentation(base_size = 14) + ggtitle(p2title) + 
      theme(
        legend.position='none'
      )

    P3 <- suppressWarnings(ggExtra::ggMarginal(P1, size=4, groupColour = TRUE))
    print(P2 + P3)
  }

  .LogProgress('Finished plotting QC')
}

PerformHashingClustering <- function(barcodeMatrix, norm) {
  barcodeMatrix <- SeuratObject::as.sparse(barcodeMatrix)
  seuratObj <- CreateSeuratObject(barcodeMatrix, assay = norm)

  # Calculate tSNE embeddings with a distance matrix
  tryCatch({
    perplexity <- .InferPerplexityFromSeuratObj(seuratObj, 100)
    print(paste0('Running tSNE with perplexity: ', perplexity, ' for normalization: ', norm))
    seuratObj[['hto_tsne']] <- RunTSNE(stats::dist(Matrix::t(barcodeMatrix)), assay = norm, perplexity = perplexity)

    .PlotClusters(barcodeMatrix, seuratObj, norm)
  }, error = function(e){
    print('Error generating tSNE, skipping')
    print(conditionMessage(e))
    traceback()
  })
}

.PlotClusters <- function(barcodeMatrix, seuratObj, norm) {
  # clara:
  ncenters <- (nrow(x = barcodeMatrix) + 1)
  init.clusters <- clara(
    x = t(x = barcodeMatrix),
    k = ncenters,
    samples = 100
  )

  if (sum(names(init.clusters$clustering) != colnames(seuratObj)) > 0) {
    stop('Expected cluster names to match cells for clara!')
  }

  seuratObj$cluster.clara <- naturalsort::naturalfactor(init.clusters$clustering)
  P <- suppressWarnings(DimPlot(seuratObj, reduction = 'hto_tsne', group.by = 'cluster.clara', label = TRUE))
  P <- P + ggtitle(paste0('Clusters: ', norm, ' (clara)'))

  Idents(seuratObj) <- 'cluster.clara'
  P2 <- .CreateClusterTilePlot(seuratObj, assay = norm)
  print(P | P2)

  # kmeans:
  init.clusters <- suppressWarnings(stats::kmeans(
    x = t(x = barcodeMatrix),
    centers = ncenters,
    nstart = 100,
    iter.max = 50
  ))

  if (sum(names(init.clusters$cluster) != colnames(seuratObj)) > 0) {
    stop('Expected cluster names to match cells for kmeans!')
  }

  seuratObj$cluster.kmeans <- naturalsort::naturalfactor(init.clusters$cluster)
  P <- suppressWarnings(DimPlot(seuratObj, group.by = 'cluster.kmeans', label = TRUE))
  P <- P + ggtitle(paste0('Clusters: ', norm, ' (kmeans)'))
  Idents(seuratObj) <- 'cluster.kmeans'
  P2 <- .CreateClusterTilePlot(seuratObj, assay = norm)
  print(P | P2)
}

.CreateClusterTilePlot <- function(seuratObj, assay) {
  average.expression <- suppressWarnings(Seurat::AverageExpression(
    object = seuratObj,
    assays = c(assay),
    verbose = FALSE,
    return.seurat = TRUE
  ))

  df <- reshape2::melt(GetAssayData(average.expression, assay = assay))
  names(df) <- c('Cluster', 'Barcode', 'AvgExpression')
  df$Cluster <- naturalsort::naturalfactor(df$Cluster)
  df$Barcode <- naturalsort::naturalfactor(df$Barcode)
  df$AvgExpression <- round(df$AvgExpression, 1)

  P2 <- ggplot(df, aes(Cluster, Barcode)) +
    geom_tile(aes(fill = AvgExpression), colour = "white") +
    geom_text(aes(label=AvgExpression)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = mean(df$AvgExpression)) +
    scale_x_discrete(position = "top") +
    labs(x = "Barcode",y = "Cluster") +
    egg::theme_presentation(base_size = 12) +
    theme(
      legend.position = 'none'
    )

  return(P2)
}

.PlotViolin <- function(df, norm) {
  if (norm == "Raw") {
    label <- 'Log Raw Counts'
    maintitle <- "Raw HTO Barcode Count Distributions"
    
  } else {
    label <- paste0(norm, ' Normalized Counts')
    maintitle <- paste0(norm, " Normalized HTO Barcode Count Distributions")
  }

  df <- data.frame(df, check.names=FALSE)
  df$cell <- rownames(df)
  df <- df %>% tidyr::pivot_longer(colnames(df)[1:length(colnames(df)) - 1], names_to = "Barcode", values_to = "count")
  df$Barcode <- naturalsort::naturalfactor(df$Barcode)
  P1 <- df %>%
    ggplot(aes( y=count, x=Barcode)) + 
    geom_violin(position="dodge", alpha=0.5) +
    xlab("") +
    ylab(label) +
    ggplot2::ggtitle(maintitle) +
    egg::theme_presentation(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90))
  print(P1)
}

