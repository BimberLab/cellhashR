#' @include Utils.R
#' @include Visualization.R

utils::globalVariables(
  names = c('classification', 'classification.global', 'HTO', 'Count', 'cellbarcode', 'Classification', 'consensuscall', 'consensuscall.global', 'topFraction', 'totalReadsPerCell'),
  package = 'cellhashR',
  add = TRUE
)

utils::globalVariables(
  names = c('sortOrder'),
  package = 'cellhashR',
  add = TRUE
)

#' @title AppendCellHashing
#'
#' @param seuratObj The seurat object
#' @param barcodeCallFile The tsv containing cell hashing calls
#' @param barcodePrefix A prefix to be applied before the cell barcodes
#' @description Appends cell hashing calls to a seurat object
#' @return A modified Seurat object.
#' @importFrom dplyr arrange
AppendCellHashing <- function(seuratObj, barcodeCallFile, barcodePrefix = NULL) {
  initialCells <- ncol(seuratObj)
  print(paste0('Initial cell barcodes in GEX data: ', ncol(seuratObj)))

  if (!file.exists(barcodeCallFile)) stop("Barcode File Not found")

  barcodeCallTable <- utils::read.table(barcodeCallFile, sep = '\t', header = T)
  if (!is.null(barcodePrefix)) {
    barcodeCallTable$cellbarcode<- paste0(barcodePrefix, '_', barcodeCallTable$cellbarcode)

    initialCells <- sum(seuratObj$BarcodePrefix == barcodePrefix)
    print(paste0('Initial cell barcodes in GEX data for prefix: ', initialCells))
  }

  if (sum(duplicated(barcodeCallTable$cellbarcode) != 0)) {
    stop('There were duplicated cell barcodes')
  }
  barcodeCallTable <- unique(barcodeCallTable)

  barcodeCallTable <- barcodeCallTable[barcodeCallTable$consensuscall != 'Negative',]
  if (nrow(barcodeCallTable)==0) stop("Something is wrong, table became 0 rows")

  print(paste0('Non-negative cell barcodes in cell hashing calls: ', nrow(barcodeCallTable)))

  # Dont overwrite in case we already added data for another dataset
  if (!('HTO' %in% names(seuratObj@meta.data))) {
    print('Adding HTO columns to seurat object')
    seuratObj$HTO <- c(NA)
    seuratObj$consensuscall.global <- c(NA)
  }

  datasetSelect <- seuratObj$BarcodePrefix == barcodePrefix
  df <- data.frame(cellbarcode = colnames(seuratObj)[datasetSelect])
  df$sortOrder = 1:nrow(df)

  bcIntersect <- barcodeCallTable[barcodeCallTable$cellbarcode %in% df$cellbarcode,]
  pct <- round(nrow(bcIntersect) / nrow(barcodeCallTable) * 100, 2)
  pct2 <- round(nrow(bcIntersect) / nrow(df) * 100, 2)

  print(paste0('Barcodes with calls: ', nrow(barcodeCallTable), ', intersecting with GEX data (total ', nrow(df),'): ', nrow(bcIntersect), " (", pct, "% / ", pct2, "%)"))
  if (nrow(bcIntersect) == 0) {
    print('no barcodes shared')
    print(paste0('first GEX barcodes:'))
    print(head(df$cellbarcode))
    print(paste0('first Hashing barcodes:'))
    print(head(barcodeCallTable$cellbarcode))
  }

  df <- merge(df, barcodeCallTable, all.x = T, all.y = F, by = c('cellbarcode'))
  df <- dplyr::arrange(df, sortOrder)
  df <- df[names(df) != 'sortOrder']

  if (sum(datasetSelect) != nrow(df)) {
    stop('Length of data select and df do not match!')
  }

  df$consensuscall <- as.character(df$consensuscall)
  df$consensuscall.global <- as.character(df$consensuscall.global)

  df$consensuscall[is.na(df$consensuscall)] <- 'ND'
  df$consensuscall.global[is.na(df$consensuscall.global)] <- 'ND'
  print(paste0('Total barcodes added: ', nrow(df[!is.na(df$consensuscall),])))

  # Check barcodes match before merge
  if (sum(df$cellbarcode != colnames(seuratObj)[seuratObj$BarcodePrefix == barcodePrefix]) > 0) {
    stop(paste0('Seurat and cell hashing barcodes do not match after merge, differences: ', sum(df$cellbarcode != colnames(seuratObj)[seuratObj$BarcodePrefix == barcodePrefix])))
  }

  HTO <- as.character(seuratObj$HTO)
  consensuscall.global <- as.character(seuratObj$consensuscall.global)

  HTO[datasetSelect] <- df$HTO
  consensuscall.global[datasetSelect] <- df$consensuscall.global

  seuratObj$consensuscall <- as.factor(consensuscall)
  seuratObj$consensuscall.global <- as.factor(consensuscall.global)

  return(seuratObj)
}

#' @title GenerateCellHashingCalls
#'
#' @param barcodeMatrix The filtered matrix of hashing count data
#' @param methods A vector of one or more calling methods to use. Currently supported are: htodemux and multiseq
#' @param cellbarcodeWhitelist A vector of expected cell barcodes. This allows reporting on the total set of expected barcodes, not just those in the filtered count matrix.
#' @param htodemux.positive.quantile The positive.quantile passed to HTODemux
#' @param metricsFile If provided, summary metrics will be written to this file.
#' @description The primary methods to generating cell hashing calls from a filtered matrix of count data
#' @return A data frame of results.
#' @export
GenerateCellHashingCalls <- function(barcodeMatrix, methods = c('htodemux', 'multiseq'), cellbarcodeWhitelist = NULL, htodemux.positive.quantile = 0.95, metricsFile = metricsFile) {
  callList <- list()
  for (method in methods) {
    if (method == 'htodemux') {
      calls <- GenerateCellHashCallsSeurat(barcodeMatrix, positive.quantile = htodemux.positive.quantile)
      if (!is.null(calls)) {
        callList[[method]] <- calls
      }
    } else if (method == 'multiseq'){
      calls <- GenerateCellHashCallsMultiSeq(barcodeMatrix)
      if (!is.null(calls)) {
        callList[[method]] <- calls
      }
    }
  }

  return(ProcessEnsemblHtoCalls(callList, cellbarcodeWhitelist = cellbarcodeWhitelist, metricsFile = metricsFile))
}

#' @title ProcessEnsemblHtoCalls
#'
#' @import ggplot2
#' @import patchwork
#' @param callList A list of data frames, produced by callers
#' @param metricsFile If provided, summary metrics will be written to this file.
#' @param cellbarcodeWhitelist A vector of expected cell barcodes. This allows reporting on the total set of expected barcodes, not just those in the filtered count matrix.
#' @importFrom dplyr %>% group_by summarise
ProcessEnsemblHtoCalls <- function(callList, cellbarcodeWhitelist = NULL, metricsFile = metricsFile) {
  print('Generating consensus calls')

  if (length(callList) == 0){
    print('No algorithms produced calls, aborting')
    return()
  }

  methods <- names(callList)
  allCalls <- NULL
  for (method in methods) {
    if (all(is.null(allCalls))) {
      allCalls <- callList[[method]]
    } else {
      allCalls <- rbind(allCalls, callList[[method]])
    }
  }

  # Because cell barcodes might have been filtered from the original whitelist (i.e. total counts), re-add:
  if (!is.null(cellbarcodeWhitelist)) {
    allCalls$classification <- as.character(allCalls$classification)
    allCalls$classification.global <- as.character(allCalls$classification.global)

    for (method in unique(allCalls$method)) {
      dat <- allCalls[allCalls$method == method,]
      toAdd <- cellbarcodeWhitelist[!(cellbarcodeWhitelist %in% unique(dat$cellbarcode))]
      if (length(toAdd) > 0) {
        print(paste0('Re-adding ', length(toAdd), ' cell barcodes to call list for: ', method))
        allCalls <- rbind(allCalls, data.frame(cellbarcode = toAdd, method = method, classification = 'Low Counts', classification.global = 'Low Counts'))
        .LogMetric(metricsFile, 'CellBarcodesFilteredForCounts', length(toAdd))
      }
    }

    allCalls$classification <- naturalsort::naturalfactor(allCalls$classification)
    allCalls$classification.global <- naturalsort::naturalfactor(allCalls$classification.global)
  }

  dataClassification <- allCalls[c('cellbarcode', 'method', 'classification')] %>% tidyr::pivot_wider(id_cols = cellbarcode, names_from = method, values_from = classification, values_fill = 'Negative')
  dataClassificationGlobal <- allCalls[c('cellbarcode', 'method', 'classification.global')] %>% tidyr::pivot_wider(id_cols = cellbarcode, names_from = method, values_from = classification.global, values_fill = 'Negative')

  for (method in methods) {
    if (!(method %in% names(dataClassification))) {
      dataClassification[method] <- 'Not Called'
    }

    if (!(method %in% names(dataClassificationGlobal))) {
      dataClassificationGlobal[method] <- 'Not Called'
    }

    if (!is.null(metricsFile)) {
      .LogMetric(metricsFile, paste0('Singlet.', method), sum(dataClassificationGlobal[[method]] == 'Singlet'))
      .LogMetric(metricsFile, paste0('Doublet.', method), sum(dataClassificationGlobal[[method]] == 'Doublet'))
      .LogMetric(metricsFile, paste0('Negative.', method), sum(dataClassificationGlobal[[method]] == 'Negative'))
    }
  }

  MakeConsensusCall <- function(x) {
    x <- unique(x)
    if (length(x) == 1) {
      return(x[1])
    }

    x <- x[x != 'Negative']
    if (length(x) == 1) {
      return(x[1])
    }

    return('Discordant')
  }

  # Concordance across all callers:
  dataClassification$consensuscall <- apply(dataClassification[,methods], 1, MakeConsensusCall)
  dataClassificationGlobal$consensuscall <- apply(dataClassificationGlobal[,methods], 1, MakeConsensusCall)

  #Summary plots:
  summary <- dataClassification[c('cellbarcode', 'consensuscall')]
  summary$method <- 'consensus'
  names(summary) <- c('cellbarcode', 'classification', 'method')
  summary <- rbind(allCalls[c('cellbarcode', 'classification', 'method')], summary)
  P1 <- ggplot(summary, aes(x = classification, group = method, fill = method)) +
    geom_bar(position = position_dodge2(preserve = 'single')) +
    egg::theme_presentation(base_size = 14) +
    labs(x = '', y = 'Cells', fill = 'Classification') +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  summary <- dataClassificationGlobal[c('cellbarcode', 'consensuscall')]
  summary$method <- 'consensus'
  names(summary) <- c('cellbarcode', 'classification.global', 'method')
  summary <- rbind(allCalls[c('cellbarcode', 'classification.global', 'method')], summary)
  P2 <- ggplot(summary, aes(x = method, group = classification.global, fill = classification.global)) +
    geom_bar() +
    egg::theme_presentation(base_size = 14) +
    labs(x = '', y = 'Cells', fill = 'Classification') +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  print(P2 + P1 + plot_layout(widths = c(1, 2)))

  print(paste0('Total concordant: ', sum(dataClassification$consensuscall != 'Discordant')))
  print(paste0('Total discordant (barcode call): ', sum(dataClassification$consensuscall == 'Discordant')))
  print(paste0('Total discordant (global classification): ', sum(dataClassificationGlobal$consensuscall == 'Discordant')))

  if (!is.null(metricsFile)) {
    .LogMetric(metricsFile, 'TotalSinglet', sum(dataClassificationGlobal$consensuscall == 'Singlet'))
    .LogMetric(metricsFile, 'TotalDoublet', sum(dataClassificationGlobal$consensuscall == 'Doublet'))
    .LogMetric(metricsFile, 'TotalNegative', sum(dataClassificationGlobal$consensuscall == 'Negative'))
    .LogMetric(metricsFile, 'TotalDiscordant', sum(dataClassificationGlobal$consensuscall == 'Discordant'))
  }

  discord <- dataClassification[dataClassification$consensuscall == 'Discordant',]
  if (nrow(discord) > 0) {
    combos <- .GetAllCombinations(methods)

    for (j in 1:nrow(combos)) {
      method1 <- as.character(combos$method1[j])
      method2 <- as.character(combos$method2[j])

      data <- dataClassification[c(method1, method2)]
      data$consensuscall <- apply(data[,c(method1, method2)], 1, MakeConsensusCall)
      data$DisagreeWithNeg <- data[[method1]] != data[[method2]]
      data <- data[data$DisagreeWithNeg,]
      if (nrow(data) == 0) {
        next
      }

      discord <- data[data$consensuscall == 'Discordant',] %>% dplyr::group_by_at(c(method1, method2)) %>% summarise(Count = dplyr::n())
      P1 <- ggplot(discord, aes_string(x = method1, y = method2, fill = 'Count')) +
        geom_tile() +
        geom_text(aes(label = Count)) +
        egg::theme_presentation(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = 'none'
        ) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
        ylab(method2) + xlab(method1)

      discordWithNeg <- data %>% dplyr::group_by_at(c(method1, method2)) %>% summarise(Count = dplyr::n())
      P2 <- ggplot(discordWithNeg, aes_string(x = method1, y = method2, fill = 'Count')) +
        geom_tile() +
        geom_text(aes(label = Count)) +
        egg::theme_presentation(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = 'none'
        ) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
        ylab(method2) + xlab(method1)

      print(P1 + P2 + plot_annotation(title = paste0('Discordant calls: ', method1, ' vs. ', method2)))
    }
  }

  toAdd <- dataClassificationGlobal[c('cellbarcode', 'consensuscall')]
  names(toAdd) <- c('cellbarcode', 'consensuscall.global')

  dataClassification <- merge(dataClassification, toAdd, by = 'cellbarcode', all.x = T)

  #final outcome
  df <- data.frame(prop.table(table(Barcode = dataClassification$consensuscall)))
  P1 <- ggplot(df, aes(x = '', y = Freq, fill = Barcode)) +
    geom_bar(width = 1, color = "black", stat = "identity") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = GetPlotColors(length(unique(dataClassification$consensuscall)))) +
    theme_minimal() +
    theme(
      axis.text.x=element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid  = element_blank()
    )

  df <- data.frame(table(Classification = dataClassification$consensuscall.global))
  P2 <- ggplot(df, aes(x = '', y = Freq, fill = Classification)) +
    geom_bar(width = 1, color = "black", stat = "identity") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = GetPlotColors(length(unique(dataClassification$consensuscall.global)))) +
    theme_minimal() +
    theme(
      axis.text.x=element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid  = element_blank()
    )

  print(P1 + P2 + plot_annotation(title = paste0('Final Calls: ', nrow(dataClassification), ' cells')))

  return(dataClassification)
}

.GetAllCombinations <- function(methods){
  ret <- data.frame(method1 = character(), method2 = character())

  encountered <- character()
  for (method1 in methods) {
    for (method2 in methods) {
      if (method2 %in% encountered || method1 == method2) {
        next
      }

      ret <- rbind(ret, data.frame(method1 = method1, method2 = method2))
    }

    encountered <- c(encountered, method1)
  }

  return(ret)
}

#' @title GetExampleMarkdown
#'
#' @description Save a template R markdown file, showing usage of this package
#' @param dest The destination filepath, where the file will be saved
#' @export
GetExampleMarkdown <- function(dest) {
  source <- system.file("rmd/cellhashR.rmd", package = "cellhashR")
  if (!file.exists(source)) {
    stop(paste0('Unable to find file: ', source))
  }

  dest <- normalizePath(dest, mustWork = F)
  success <- file.copy(source, dest, overwrite = TRUE)

  if (!success) {
    stop(paste0('Unable to copy file to: ', dest))
  }
}

#' @title CallAndGenerateReport
#'
#' @description Runs the default processing pipeline
#' @param rawCountData The input barcode file or umi_count folder
#' @param reportFile The file to which the HTML report will be written
#' @param callFile The file to which the table of calls will be written
#' @param barcodeWhitelist A vector of barcode names to retain.
#' @param cellbarcodeWhitelist Either a vector of expected barcodes (such as all cells with passing gene expression data), or the string 'inputMatrix'. If the latter is provided, the set of cellbarcodes present in the original unfiltered count matrix will be stored and used for reporting. This allows the report to count cells that were filtered due to low counts separately from negative/non-callable cells.
#' @param methods The set of methods to use for calling. See GenerateCellHashingCalls for options.
#' @param citeSeqCountDir This is the root folder of the Cite-seq-Count output, containing umi_count and read_count folders. If provided, this will be used to generate a library saturation plot
#' @param minCountPerCell Cells (columns) will be dropped if their total count is less than this value.
#' @param metricsFile If provided, summary metrics will be written to this file.
#' @param title A title for the HTML report
#' @export
CallAndGenerateReport <- function(rawCountData, reportFile, callFile, barcodeWhitelist = NULL, cellbarcodeWhitelist = 'inputMatrix', methods = c('multiseq', 'htodemux'), citeSeqCountDir = NULL, minCountPerCell = 5, title = NULL, metricsFile = NULL) {
  rmd <- system.file("rmd/cellhashR.rmd", package = "cellhashR")
  if (!file.exists(rmd)) {
    stop(paste0('Unable to find file: ', rmd))
  }

  paramList <- list()
  if (!is.null(title)) {
    paramList[['doc_title']] <- title
  }

  rawCountData <- normalizePath(rawCountData)
  if (!is.null(citeSeqCountDir)) {
    citeSeqCountDir <- normalizePath(citeSeqCountDir)
  }

  reportFile <- normalizePath(reportFile, mustWork = F)
  callFile <- normalizePath(callFile, mustWork = F)

  # Use suppressWarnings() to avoid 'MathJax doesn't work with self_contained' warning:
  suppressWarnings(rmarkdown::render(output_file = reportFile, input = rmd, params = paramList))

  return(reportFile)
}

#' @title SummarizeCellsByClassification
#'
#' @description Create summary plots to contrast cells based on call-status. This is designed to help inform why specific cells were not called.
#' @param calls The data frame of calls, produced by GenerateCellHashingCalls
#' @param barcodeMatrix The filtered matrix of hashing count data
#' @export
SummarizeCellsByClassification <- function(calls, barcodeMatrix) {
  df <- data.frame(cellbarcode = colnames(barcodeMatrix), totalReadsPerCell = colSums(barcodeMatrix))
  df$topFraction <- apply(sweep(barcodeMatrix, 2, colSums(barcodeMatrix),`/`), 2, function(x){
    max(x)
  })

  df <- merge(calls, df, by = 'cellbarcode', all.x = F)
  df$topFraction[is.na(df$topFraction)] <- 0
  df$totalReadsPerCell[is.na(df$totalReadsPerCell)] <- 0

  P1 <- ggplot(df, aes(x = consensuscall.global, y = totalReadsPerCell, color = consensuscall)) +
    geom_boxplot() +
    geom_jitter() +
    xlab('') +
    ylab('Counts/Cell') + labs(color = 'Call') +
    egg::theme_presentation(base_size = 14) +
    theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

  P2 <- ggplot(df, aes(x = consensuscall.global, y = topFraction, color = consensuscall)) +
    geom_boxplot() +
    geom_jitter() +
    xlab('') +
    ylab('Top Barcode Fraction') + labs(color = 'Call') +
    egg::theme_presentation(base_size = 14) +
    theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

  print(P1 + P2 + plot_layout(guides = 'collect'))

  df2 <- df[df$consensuscall.global %in% c('Singlet', 'Doublet', 'Negative'),]
  out <- grDevices::boxplot.stats(df2$totalReadsPerCell)$out
  out <- out[out > mean(df2$totalReadsPerCell[df2$totalReadsPerCell > 0])]
  df2 <- df2[df2$totalReadsPerCell < min(out),]

  P1 <- ggplot(df2, aes(x = totalReadsPerCell, color = consensuscall)) +
    geom_density() +
    xlab('Counts/Cell') + labs(color = 'Call') +
    egg::theme_presentation(base_size = 14) +
    facet_wrap(. ~ consensuscall.global, ncol = 1)

  print(P1)


  df2$Category <- 'All'

  df3 <- df2[df2$consensuscall.global == 'Negative',]
  df3$Category <- 'Negatives'
  df3 <- rbind(df2, df3)

  P2 <- ggplot(df3, aes(x = totalReadsPerCell, y = topFraction, color = consensuscall)) +
    geom_point() +
    xlab('Counts/Cell') + ylab('Top Barcode Fraction') + labs(color = 'Call') +
    egg::theme_presentation(base_size = 14) +
    facet_grid(. ~ Category)

  print(P2)
}