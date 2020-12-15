#' @include Utils.R
#' @include Visualization.R

utils::globalVariables(
  names = c('classification', 'classification.global', 'HTO', 'Count', 'cellbarcode', 'Classification', 'consensuscall'),
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
#' @description The primary methods to generating cell hashing calls from a filtered matrix of count data
#' @return A data frame of results.
#' @export
GenerateCellHashingCalls <- function(barcodeMatrix, methods = c('htodemux', 'multiseq'), cellbarcodeWhitelist = NULL, htodemux.positive.quantile = 0.95) {
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

  return(ProcessEnsemblHtoCalls(callList, barcodeMatrix, cellbarcodeWhitelist))
}

#' @title ProcessEnsemblHtoCalls
#'
#' @import ggplot2
#' @param callList A list of dataframes, produced by callers
#' @param barcodeMatrix The barcode count matrix
#' @param cellbarcodeWhitelist A vector of expected cell barcodes. This allows reporting on the total set of expected barcodes, not just those in the filtered count matrix.
#' @importFrom dplyr %>% group_by summarise
ProcessEnsemblHtoCalls <- function(callList, barcodeMatrix, cellbarcodeWhitelist = NULL) {
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

  dataClassification <- allCalls[c('cellbarcode', 'method', 'classification')] %>% tidyr::pivot_wider(id_cols = cellbarcode, names_from = method, values_from = classification, values_fill = 'Negative')
  dataClassificationGlobal <- allCalls[c('cellbarcode', 'method', 'classification.global')] %>% tidyr::pivot_wider(id_cols = cellbarcode, names_from = method, values_from = classification.global, values_fill = 'Negative')

  # Because cell barcodes might have been filtered from the original whitelist (i.e. total counts), re-add:
  if (!is.null(cellbarcodeWhitelist)) {
    toAdd <- cellbarcodeWhitelist[!(cellbarcodeWhitelist %in% unique(allCalls$cellbarcode))]
    if (length(toAdd) > 0) {
      print(paste0('Re-adding ', length(toAdd), ' cell barcodes to call list'))

      for (method in methods) {
        if (method %in% names(dataClassification)) {
          dataClassification[[method]] <- naturalsort::naturalfactor(dataClassification[[method]], levels = c(levels(dataClassification[[method]]), 'No Counts'))
        }

        if (method %in% names(dataClassificationGlobal)) {
          dataClassificationGlobal[[method]] <- naturalsort::naturalfactor(dataClassificationGlobal[[method]], levels = c(levels(dataClassificationGlobal[[method]]), 'No Counts'))
        }
      }

      merged <- merge(data.frame(cellbarcode = toAdd), dataClassification[FALSE,], by = c('cellbarcode'), all.x = TRUE)
      merged[is.na(merged)] <- 'No Counts'
      dataClassification <- rbind(dataClassification, merged)


      merged <- merge(data.frame(cellbarcode = toAdd), dataClassificationGlobal[FALSE,], by = c('cellbarcode'), all.x = TRUE)
      merged[is.na(merged)] <- 'No Counts'
      dataClassificationGlobal <- rbind(dataClassificationGlobal, merged)
      dataClassificationGlobal[is.na(dataClassificationGlobal)] <- 'No Counts'
    }
  }

  for (method in methods) {
    if (!(method %in% names(dataClassification))) {
      dataClassification[method] <- 'Not Called'
    }

    if (!(method %in% names(dataClassificationGlobal))) {
      dataClassificationGlobal[method] <- 'Not Called'
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

  discord <- dataClassification[dataClassification$consensuscall == 'Discordant',]
  if (nrow(discord) > 0) {
    combos <- .GetAllCombinations(methods)

    for (j in 1:nrow(combos)) {
      method1 <- as.character(combos$method1[j])
      method2 <- as.character(combos$method2[j])

      data <- dataClassification[c(method1, method2)]
      data$consensuscall <- apply(data[,c(method1, method2)], 1, MakeConsensusCall)
      data$DisagreeWithNeg <- data[method1] != data[method2]
print(typeof(data$DisagreeWithNeg))
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

  print(P1 + P2 + plot_annotation(title = 'Final Calls'))

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
  source <- system.file("rmd/cellhashR.Rmd", package = "cellhashR")
  invisible(file.copy(source, dest, overwrite = TRUE))
}

#' @title CallAndGenerateReport
#'
#' @description Runs the default processing pipeline
#' @param rawCountData The input barcode file or umi_count folder
#' @param reportFile The file to which the HTML report will be written
#' @param callFile The file to which the table of calls will be written
#' @param minCountPerCell Cells (columns) will be dropped if their total count is less than this value.
#' @param barcodeWhitelist A vector of barcode names to retain.
#' @param cellbarcodeWhitelist Either a vector of expected barcodes (such as all cells with passing gene expression data), or the string 'matrix'. If the latter is provided, the set of cellbarcodes present in the original unfiltered count matrix will be stored and used for reporting. This allows the report to count cells that were filtered due to low counts separately from negative/non-callable cells.
CallAndGenerateReport <- function(rawCountData, reportFile, callFile, barcodeWhitelist = NULL, cellbarcodeWhitelist = 'matrix') {
  rmd <- system.file("rmd/cellhashR.Rmd", package = "cellhashR")

  rmarkdown::render(output_file = reportFile, input = rmd)

  return(reportFile)
}

