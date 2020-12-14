#' @include Utils.R
#' @include Visualization.R

utils::globalVariables(
  names = c('classification', 'classification.global', 'HTO', 'Count'),
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
#' @description Appends cell hashing calls to a seurat object
#' @param seuratObj, A Seurat object.
#' @return A modified Seurat object.
#' @importFrom dplyr arrange
AppendCellHashing <- function(seuratObj, barcodeCallFile, barcodePrefix = NULL) {
  initialCells <- ncol(seuratObj)
  print(paste0('Initial cell barcodes in GEX data: ', ncol(seuratObj)))

  if (!file.exists(barcodeCallFile)) stop("Barcode File Not found")

  barcodeCallTable <- utils::read.table(barcodeCallFile, sep = '\t', header = T)
  if (!is.null(barcodePrefix)) {
    barcodeCallTable$CellBarcode <- paste0(barcodePrefix, '_', barcodeCallTable$CellBarcode)

    initialCells <- sum(seuratObj$BarcodePrefix == barcodePrefix)
    print(paste0('Initial cell barcodes in GEX data for prefix: ', initialCells))
  }

  #Hack until we figure this out upstream
  #TODO: Find discordant duplicates add as second col or convert to doublets or some error
  barcodeCallTable <- unique(barcodeCallTable)

  barcodeCallTable <- barcodeCallTable[barcodeCallTable$HTO != 'Negative',]
  if (nrow(barcodeCallTable)==0) stop("Something is wrong, table became 0 rows")

  print(paste0('Non-negative cell barcodes in HTO calls: ', nrow(barcodeCallTable)))

  #dup <- barcodeCallTable[barcodeCallTable$CellBarcode %in% barcodeCallTable$CellBarcode[duplicated(barcodeCallTable$CellBarcode)],]
  #write.table(dup, file='duplicates.txt', sep = '\t', quote = F, row.names = F)

  # Dont overwrite in case we already added data for another dataset
  if (!('HTO' %in% names(seuratObj@meta.data))) {
    print('Adding HTO columns to seurat object')
    seuratObj$HTO <- c(NA)
    seuratObj$HTO_Classification <- c(NA)
  }

  datasetSelect <- seuratObj$BarcodePrefix == barcodePrefix
  df <- data.frame(CellBarcode = colnames(seuratObj)[datasetSelect])
  df$sortOrder = 1:nrow(df)

  bcIntersect <- barcodeCallTable[barcodeCallTable$CellBarcode %in% df$CellBarcode,]
  pct <- round(nrow(bcIntersect) / nrow(barcodeCallTable) * 100, 2)
  pct2 <- round(nrow(bcIntersect) / nrow(df) * 100, 2)

  print(paste0('Barcodes with calls: ', nrow(barcodeCallTable), ', intersecting with GEX data (total ', nrow(df),'): ', nrow(bcIntersect), " (", pct, "% / ", pct2, "%)"))
  if (nrow(bcIntersect) == 0) {
    print('no barcodes shared')
    print(paste0('first GEX barcodes:'))
    print(head(df$CellBarcode))
    print(paste0('first Hashing barcodes:'))
    print(head(barcodeCallTable$CellBarcode))
  }

  df <- merge(df, barcodeCallTable, all.x = T, all.y = F, by = c('CellBarcode'))
  df <- dplyr::arrange(df, sortOrder)
  df <- df[names(df) != 'sortOrder']

  if (sum(datasetSelect) != nrow(df)) {
    stop('Length of data select and df do not match!')
  }

  df$HTO <- as.character(df$HTO)
  df$HTO_Classification <- as.character(df$HTO_Classification)

  df$HTO[is.na(df$HTO)] <- 'ND'
  df$HTO_Classification[is.na(df$HTO_Classification)] <- 'ND'
  print(paste0('Total HTOs added: ', nrow(df[!is.na(df$HTO),])))

  # Check barcodes match before merge
  if (sum(df$cellBarcode != colnames(seuratObj)[seuratObj$BarcodePrefix == barcodePrefix]) > 0) {
    stop(paste0('Seurat and HTO barcodes do not match after merge, differences: ', sum(df$cellBarcode != colnames(seuratObj)[seuratObj$BarcodePrefix == barcodePrefix])))
  }

  HTO <- as.character(seuratObj$HTO)
  HTO_Classification <- as.character(seuratObj$HTO_Classification)

  HTO[datasetSelect] <- df$HTO
  HTO_Classification[datasetSelect] <- df$HTO_Classification

  seuratObj$HTO <- as.factor(HTO)
  seuratObj$HTO_Classification <- as.factor(HTO_Classification)

  return(seuratObj)
}

#' @title GenerateCellHashingCalls
#'
#' @description A description
#' @return A data table of results.
#' @export
GenerateCellHashingCalls <- function(barcodeMatrix, htodemux.positive.quantile = 0.95, methods = c('htodemux', 'multiseq')) {
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

  return(callList)
  #return(ProcessEnsemblHtoCalls(callList, barcodeMatrix, outFile = outFile))
}

utils::globalVariables(
  names = c('HTO_classification.global.MultiSeq', 'HTO_classification.global.Seurat', 'n', 'HTO_classification.MultiSeq', 'HTO_classification.Seurat'),
  package = 'cellhashR',
  add = TRUE
)

#' @title ProcessEnsemblHtoCalls
#'
#' @import ggplot2
#' @param callList A list of dataframes, produced by callers
#' @param barcodeMatrix The barcode count matrix
#' @importFrom dplyr %>% group_by summarise
ProcessEnsemblHtoCalls <- function(callList, barcodeMatrix, cellbarcodeWhitelist = NA) {
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
  if (!is.na(cellbarcodeWhitelist)) {
    toAdd <- cellbarcodeWhitelist[!(cellbarcodeWhitelist %in% unique(allCalls$cellbarcode))]
    if (length(toAdd) > 0) {
      print(paste0('Re-adding ', length(toAdd), ' cell barcodes to call list'))
      toAdd <- merge(data.frame(cellbarcode = toAdd), dataClassification[FALSE,], by = 'cellbarcode', all.x = TRUE)
      dataClassification <- rbind(dataClassification, toAdd)

      toAdd <- merge(data.frame(cellbarcode = toAdd), dataClassificationGlobal[FALSE,], by = 'cellbarcode', all.x = TRUE)
      dataClassificationGlobal <- rbind(dataClassificationGlobal, toAdd)
    }
  }

  for (method in methods) {
    if (!(method %in% names(dataClassification))) {
      dataClassification[method] <- 'Negative'
    }

    if (!(method %in% names(dataClassificationGlobal))) {
      dataClassificationGlobal[method] <- 'Negative'
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
  dataClassification$ConsensusCall <- apply(dataClassification[,methods], 1, MakeConsensusCall)
  dataClassificationGlobal$ConsensusCall <- apply(dataClassificationGlobal[,methods], 1, MakeConsensusCall)

  #TODO: make plots of this:
  print(paste0('Total concordant: ', sum(dataClassification$ConsensusCall != 'Discordant')))
  print(paste0('Total discordant (barcode call): ', sum(dataClassification$ConsensusCall == 'Discordant')))
  print(paste0('Total discordant (global classification): ', sum(dataClassificationGlobal$ConsensusCall != 'Discordant')))

  discord <- dataClassification[dataClassification$ConsensusCall == 'Discordant',]
  if (nrow(discord) > 0) {
    combos <- .GetAllCombinations(methods)

    for (j in 1:nrow(combos)) {
      method1 <- as.character(combos$method1[j])
      method2 <- as.character(combos$method2[j])

      data <- dataClassification[c(method1, method2)]
      data$ConsensusCall <- apply(data[,c(method1, method2)], 1, MakeConsensusCall)
      data$DisagreeWithNeg <- data[method1] != data[method2]
      data <- data[data$DisagreeWithNeg,]
      if (nrow(data) == 0) {
        next
      }

      discord <- data[data$ConsensusCall == 'Discordant',] %>% dplyr::group_by_at(c(method1, method2)) %>% summarise(Count = dplyr::n())
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

  return(dataClassification)

  # df <- data.frame(
  #   TotalSinglet = c(sum(merged$HTO_classification.global.Seurat == 'Singlet'), sum(merged$HTO_classification.global.MultiSeq == 'Singlet'), sum(merged$FinalClassification == 'Singlet')),
  #   ConcordantSinglet = c(sum(merged$Concordant & merged$HTO_classification.global.Seurat == 'Singlet'), sum(merged$Concordant & merged$HTO_classification.global.MultiSeq == 'Singlet'), sum(merged$FinalClassification == 'Singlet'))
  # )
  # rownames(df) <- c('Seurat', 'MultiSeq', 'Final')
  # df <- t(df)
  # print(kableExtra::kbl(df) %>% kableExtra::kable_styling())
  #
  # if (!is.na(allCallsOutFile) && nrow(merged) > 0) {
  #   write.table(merged, file = allCallsOutFile, row.names = F, sep = '\t', quote = F)
  # }
  #
  # if (nrow(merged) > 0){
  #   dt <- data.frame(CellBarcode = merged$Barcode, HTO = merged$FinalCall, HTO_Classification = merged$FinalClassification, key = 'CellBarcode', Seurat = merged$HasSeuratCall, MultiSeq = merged$HasMultiSeqCall)
  #   dt <- PrintFinalSummary(dt, barcodeMatrix)
  #   write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)
  #
  #   return(dt)
  #
  # } else {
  #   print('No rows, not saving ')
  # }
}


utils::globalVariables(
  names = c('HTO_Classification', 'TotalCounts', 'Freq', 'HTO'),
  package = 'cellhashR',
  add = TRUE
)

#' @title PrintFinalSummary
#'
#' @return A modified Seurat object.
#' @importFrom naturalsort naturalfactor
#' @param The data table with calls
#' @param The barcode counts matrix
#' @import ggplot2
PrintFinalSummary <- function(df, barcodeMatrix){
  #Append raw counts:
  bc <- t(barcodeMatrix)
  x <- reshape2::melt(bc)
  names(x) <- c('CellBarcode', 'HTO', 'Count')

  merged <- merge(df, x, by = c('CellBarcode', 'HTO'), all.x = T, all.y = F)

  bc <- as.data.frame(bc)
  bc$CellBarcode <- rownames(bc)
  merged <- merge(merged, bc, by = c('CellBarcode'), all.x = T, all.y = F)

  merged$HTO <- as.character(merged$HTO)
  merged$HTO[is.na(merged$HTO)] <- c('Negative')
  merged$HTO <- as.factor(merged$HTO)

  merged$HTO_Classification <- as.character(merged$HTO_Classification)
  merged$HTO_Classification[is.na(merged$HTO_Classification)] <- c('Negative')
  merged$HTO_Classification <- as.factor(merged$HTO_Classification)

  #summarise reads by type:
  barcodeMatrix <- as.matrix(barcodeMatrix)
  cs <- colSums(barcodeMatrix)
  cs <- cs[as.character(merged$CellBarcode)]
  merged$TotalCounts <- cs

  htoNames <- SimplifyHtoNames(as.character(merged$HTO))

  merged$HTO <- naturalfactor(as.character(htoNames))

  tbl <- table(SeuratCall = merged$Seurat, MultiSeqCall = merged$MultiSeq)

  colnames(tbl)[colnames(tbl) == T] <- c('MultiSeq Call')
  colnames(tbl)[colnames(tbl) == F] <- c('MultiSeq No Call')

  rownames(tbl)[rownames(tbl) == T] <- c('Seurat Call')
  rownames(tbl)[rownames(tbl) == F] <- c('Seurat No Call')

  print(kableExtra::kbl(tbl) %>% kableExtra::kable_styling())

  print(ggplot(merged, aes(x = HTO)) +
          geom_bar(stat = 'count') +
          xlab('HTO') +
          ylab('Count') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )

  tbl <- table(HTO = merged$HTO)
  df <- data.frame(tbl)
  df$Pct <- round((df$Freq / sum(df$Freq)) * 100, 2)
  print(kableExtra::kbl(df) %>% kableExtra::kable_styling())

  print(ggplot(df, aes(x = '', y=Freq, fill=HTO)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = GetPlotColors(length(unique(df$HTO)))) +
    theme_minimal() +
    theme(
      axis.text.x=element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid  = element_blank()
    ) +
    ggtitle('HTO')
  )

  print(ggplot(merged, aes(x = HTO)) +
          geom_bar(stat = 'count') +
          xlab('HTO') +
          ylab('Count') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )

  cellData <- merged[c('TotalCounts', 'Count', 'HTO_Classification', 'HTO')]
  cellData$TotalCounts <- log(cellData$TotalCounts + 0.5)
  cellData$Count <- log(cellData$Count + 0.5)

  basePlot <- ggplot(cellData, aes(x = TotalCounts, fill = HTO_Classification)) +
    geom_density() +
    scale_fill_manual(values = GetPlotColors(length(unique(cellData$HTO_Classification)))) +
    xlab('Total Counts/Cell (log)') +
    ylab('Density') +
    ggtitle('Total Counts By HTO')


  for (i in 1:GetTotalPlotPages(length(unique(cellData$HTO_Classification)))) {
    print(basePlot + ggforce::facet_grid_paginate(HTO_Classification ~ ., scales = 'free', rows = 4, page = i))
  }

  basePlot <- ggplot(cellData[!(cellData$HTO %in% c('Negative', 'Doublet', 'Discordant')),], aes(x = Count, fill = HTO)) +
    geom_density() +
    scale_fill_manual(values = GetPlotColors(length(unique(cellData$HTO)))) +
    xlab('HTO Counts/Cell (log)') +
    ylab('Density') +
    ggtitle('Counts By HTO')

  for (i in 1:GetTotalPlotPages(length(unique(cellData$HTO)))) {
    print(basePlot + ggforce::facet_grid_paginate(HTO ~ ., scales = 'free', rows = 4, page = i))
  }

  if (sum(merged$HTO == 'Negative') == 0) {
    print('There were no negative cells')
  } else {
    #Melt data:
    melted <- as.data.frame(merged)
    melted <- melted[melted$HTO == 'Negative', !(colnames(melted) %in% c('HTO_Classification', 'HTO', 'key', 'Seurat', 'MultiSeq', 'Count', 'TotalCounts')), drop = FALSE]
    print(utils::str(melted))
    melted <- tidyr::gather(melted, key = 'HTO', value = 'Count', -CellBarcode)

    htoNames <- SimplifyHtoNames(as.character(melted$HTO))
    melted$HTO <- naturalfactor(as.character(htoNames))
    melted$Count <- log10(melted$Count + 0.5)

    basePlot <- ggplot(melted, aes(x = Count, fill = HTO)) +
      geom_density() +
      scale_fill_manual(values = GetPlotColors(length(unique(melted$HTO)))) +
      xlab('HTO Counts/Cell (log)') +
      ylab('Density') +
      ggtitle('Raw HTO Counts For Negative Cells (log10)')

    for (i in 1:GetTotalPlotPages(length(unique(melted$HTO)))) {
      print(basePlot + ggforce::facet_grid_paginate(HTO ~ ., scales = 'free', rows = 4, page = i))
    }
  }

  tbl <- table(HTO_Classification = merged$HTO_Classification)
  df <- data.frame(tbl)
  df$Pct <- round((df$Freq / sum(df$Freq)) * 100, 2)
  print(kableExtra::kbl(df) %>% kableExtra::kable_styling())

  print(ggplot(df, aes(x = '', y=Freq, fill=HTO_Classification)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = GetPlotColors(length(names(tbl)))) +
    theme_minimal() +
    theme(
      axis.text.x=element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid  = element_blank()
    ) +
    ggtitle('HTO Classification')
  )

  return(merged)
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