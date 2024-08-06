#' @include Utils.R
#' @include Visualization.R
#' @include Multiseq.R
#' @include Seurat_HTO_Demux.R
#' @include SeqND_Demux.R
#' @include BFF_Demux.R
#' @include Threshold_Demux.R
#' @include DropletUtils_Demux.R
#' @include GMM_Demux.R
#' @include DemuxEM.R
#' @include demuxmix.R
#'
utils::globalVariables(
  names = c('classification', 'classification.global', 'HTO', 'Count', 'cellbarcode', 'Classification', 'consensuscall', 'consensuscall.global', 'topFraction', 'totalReadsPerCell', 'Method', 'DisagreementRate', 'Label', 'TotalPerCell'),
  package = 'cellhashR',
  add = TRUE
)

utils::globalVariables(
  names = c('sortOrder'),
  package = 'cellhashR',
  add = TRUE
)

#' @title Append Cell Hashing to a Seurat Object
#'
#' @param seuratObj The seurat object
#' @param barcodeCallFile The tsv containing cell hashing calls
#' @param barcodePrefix A prefix to be applied before the cell barcodes
#' @importFrom magrittr %>%
#' @description Appends cell hashing calls to a seurat object
#' @return A modified Seurat object.
#' @export
#' @importFrom dplyr arrange
AppendCellHashing <- function(seuratObj, barcodeCallFile, barcodePrefix) {
  if (!file.exists(barcodeCallFile)) {
    stop("Barcode File Not found")
  }

  if (is.null(barcodePrefix)) {
    stop('Must provide the barcodePrefix')
  }

  print(paste0('Initial cell barcodes in GEX data: ', ncol(seuratObj)))

  initialCells <- sum(seuratObj$BarcodePrefix == barcodePrefix)
  print(paste0('Initial cell barcodes in GEX data for prefix: ', initialCells))

  barcodeCallTable <- utils::read.table(barcodeCallFile, sep = '\t', header = T)
  barcodeCallTable$cellbarcode<- paste0(barcodePrefix, '_', barcodeCallTable$cellbarcode)
  if (sum(duplicated(barcodeCallTable$cellbarcode) != 0)) {
    stop('There were duplicated cell barcodes')
  }

  barcodeCallTable <- barcodeCallTable[barcodeCallTable$consensuscall != 'Negative',]
  if (nrow(barcodeCallTable)== 0) {
    stop("There are zero rows with non-negative HTO calls")
	}

  print(paste0('Non-negative cell barcodes in cell hashing calls: ', nrow(barcodeCallTable)))

  # Dont overwrite in case we already added data for another dataset
  if (!('HTO' %in% names(seuratObj@meta.data))) {
    print('Adding HTO columns to seurat object')
    seuratObj$HTO <- c(NA)
    seuratObj$HTO.Classification <- c(NA)
    seuratObj$Saturation.HTO <- c(NA)
  }

  datasetSelect <- seuratObj$BarcodePrefix == barcodePrefix
  df <- data.frame(cellbarcode = colnames(seuratObj)[datasetSelect])
  df$sortOrder <- 1:nrow(df)

  bcIntersect <- barcodeCallTable[barcodeCallTable$cellbarcode %in% df$cellbarcode,]
  pct <- round(nrow(bcIntersect) / nrow(barcodeCallTable) * 100, 2)
  pct2 <- round(nrow(bcIntersect) / nrow(df) * 100, 2)

  print(paste0('Barcodes with calls: ', nrow(barcodeCallTable), ', intersecting with GEX data (total ', nrow(df),'): ', nrow(bcIntersect), " (", pct, "% / ", pct2, "%)"))
  if (nrow(bcIntersect) == 0) {
    print(paste0('First GEX cell barcodes:'))
    print(head(df$cellbarcode))
    print(paste0('First hashing cell barcodes:'))
    print(head(barcodeCallTable$cellbarcode))
    stop('Error: no barcodes shared between GEX and hashing data')
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

  consensuscall <- as.character(seuratObj$HTO)
  consensuscall.global <- as.character(seuratObj$HTO.Classification)

  consensuscall[datasetSelect] <- df$consensuscall
  consensuscall.global[datasetSelect] <- df$consensuscall.global

  seuratObj$HTO <- naturalsort::naturalfactor(consensuscall)
  seuratObj$HTO.Classification <- naturalsort::naturalfactor(consensuscall.global)

  if ('saturation' %in% names(df)) {
    seuratObj$Saturation.HTO <- df$saturation
  }

  return(seuratObj)
}


#' @title Generate Cell Hashing Calls
#'
#' @param barcodeMatrix The filtered matrix of hashing count data
#' @param methods A vector of one or more calling methods to use. Currently supported are: htodemux, multiseq, dropletutils, gmm_demux, demuxem, demuxmix, bff_raw, and bff_cluster
#' @param methodsForConsensus By default, a consensus call will be generated using all methods; however, if this parameter is provided, all algorithms specified by methods will be run, but only the list here will be used for the final consensus call. This allows one to see the results of a given caller without using it for the final calls.
#' @param cellbarcodeWhitelist A vector of expected cell barcodes. This allows reporting on the total set of expected barcodes, not just those in the filtered count matrix.
#' @param metricsFile If provided, summary metrics will be written to this file.
#' @param doTSNE If true, tSNE will be run on the resulting hashing calls after each caller. This can be useful as a sanity check; however, adds time.
#' @param doHeatmap If true, Seurat::HTOHeatmap will be run on the results of each caller. Not supported by all callers.
#' @param perCellSaturation An optional dataframe with the columns cellbarcode and saturation. This will be merged into the final output.
#' @param majorityConsensusThreshold This applies to calculating a consensus call when multiple algorithms are used. If NULL, then all non-negative calls must agree or that cell is marked discordant. If non-NULL, then the number of algorithms returning the top call is divided by the total number of non-negative calls. If this ratio is above the majorityConsensusThreshold, that value is selected. For example, when majorityConsensusThreshold=0.6 and the calls are: HTO-1,HTO-1,Negative,HTO-2, then 2/3 calls are for HTO-1, giving 0.66. This is greater than the majorityConsensusThreshold of 0.6, so HTO-1 is returned. This can be useful for situations where most algorithms agree, but a single caller fails.
#' @param chemistry This string is passed to EstimateMultipletRate. Should be either 10xV2 or 10xV3. This is used to calculate and present the expected doublet rate and does not influence the actual calls.
#' @param callerDisagreementThreshold If provided, the agreement rate will be calculated between each caller and the simple majority call, ignoring discordant and no-call cells. If any caller has an disagreement rate above this threshold, it will be dropped and the consensus call re-calculated. The general idea is to drop a caller that is systematically discordant.
#' @param rawFeatureMatrixH5 If either demuxem or demuxmix are used, you must provide the filepath to the 10x h5 gene expression counts file
#' @param maxAllowableDoubletRate Per caller, the doublet rate will be computed as the total doublets / total droplets (including negatives). Any individual caller with a doublet rate above this value will be converted to NoCall. Note: if 'auto' is chosen, the value will be selected as twice the theoretical doublet rate.
#' @param \dots Caller-specific arguments can be passed by prefixing with the method name. For example, htodemux.positive.quantile = 0.95, will be passed to the htodemux positive.quantile argument).
#' @description The primary methods to generating cell hashing calls from a filtered matrix of count data.
#' @return A data frame of results.
#' @export
GenerateCellHashingCalls <- function(barcodeMatrix, methods = c('bff_cluster', 'gmm_demux', 'dropletutils'), methodsForConsensus = NULL, cellbarcodeWhitelist = NULL, metricsFile = NULL, doTSNE = TRUE, doHeatmap = TRUE, perCellSaturation = NULL, majorityConsensusThreshold = NULL, chemistry = '10xV3', callerDisagreementThreshold = NULL, rawFeatureMatrixH5 = NULL, maxAllowableDoubletRate = 'auto', ...) {
  .LogProgress('Generating calls')
  if (is.data.frame(barcodeMatrix)) {
    print('Converting input data.frame to a matrix')
    barcodeMatrix <- as.matrix(barcodeMatrix)
  }
  else if (!inherits(x = barcodeMatrix, what = 'dgCMatrix') && !is.matrix(barcodeMatrix)) {
    stop(paste0('barcodeMatrix must be a matrix, was: ', class(barcodeMatrix)[1]))
  }

  if (all(is.null(methodsForConsensus))) {
    methodsForConsensus <- methods
  } else {
    methodsForConsensus2 <- intersect(methodsForConsensus, methods)
    if (length(methodsForConsensus2) != length(methodsForConsensus)) {
      stop('All methods in methodsForConsensus must be included in methods')
    }
  }

  if (length(methodsForConsensus) == 0) {
    stop('No methodsForConsensus were provided for consensus calling! Either leave blank (NULL), or use a subset of the methods provided under the methods argument')
  }

  if (('demuxem' %in% methods || 'demuxmix' %in% methods) && is.null(rawFeatureMatrixH5)) {
    stop('If either demuxem or demuxmix are used, you must supply rawFeatureMatrixH5')
  }

  callList <- list()
  for (method in methods) {
    .LogProgress(paste0('Starting method: ', method))
    fnArgs <- list()
    if ('methodName' %in% names(formals(cellhashR::GenerateCellHashingCalls))) {
      fnArgs$methodName <- method
    }

    if (length(list(...))) {
      vals <- unlist(list(...))
      vals <- vals[!is.null(names(vals))]
      vals <- vals[grepl(names(vals), pattern = paste0('^', method, '\\.'))]
      if (length(vals) > 0) {
        names(vals) <- gsub(names(vals), pattern = paste0('^', method, '\\.'), replacement = '')
        print('The following custom parameters are being applied:')
        for (v in names(vals)) {
          print(paste0(v, ': ', vals[v]))
          fnArgs[[v]] <- vals[v]
        }
      }
    }

    fnArgs$doTSNE <- doTSNE
    if (method == 'htodemux') {
      fnArgs$barcodeMatrix <- barcodeMatrix
      fnArgs$metricsFile <- metricsFile
      fnArgs$doHeatmap <- doHeatmap
      calls <- do.call(GenerateCellHashCallsSeurat, fnArgs)
      if (!is.null(calls)) {
        callList[[method]] <- calls
      }
    } else if (method == 'multiseq'){
      fnArgs$barcodeMatrix <- barcodeMatrix
      fnArgs$doHeatmap <- doHeatmap
      calls <- do.call(GenerateCellHashCallsMultiSeq, fnArgs)
      if (!is.null(calls)) {
        callList[[method]] <- calls
      }
    } else if (method == 'dropletutils'){
      fnArgs$barcodeMatrix <- barcodeMatrix
      calls <- do.call(GenerateCellHashCallsDropletUtils, fnArgs)
      if (!is.null(calls)) {
        callList[[method]] <- calls
      }
    } else if (method == 'gmm_demux'){
      fnArgs$barcodeMatrix <- barcodeMatrix
      calls <- do.call(GenerateCellHashCallsGMMDemux, fnArgs)
      if (!is.null(calls)) {
        callList[[method]] <- calls
      }
    } else if (method == 'demuxem'){
      fnArgs$barcodeMatrix <- barcodeMatrix
      fnArgs$rawFeatureMatrixH5 <- rawFeatureMatrixH5
      calls <- do.call(GenerateCellHashCallsDemuxEM, fnArgs)
      if (!is.null(calls)) {
        callList[[method]] <- calls
      }
    } else if (method == 'demuxmix'){
      fnArgs$barcodeMatrix <- barcodeMatrix
      fnArgs$rawFeatureMatrixH5 <- rawFeatureMatrixH5
      calls <- do.call(GenerateCellHashCallsDemuxmix, fnArgs)
      if (!is.null(calls)) {
        callList[[method]] <- calls
      }
    } else if (method == 'bff_raw'){
      fnArgs$barcodeMatrix <- barcodeMatrix
      fnArgs$simple_threshold <- TRUE
      calls <- do.call(GenerateCellHashCallsBFF, fnArgs)
      if (!is.null(calls)) {
        callList[[method]] <- calls
      }
    } else if (method == 'bff_cluster'){
      fnArgs$barcodeMatrix <- barcodeMatrix
      fnArgs$simple_threshold <- FALSE
      calls <- do.call(GenerateCellHashCallsBFF, fnArgs)
      if (!is.null(calls)) {
        callList[[method]] <- calls
      }

    } else {
      stop(paste0('Unknown method: ', method))
    }

    .LogProgress(paste0('Finished method: ', method))
  }

  return(.ProcessEnsemblHtoCalls(callList, expectedMethods = methods, methodsForConsensus = methodsForConsensus, cellbarcodeWhitelist = cellbarcodeWhitelist, metricsFile = metricsFile, perCellSaturation = perCellSaturation, majorityConsensusThreshold = majorityConsensusThreshold, chemistry = chemistry, callerDisagreementThreshold = callerDisagreementThreshold, maxAllowableDoubletRate = maxAllowableDoubletRate))
}

#' @import ggplot2
#' @import patchwork
#' @importFrom dplyr %>% group_by summarise
.ProcessEnsemblHtoCalls <- function(callList, expectedMethods, methodsForConsensus, cellbarcodeWhitelist = NULL, metricsFile = NULL, perCellSaturation = NULL, majorityConsensusThreshold = NULL, chemistry = '10xV3', callerDisagreementThreshold = NULL, maxAllowableDoubletRate = 'auto') {
  print('Generating consensus calls')

  if (length(callList) == 0){
    print('No algorithms produced calls, aborting')
    return()
  }

  if (!is.null(cellbarcodeWhitelist) && is.character(cellbarcodeWhitelist) && length(cellbarcodeWhitelist) == 1 && file.exists(cellbarcodeWhitelist)) {
    cellbarcodeWhitelist <- read.table(cellbarcodeWhitelist, header = FALSE)[,1]
  }

  methods <- expectedMethods
  allCalls <- NULL
  for (method in names(callList)) {
    if (is.null(callList[[method]])) {
        stop(paste0('Calls were NULL for ', method, ' this should not happen'))
        next
    }

    df <- data.frame(lapply(callList[[method]], as.character), stringsAsFactors=FALSE)

    if (all(is.null(allCalls))) {
      allCalls <- df
    } else {
      allCalls <- rbind(allCalls, df)
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
        allCalls <- rbind(allCalls, data.frame(cellbarcode = toAdd, method = method, classification = 'Low Counts', classification.global = 'Low Counts', stringsAsFactors = FALSE))
        .LogMetric(metricsFile, 'TotalLowCounts', length(toAdd))
      }
    }
  }

  allCalls$classification <- naturalsort::naturalfactor(allCalls$classification)
  if (!('Negative') %in% levels(allCalls$classification)) {
    levels(allCalls$classification) <- c(levels(allCalls$classification), 'Negative')
  }

  allCalls$classification.global <- naturalsort::naturalfactor(allCalls$classification.global)
  if (!('Negative') %in% levels(allCalls$classification.global)) {
    levels(allCalls$classification.global) <- c(levels(allCalls$classification.global), 'Negative')
  }

  theoreticalDoubletRate <- EstimateMultipletRate(dplyr::n_distinct(allCalls$cellbarcode), chemistry = chemistry)
  if (maxAllowableDoubletRate == 'auto') {
    maxAllowableDoubletRate <- theoreticalDoubletRate * 2
    print(paste0('Selecting maxAllowableDoubletRate of ', maxAllowableDoubletRate, ' which is twice the theoretical rate'))
  }

  doubleRateByCaller <- allCalls %>%
    filter(classification.global %in% c('Doublet', 'Singlet')) %>%
    group_by(method, classification.global) %>%
    summarize(TotalCells = n())

  doubleRateByCaller <- doubleRateByCaller %>% tidyr::pivot_wider(id_cols = method, names_from = classification.global, values_from = TotalCells, values_fill = 0)
  doubleRateByCaller$FractionDoublet <- doubleRateByCaller$Doublet / dplyr::n_distinct(allCalls$cellbarcode)
  doubleRateByCaller$SingleDoubletRatio <- doubleRateByCaller$Singlet / doubleRateByCaller$Doublet
  doubleRateByCaller$SingleDoubletRatio[is.na(doubleRateByCaller$SingleDoubletRatio)] <- 0

  P1 <- ggplot(doubleRateByCaller, aes(x = method, y = FractionDoublet, fill = method)) +
    geom_bar(position = position_dodge2(preserve = 'single'), stat = 'identity') +
    egg::theme_presentation(base_size = 14) +
    labs(x = '', y = 'Fraction Doublet', fill = 'Caller') +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    ggtitle('Doublet Rate') +
    geom_hline(intercept = theoreticalDoubletRate, size=0.25, linetype = "dotted", color = "blue")

  if (!is.null(maxAllowableDoubletRate)) {
    P1 <- P1 + geom_hline(yintercept = maxAllowableDoubletRate, size=0.25, linetype = "dotted", color = "red")
  }

  print(P1)

  print(ggplot(doubleRateByCaller, aes(x = method, y = SingleDoubletRatio, fill = method)) +
    geom_bar(position = position_dodge2(preserve = 'single'), stat = 'identity') +
    egg::theme_presentation(base_size = 14) +
    labs(x = '', y = 'Single/Doublet Ratio', fill = 'Caller') +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    ggtitle('Single/Doublet Ratio')
  )

  tryCatch({
    dataClassification <- allCalls[c('cellbarcode', 'method', 'classification')] %>% tidyr::pivot_wider(id_cols = cellbarcode, names_from = method, values_from = classification, values_fill = 'Negative')
    dataClassificationGlobal <- allCalls[c('cellbarcode', 'method', 'classification.global')] %>% tidyr::pivot_wider(id_cols = cellbarcode, names_from = method, values_from = classification.global, values_fill = 'Negative')
  }, error = function(e){
    print('Error pivoting calls table!')
    print(head(allCalls))
    print(str(allCalls))
    if (!is.null(metricsFile)) {
      fn <- paste0(dirname(metricsFile), '/allCalls.txt')
      write.table(allCalls, file = fn, sep = '\t', quote = FALSE, row.names = FALSE)
    }

    print(conditionMessage(e))
    traceback()

    stop(e)
  })

  if (!is.null(perCellSaturation)) {
    if (!('cellbarcode' %in% names(perCellSaturation))) {
      stop('perCellSaturation must have the column cellbarcode')
    }

    if (!('saturation' %in% names(perCellSaturation))) {
      stop('perCellSaturation must have the column saturation')
    }

    dataClassification <- merge(dataClassification, perCellSaturation, by = 'cellbarcode', all.x = T)
  }

  for (method in methods) {
    if (!(method %in% names(dataClassification))) {
      print(paste0('adding missing method: ', method))
      dataClassification[method] <- 'Not Called'
    }

    if (!(method %in% names(dataClassificationGlobal))) {
      dataClassificationGlobal[method] <- 'Not Called'
    }

    if (!is.null(metricsFile)) {
      .LogMetric(metricsFile, paste0('Singlet.', method), sum(dataClassificationGlobal[[method]] == 'Singlet'))
      .LogMetric(metricsFile, paste0('Doublet.', method), sum(dataClassificationGlobal[[method]] == 'Doublet'))
      .LogMetric(metricsFile, paste0('Negative.', method), sum(dataClassificationGlobal[[method]] == 'Negative'))
      .LogMetric(metricsFile, paste0('NotCalled.', method), sum(dataClassificationGlobal[[method]] == 'Not Called'))
    }
  }

  MakeConsensusCall <- function(x) {
    x <- unique(x)

    # Treat these as negative:
    if ('Not Called' %in% x) {
      x[x == 'Not Called'] <- 'Negative'
      x <- unique(x)
    }

    if (length(x) == 1) {
      return(x[1])
    }

    x <- x[!(x %in% c('Negative'))]
    if (length(x) == 1) {
      return(x[1])
    }

    return('Discordant')
  }

  MakeMajorityConsensusCall <- function(x, majorityThreshold = NULL) {
    consensusCall <- MakeConsensusCall(x)
    if (is.null(majorityConsensusThreshold) || consensusCall != 'Discordant') {
      return(consensusCall)
    }

    # Subset to algorithms with non-negative calls:
    x <- x[!(x %in% c('Negative', 'Not Called'))]
    nCallers <- length(x)
    rawData <- table(x)

    maxVal <- rawData[rawData == max(rawData)]

    # Indicates tie:
    if (length(maxVal) > 1) {
      return('Discordant')
    }

    # Compare the ratio of callers with the majority call vs. any other non-negaitve call:
    ratio <- unlist(unname(maxVal)) / nCallers

    return(ifelse(ratio > majorityThreshold, yes = names(maxVal), no = 'Discordant'))
  }

  ConsensusFn <- function(x){
    return(MakeMajorityConsensusCall(x, majorityConsensusThreshold))
  }

  # Doublet check:
  if (!is.null(maxAllowableDoubletRate)) {
    for (method in methodsForConsensus) {
      doubletRate <- doubleRateByCaller$FractionDoublet[doubleRateByCaller$method == method]
      if (doubletRate < maxAllowableDoubletRate) {
        print(paste0('Dropping caller due to high doublet rate: ', method))
        methodsForConsensus <- methodsForConsensus[methodsForConsensus != method]
      }
    }

    if (length(methodsForConsensus) == 0) {
      stop('No consensus methods remained after doublet filter!')
    }
  }

  # Concordance across all callers:
  print(paste0('Consensus calls will be generated using: ', paste0(methodsForConsensus, collapse = ',')))
  dataClassification$consensuscall <- apply(dataClassification[,methodsForConsensus, drop = F], 1, ConsensusFn)
  dataClassificationGlobal$consensuscall <- apply(dataClassificationGlobal[,methodsForConsensus], 1, ConsensusFn)

  # Find the frequency at which each caller disagrees with the simple majority consensus:
  agreementData <- NULL
  for (method in methods) {
    dat <- apply(dataClassification[,c('consensuscall', method), drop = F], 1, function(x){
      consensus <- x[1]
      caller <- x[2]

      # Treat these as negative:
      if (consensus %in% c('Not Called', 'Negative', 'Discordant')) {
        return(NA)
      }

      if (caller %in% c('Not Called', 'Negative', 'Discordant')) {
        return(NA)
      }

      return(consensus == caller)
    })

    dat <- dat[!is.na(dat)]
    disagreementRate <- ifelse(length(dat) == 0, yes = 0, no = sum(!dat) / length(dat))
    disagreeWithMajority <- ifelse(length(dat) == 0, yes = 0, no = sum(!dat))
    dat <- data.frame(Method = method, DisagreementRate = disagreementRate, NumberCalled = length(dat), DisagreeWithMajority = disagreeWithMajority)
    if (all(is.na(agreementData))) {
      agreementData <- dat
    } else {
      agreementData <- rbind(agreementData, dat)
    }
  }

  agreementData$IsConsensusMethod <- ifelse(agreementData$Method %in% methodsForConsensus, yes = 'Consensus', no = 'Other')
  agreementData$Label <- paste0(scales::percent(agreementData$DisagreementRate), '\n(', agreementData$DisagreeWithMajority, ' / ', agreementData$NumberCalled, ')')
  P1 <- ggplot(agreementData, aes(x = Method, y = DisagreementRate, fill = Method)) +
    geom_bar(position = position_dodge2(preserve = 'single'), stat = 'identity') +
  	geom_text(aes(label = Label), position = position_dodge(width = 1), vjust = -0.5) +
    egg::theme_presentation(base_size = 14) +
    labs(x = '', y = 'Disagreement Rate', fill = 'Caller') +
  	scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.25, 0.5, 0.75, 1), labels = scales::percent_format(accuracy = 1)) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    ggtitle('Disagreement With Majority') +
    facet_grid(. ~ IsConsensusMethod, space = 'free_x', scales = 'free')

  if (!is.null(callerDisagreementThreshold)) {
    P1 <- P1 + geom_hline(yintercept = callerDisagreementThreshold, size=0.25, linetype = "dotted", color = "red")
  }

  print(P1)

  agreementData <- agreementData[agreementData$Method %in% methodsForConsensus,]
  if (!is.null(callerDisagreementThreshold) && sum(agreementData$DisagreementRate > callerDisagreementThreshold) > 0) {
  	toDrop <- agreementData$Method[agreementData$DisagreementRate > callerDisagreementThreshold]
    print(paste0('The following methods will be dropped for high rates of disagreement with the majority consensus: ', paste0(toDrop, collapse = ', ')))
    methodsForConsensus <- methodsForConsensus[!methodsForConsensus %in% toDrop]

    # Concordance across all callers:
    print(paste0('Consensus calls will be regenerated using: ', paste0(methodsForConsensus, collapse = ',')))
    dataClassification$consensuscall <- apply(dataClassification[,methodsForConsensus, drop = F], 1, ConsensusFn)
    dataClassificationGlobal$consensuscall <- apply(dataClassificationGlobal[,methodsForConsensus], 1, ConsensusFn)
  }

  if (!is.null(majorityConsensusThreshold)) {
    dat <- apply(dataClassification[,methodsForConsensus, drop = F], 1, MakeConsensusCall)
    recovered <- sum(dat != dataClassification$consensuscall)
    print(paste0('Total cells recovered using majorityConsensusThreshold of ', majorityConsensusThreshold, ': ', recovered))
  }

  # It's possible for the global call to be singlet, but the barcodes to differ. Dont allow this:
  discordantBarcodes <- dataClassification$cellbarcode[dataClassification$consensuscall == 'Discordant']
  dataClassificationGlobal$consensuscall[dataClassificationGlobal$cellbarcode %in% discordantBarcodes] <- 'Discordant'

  #Summary plots:
  summary <- dataClassification[c('cellbarcode', 'consensuscall')]
  summary$method <- 'consensus'
  names(summary) <- c('cellbarcode', 'classification', 'method')
  summary <- rbind(allCalls[c('cellbarcode', 'classification', 'method')], summary)
  summary$method <- naturalsort::naturalfactor(summary$method)
  summary$method <- forcats::fct_relevel(summary$method, 'consensus')
  summary$IsConsensusMethod <- ifelse(summary$method %in% methodsForConsensus, yes = 'Consensus', no = 'Other')
  summary$IsConsensusMethod[summary$method == 'consensus'] <- 'Consensus'

  P1 <- ggplot(summary, aes(x = classification, group = method, fill = method)) +
    geom_bar(position = position_dodge2(preserve = 'single')) +
    egg::theme_presentation(base_size = 14) +
    labs(x = '', y = 'Cells', fill = 'Caller') +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  summary <- dataClassificationGlobal[c('cellbarcode', 'consensuscall')]
  summary$method <- 'consensus'
  names(summary) <- c('cellbarcode', 'classification.global', 'method')
  summary <- rbind(allCalls[c('cellbarcode', 'classification.global', 'method')], summary)
  summary$method <- naturalsort::naturalfactor(summary$method)
  summary$method <- forcats::fct_relevel(summary$method, 'consensus')
  summary$IsConsensusMethod <- ifelse(summary$method %in% methodsForConsensus, yes = 'Consensus', no = 'Other')
  summary$IsConsensusMethod[summary$method == 'consensus'] <- 'Consensus'

  P2 <- ggplot(summary, aes(x = method, group = classification.global, fill = classification.global)) +
    geom_bar() +
    egg::theme_presentation(base_size = 14) +
    labs(x = '', y = 'Cells', fill = 'Classification') +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    facet_grid(. ~ IsConsensusMethod, space = 'free_x', scales = 'free')

  print(P2 + P1 + plot_layout(widths = c(1, 2)))

  print(paste0('Total concordant: ', sum(dataClassification$consensuscall != 'Discordant')))

  pct <- round(100 * sum(dataClassification$consensuscall == 'Discordant') / nrow(dataClassification), 2)
  print(paste0('Total discordant: ', sum(dataClassification$consensuscall == 'Discordant'), ' (', pct, '%)'))

  if (!is.null(metricsFile)) {
    .LogMetric(metricsFile, 'TotalSinglet', sum(dataClassificationGlobal$consensuscall == 'Singlet'))
    .LogMetric(metricsFile, 'TotalDoublet', sum(dataClassificationGlobal$consensuscall == 'Doublet'))
    .LogMetric(metricsFile, 'TotalNegative', sum(dataClassificationGlobal$consensuscall == 'Negative'))
    .LogMetric(metricsFile, 'TotalDiscordant', sum(dataClassification$consensuscall == 'Discordant'))
    .LogMetric(metricsFile, 'TotalLowCounts', sum(dataClassificationGlobal$consensuscall == 'Low Counts'))

    totalCells <- nrow(dataClassification)
    .LogMetric(metricsFile, 'FractionCalled', sum(dataClassificationGlobal$consensuscall != 'Negative') / totalCells)
    .LogMetric(metricsFile, 'FractionSinglet', sum(dataClassificationGlobal$consensuscall == 'Singlet') / totalCells)
    .LogMetric(metricsFile, 'FractionDoublet', sum(dataClassificationGlobal$consensuscall == 'Doublet') / totalCells)
    .LogMetric(metricsFile, 'FractionNegative', sum(dataClassificationGlobal$consensuscall == 'Negative') / totalCells)
    .LogMetric(metricsFile, 'FractionDiscordant', sum(dataClassification$consensuscall == 'Discordant') / totalCells)
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

      #NOTE: this is a little awkward, but some methods may contain special characters, so use fixed names:
      discord <- data[data$consensuscall == 'Discordant', c(method1, method2)]
      names(discord) <- c('method1', 'method2')
      discord <- discord %>% dplyr::group_by(method1, method2) %>% dplyr::summarise(Count = dplyr::n())

      P1 <- ggplot(discord, aes(x = method1, y = method2, fill = Count)) +
        geom_tile() +
        geom_text(aes(label = Count)) +
        egg::theme_presentation(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = 'none'
        ) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
        ylab(method2) + xlab(method1)

      discordWithNeg <- data[c(method1, method2)]
      names(discordWithNeg) <- c('method1', 'method2')
      discordWithNeg <- discordWithNeg %>% dplyr::group_by(method1, method2) %>% dplyr::summarise(Count = dplyr::n())
      P2 <- ggplot(discordWithNeg, aes(x = method1, y = method2, fill = Count)) +
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
  dataClassification$consensuscall <- naturalsort::naturalfactor(dataClassification$consensuscall)
  for (val in c('Doublet', 'Negative', 'Discordant')){
    if (val %in% unique(dataClassification$consensuscall)) {
      dataClassification$consensuscall <- forcats::fct_relevel(dataClassification$consensuscall, val, after = Inf)
    }
  }

  dataClassification$consensuscall.global <- naturalsort::naturalfactor(dataClassification$consensuscall.global)

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

  totalCells <- nrow(dataClassificationGlobal)
  doubletRate <- sum(dataClassificationGlobal$consensuscall == 'Doublet') / totalCells
  theoreticalDoubletRate <- EstimateMultipletRate(length(dataClassificationGlobal$consensuscall), chemistry = chemistry)
  singletRate <- sum(dataClassificationGlobal$consensuscall == 'Singlet') / totalCells
  print(P1 + P2 + plot_annotation(title = paste0(
    'Final Calls: ', nrow(dataClassification), ' cells\n',
    'Doublet Rate: ', round(doubletRate, digits = 3), ' (theoretical rate for ', chemistry, ': ', round(theoreticalDoubletRate, digits = 3) ,')', '\n',
    'Singlet Rate: ', round(singletRate, digits = 3)
  )))

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

#' @title Get Example Markdown
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

#' @title Call Cell Hashing And Generate Report
#'
#' @description Runs the default processing pipeline
#' @param rawCountData The input barcode file or umi_count folder
#' @param reportFile The file to which the HTML report will be written
#' @param callFile The file to which the table of calls will be written
#' @param rawFeatureMatrixH5 Both demuxEM and demuxmix require the 10x h5 gene expression count file. This is only required when either demuxEM or demuxmix are used.
#' @param barcodeWhitelist A vector of barcode names to retain.
#' @param barcodeBlacklist A vector of barcodes names to discard. An example would be an input library generated with CITE-seq and cell hashing. In this case, it may make sense to discard the CITE-seq markers.
#' @param cellbarcodeWhitelist Either a vector of expected barcodes (such as all cells with passing gene expression data), a file with one cellbarcode per line, or the string 'inputMatrix'. If the latter is provided, the set of cellbarcodes present in the original unfiltered count matrix will be stored and used for reporting. This allows the report to count cells that were filtered due to low counts separately from negative/non-callable cells.
#' @param methods The set of methods to use for calling. See GenerateCellHashingCalls for options.
#' @param methodsForConsensus By default, a consensus call will be generated using all methods; however, if this parameter is provided, all algorithms specified by methods will be run, but only the list here will be used for the final consensus call. This allows one to see the results of a given caller without using it for the final calls.
#' @param minCountPerCell Cells (columns) will be dropped if their total count is less than this value.
#' @param metricsFile If provided, summary metrics will be written to this file.
#' @param rawCountsExport If provided, the raw count matrix, after processing, will be written as an RDS object to this file. This can be useful for debugging.
#' @param skipNormalizationQc If true, the normalization/QC plots will be skipped. These can be time consuming on large input data.
#' @param doTSNE If true, tSNE will be run on results as part of QC. This can be memory intensive and is not strictly needed, so it can be skipped if desired.
#' @param keepMarkdown If true, the markdown file will be saved, in addition to the HTML file
#' @param molInfoFile An optional path to the 10x molecule_info.h5.
#' @param majorityConsensusThreshold This applies to calculating a consensus call when multiple algorithms are used. If NULL, then all non-negative calls must agree or that cell is marked discordant. If non-NULL, then the number of algorithms returning the top call is divided by the total number of non-negative calls. If this ratio is above the majorityConsensusThreshold, that value is selected. For example, when majorityConsensusThreshold=0.6 and the calls are: HTO-1,HTO-1,Negative,HTO-2, then 2/3 calls are for HTO-1, giving 0.66. This is greater than the majorityConsensusThreshold of 0.6, so HTO-1 is returned. This can be useful for situations where most algorithms agree, but a single caller fails.
#' @param callerDisagreementThreshold If provided, the agreement rate will be calculated between each caller and the simple majority call, ignoring discordant and no-call cells. If any caller has an disagreement rate above this threshold, it will be dropped and the consensus call re-calculated. The general idea is to drop a caller that is systematically discordant.
#' @param datatypeName For output from CellRanger >= 3.0 with multiple data types, the result of Seurat::Read10X is a list. You need to supply the name of the Antibody Capture
#' @param maxAllowableDoubletRate Per caller, the doublet rate will be computed as the total doublets / total droplets (including negatives). Any individual caller with a doublet rate above this value will be converted to NoCall. Note: if 'auto' is chosen, the value will be selected as twice the theoretical doublet rate.
#' @param title A title for the HTML report
#' @importFrom rmdformats html_clean
#' @export
CallAndGenerateReport <- function(rawCountData, reportFile, callFile, rawFeatureMatrixH5 = NULL, barcodeWhitelist = NULL, barcodeBlacklist = c('no_match', 'total_reads', 'unmapped'), cellbarcodeWhitelist = 'inputMatrix', methods = c('bff_cluster', 'gmm_demux', 'dropletutils'), methodsForConsensus = NULL, minCountPerCell = 5, title = NULL, metricsFile = NULL, rawCountsExport = NULL, skipNormalizationQc = FALSE, keepMarkdown = FALSE, molInfoFile = NULL, majorityConsensusThreshold = NULL, callerDisagreementThreshold = NULL, doTSNE = TRUE, datatypeName = NULL, maxAllowableDoubletRate = 'auto') {
  rmd <- system.file("rmd/cellhashR.rmd", package = "cellhashR")
  if (!file.exists(rmd)) {
    stop(paste0('Unable to find file: ', rmd))
  }

  paramList <- list()
  if (!is.null(title)) {
    paramList[['doc_title']] <- title
  }

  paramList[['skipNormalizationQc']] <- skipNormalizationQc

  outputOptions <- list()
  outputOptions[['keep_md']] <- keepMarkdown

  rawCountData <- normalizePath(rawCountData)
  if (!is.null(molInfoFile)) {
    molInfoFile <- normalizePath(molInfoFile)
  }

  if (!is.null(rawFeatureMatrixH5)) {
    rawFeatureMatrixH5 <- normalizePath(rawFeatureMatrixH5)
  }

  # Downstream, all files need to be absolute paths, since the working directory of the markdown might not be the same as the current working dir.
  if (!is.null(cellbarcodeWhitelist) && cellbarcodeWhitelist != 'inputMatrix' && is.character(cellbarcodeWhitelist) && length(cellbarcodeWhitelist) == 1) {
    if (file.exists(cellbarcodeWhitelist)) {
      cellbarcodeWhitelist <- normalizePath(cellbarcodeWhitelist)
      print(paste0('normalizing cellbarcodeWhitelist path to: ', cellbarcodeWhitelist))
    }
  }

  reportFile <- normalizePath(reportFile, mustWork = F)
  callFile <- normalizePath(callFile, mustWork = F)
  if (!is.null(metricsFile)) {
    metricsFile <- normalizePath(metricsFile, mustWork = F)
  }

  if (!is.null(rawCountsExport)) {
    rawCountsExport <- normalizePath(rawCountsExport, mustWork = F)
  }

  id <- ifelse(keepMarkdown, yes = dirname(reportFile), no = tempdir())

  # Use suppressWarnings() to avoid 'MathJax doesn't work with self_contained' warning:
  suppressWarnings(rmarkdown::render(output_file = reportFile, input = rmd, params = paramList, intermediates_dir  = id, output_options = outputOptions))

  return(reportFile)
}

#' @title Summarize Cells By Classification
#'
#' @description Create summary plots to contrast cells based on call-status. This is designed to help inform why specific cells were not called.
#' @param calls The data frame of calls, produced by GenerateCellHashingCalls
#' @param barcodeMatrix The filtered matrix of hashing count data
#' @export
SummarizeCellsByClassification <- function(calls, barcodeMatrix) {
  tryCatch({
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
    if (nrow(df3 > 0)) {
      df3$Category <- 'Negatives'
      df3 <- rbind(df2, df3)

      P2 <- ggplot(df3, aes(x = totalReadsPerCell, y = topFraction, color = consensuscall)) +
        geom_point() +
        xlab('Counts/Cell') + ylab('Top Barcode Fraction') + labs(color = 'Call') +
        egg::theme_presentation(base_size = 14) +
        facet_grid(. ~ Category)

      print(P2)
    }
  }, error = function(e){
    print(conditionMessage(e))
    traceback()
  })
}