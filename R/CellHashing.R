utils::globalVariables(
  names = c('p_val_adj', 'avg_logFC', 'cluster'),
  package = 'cellhashR',
  add = TRUE
)


utils::globalVariables(
  names = c('HTO_Classification', 'HTO_classification.global', 'HTO', 'Count'),
  package = 'cellhashR',
  add = TRUE
)

#' @title GenerateCellHashCallsSeurat
#'
#' @description Generates final cell hashing calls using Seurat3 HTODemux
#' @return A data frame of results
GenerateCellHashCallsSeurat <- function(barcodeMatrix, positive.quantile = 0.99, attemptRecovery = F, minCellsForRecovery = 500) {
  seuratObj <- CreateSeuratObject(barcodeMatrix, assay = 'HTO')

  tryCatch({
    seuratObj <- DoHtoDemux(seuratObj, positive.quantile = positive.quantile)
    dt <- data.frame(Barcode = as.factor(colnames(seuratObj)), HTO_classification = seuratObj$hash.ID, HTO_classification.all = seuratObj$HTO_classification, HTO_classification.global = seuratObj$HTO_classification.global, key = c('Barcode'), stringsAsFactors = F)

    #attempt recovery:
    if (attemptRecovery && length(seuratObj$HTO_classification.global == 'Negative') > minCellsForRecovery) {
      print('Attempting recovery of negatives')
      seuratObj2 <- subset(seuratObj, subset = HTO_classification.global == 'Negative')
      seuratObj2 <- CreateSeuratObject(seuratObj2@assays$HTO@counts, assay = 'HTO')
      print(paste0('Initial negative cells : ', ncol(seuratObj2)))

      tryCatch({
        seuratObj2 <- DoHtoDemux(seuratObj2, positive.quantile = positive.quantile, label = 'Seurat HTODemux (2nd Round)')
        seuratObj2 <- subset(seuratObj2, subset = HTO_classification.global == 'Singlet')
        print(paste0('Total cells rescued by second HTODemux call: ', ncol(seuratObj2)))
        if (ncol(seuratObj2) > 0) {
          dt2 <- data.frame(Barcode = as.factor(colnames(seuratObj2)), HTO_classification = seuratObj2$hash.ID, HTO_classification.all = seuratObj2$HTO_classification, HTO_classification.global = seuratObj2$HTO_classification.global, key = c('Barcode'), stringsAsFactors = F)
          dt[dt$Barcode %in% dt2$Barcode]$HTO_classification <- dt2$HTO_classification
          dt[dt$Barcode %in% dt2$Barcode]$HTO_classification.all <- dt2$HTO_classification.all
          dt[dt$Barcode %in% dt2$Barcode]$HTO_classification.global <- dt2$HTO_classification.global
        }
      }, error = function(e){
        print('Error generating second round of seurat calls, using first round')
        return(dt)
      })
    } else if (attemptRecovery) {
      print(paste0('Attempt recovery was selected but will not be performed because there were insufficient negative cells after first round: ', length(seuratObj$HTO_classification.global == 'Negative')))
    }

    return(dt)
  }, error = function(e){
    print(e)
    print('Error generating seurat calls, aborting')
    saveRDS(seuratObj, file = './seuratHashingError.rds')
    return(NA)
  })
}

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

#' @title A Title
#'
#' @description A description
#' @param seuratObj, A Seurat object.
#' @importFrom cluster clara
#' @importFrom Matrix t
#' @return A modified Seurat object.
DebugDemux <- function(seuratObj, assay = 'HTO', reportKmeans = FALSE) {
  print('Debugging information for Seurat HTODemux:')
  data <- GetAssayData(object = seuratObj, assay = assay)
  ncenters <- (nrow(x = data) + 1)

  init.clusters <- clara(
    x = t(x = data),
    k = ncenters,
    samples = 100
  )
  Idents(object = seuratObj, cells = names(x = init.clusters$clustering), drop = TRUE) <- init.clusters$clustering

  # Calculate tSNE embeddings with a distance matrix
	tryCatch({
		perplexity <- .InferPerplexityFromSeuratObj(seuratObj, 100)
		seuratObj[['hto_tsne']] <- RunTSNE(dist(t(data)), assay = assay, perplexity = perplexity)
		P <- DimPlot(seuratObj, reduction = 'hto_tsne', label = TRUE)
		P <- P + ggtitle('Clusters: clara')
		print(P)
	}, error = function(e){
		print(e)
		print('Error generating tSNE, skipping')
	})

  average.expression <- AverageExpression(
    object = seuratObj,
    assays = c(assay),
    verbose = FALSE
  )[[assay]]

  print(kableExtra::kbl(average.expression, label = 'clara') %>% kableExtra::kable_styling())

  if (reportKmeans) {
    print('kmeans:')
    init.clusters <- kmeans(
      x = t(x = data),
      centers = ncenters,
      nstart = 100
    )
    Idents(object = seuratObj, cells = names(x = init.clusters$cluster), drop = TRUE) <- init.clusters$cluster

    # Calculate tSNE embeddings with a distance matrix
    P <- DimPlot(seuratObj, label = TRUE)
    P <- P + ggtitle('Clusters: kmeans')
    print(P)

    average.expression <- AverageExpression(
      object = seuratObj,
      assays = c(assay),
      verbose = FALSE
    )[[assay]]

    print(kableExtra::kbl(average.expression, label = 'kmeans') %>% kableExtra::kable_styling())
  }
}

#' @title A Title
#'
#' @description A description
#' @param seuratObj, A Seurat object.
#' @return A modified Seurat object.
DoHtoDemux <- function(seuratObj, positive.quantile = 0.99, label = 'Seurat HTODemux', plotDist = FALSE) {
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  seuratObj <- NormalizeData(seuratObj, assay = "HTO", normalization.method = "CLR", verbose = FALSE)

  DebugDemux(seuratObj)
  PlotHtoCountData(seuratObj, label = 'Seurat HTO Demux')

  seuratObj <- HTODemux2(seuratObj, positive.quantile =  positive.quantile, plotDist = plotDist)
  seuratObj$hash.ID <- naturalsort::naturalfactor(as.character(seuratObj$hash.ID))

  HtoSummary(seuratObj, label = label, htoClassificationField = 'hash.ID', globalClassificationField = 'HTO_classification.global')

  return(seuratObj)
}

#' @title A Title
#'
#' @description A description
#' @return A modified Seurat object.
GenerateCellHashCallsMultiSeq <- function(barcodeMatrix, method = 'seurat') {
  seuratObj <- CreateSeuratObject(barcodeMatrix, assay = 'HTO')

  tryCatch({
    if (method == 'seurat') {
      seuratObj <- DoMULTIseqDemux(seuratObj)
    } else {
      stop(paste0('Unknown method: ', method))
    }

    return(data.frame(Barcode = as.factor(colnames(seuratObj)), HTO_classification = seuratObj$MULTI_ID, HTO_classification.global = seuratObj$MULTI_classification.global, key = c('Barcode')))
  }, error = function(e){
    print(e)
    print('Error generating multiseq calls, aborting')
    saveRDS(seuratObj, file = './multiseqHashingError.rds')
    return(NA)
  })
}

#' @title GenerateCellHashingCalls
#'
#' @description A description
#' @return A data table of results.
#' @export
GenerateCellHashingCalls <- function(barcodeMatrix, positive.quantile = 0.99, attemptRecovery = FALSE, useSeurat = TRUE, useMultiSeq = TRUE, outFile = 'combinedHtoCalls.txt', allCallsOutFile = NA) {
  sc <- NA
  if (useSeurat) {
    sc <- GenerateCellHashCallsSeurat(barcodeMatrix, positive.quantile = positive.quantile, attemptRecovery = attemptRecovery )
  }

  mc <- NA
  if (useMultiSeq) {
    mc <- GenerateCellHashCallsMultiSeq(barcodeMatrix)
  }

  dt <- ProcessEnsemblHtoCalls(mc, sc, barcodeMatrix, outFile = outFile, allCallsOutFile = allCallsOutFile)

  return(dt)
}

#' @title Perform HTO classification using Seurat's implementation of MULTI-seq classification
#'
#' @description Perform HTO classification using Seurat's implementation of MULTI-seq classification
#' @param seuratObj, A Seurat object.
#' @return A modified Seurat object.
DoMULTIseqDemux <- function(seuratObj, assay = 'HTO', autoThresh = TRUE, quantile = NULL, maxiter = 20, qrange = seq(from = 0.2, to = 0.95, by = 0.05)) {
  ## Normalize Data: Log2 Transform, mean-center
  counts <- GetAssayData(seuratObj, assay = assay, slot = 'counts')
  log2Scaled <- as.data.frame(log2(Matrix::t(counts)))
  for (i in 1:ncol(log2Scaled)) {
    ind <- which(is.finite(log2Scaled[,i]) == FALSE)
    log2Scaled[ind,i] <- 0
    log2Scaled[,i] <- log2Scaled[,i]-mean(log2Scaled[,i])
  }
  seuratObjMS <- CreateSeuratObject(counts, assay = 'MultiSeq')
  seuratObjMS[['MultiSeq']]@data <- Matrix::t(as.matrix(log2Scaled))

  PlotHtoCountData(seuratObjMS, label = 'MultiSeq', assay = 'MultiSeq')

  seuratObjMS <- MULTIseqDemux2(seuratObjMS, assay = "MultiSeq", quantile = quantile, verbose = TRUE, autoThresh = autoThresh, maxiter = maxiter, qrange = qrange)

  seuratObjMS$MULTI_classification.global <- as.character(seuratObjMS$MULTI_ID)
  seuratObjMS$MULTI_classification.global[!(seuratObjMS$MULTI_ID %in% c('Negative', 'Doublet'))] <- 'Singlet'
  seuratObjMS$MULTI_classification.global <- as.factor(seuratObjMS$MULTI_classification.global)

  seuratObjMS$MULTI_ID <- naturalsort::naturalfactor(as.character(seuratObjMS$MULTI_ID))

  HtoSummary(seuratObjMS, label = 'MULTI-SEQ', htoClassificationField = 'MULTI_ID', globalClassificationField = 'MULTI_classification.global', assay = 'MultiSeq')

  seuratObj$MULTI_ID <- as.character(seuratObjMS$MULTI_ID)
  seuratObj$MULTI_classification.global <- seuratObjMS$MULTI_classification.global

  return(seuratObj)
}

PlotHtoCountData <- function(seuratObj, label, assay = 'HTO') {
  #Plot raw data by HTO:
  data <- GetAssayData(seuratObj, assay = assay, slot = 'counts')
  df2 <- as.data.frame(data)
  df2$HTO <- row.names(data)
  df2 <- tidyr::gather(df2, key = 'CellBarcode', value = 'Count', -HTO)
  df2$HTO <- simplifyHtoNames(as.character(df2$HTO))
  df2$HTO <- naturalsort::naturalfactor(df2$HTO)

  df2$Count <- log10(df2$Count + 0.5)
  print(ggplot(df2, aes(y = Count, x = HTO, fill = HTO)) +
    geom_boxplot() +
    ggtitle(paste0(label, ': Raw counts by HTO (log10)'))
  )

  #Plot normalized data by HTO:
  data <- GetAssayData(seuratObj, assay = assay, slot = 'data')
  df2 <- as.data.frame(data)
  df2$HTO <- row.names(data)
  df2 <- tidyr::gather(df2, key = 'CellBarcode', value = 'Count', -HTO)
  df2$HTO <- simplifyHtoNames(as.character(df2$HTO))
  df2$HTO <- naturalsort::naturalfactor(df2$HTO)

  print(ggplot(df2, aes(y = Count, x = HTO, fill = HTO)) +
    geom_boxplot() +
    ggtitle(paste0(label, ': Normalized data by HTO'))
  )
}

#' @title A Title
#'
#' @description A description
#' @param seuratObj, A Seurat object.
#' @return A modified Seurat object.
HtoSummary <- function(seuratObj, htoClassificationField, globalClassificationField, label, doHeatmap = T, doTSNE = T, assay = 'HTO') {
  #report outcome
  print(table(seuratObj[[htoClassificationField]]))
  print(table(seuratObj[[globalClassificationField]]))

  # Group cells based on the max HTO signal
  seuratObj_hashtag <- seuratObj
  Idents(seuratObj_hashtag) <- htoClassificationField
  htos <- rownames(GetAssayData(seuratObj_hashtag, assay = assay))
  for (hto in naturalsort::naturalsort(htos)){
    print(VlnPlot(seuratObj_hashtag, features = c(hto), assay = assay, ncol = 1, log = F) + ggtitle(paste0(label, ": ", hto)))
  }

  if (doTSNE) {
    perplexity <- .InferPerplexityFromSeuratObj(seuratObj, 100)
    tryCatch({
      seuratObj[['hto_tsne']] <- RunTSNE(dist(Matrix::t(GetAssayData(seuratObj, slot = "data", assay = assay))), assay = assay, perplexity = perplexity)
      print(DimPlot(seuratObj, reduction = 'hto_tsne', group.by = htoClassificationField, label = FALSE) + ggtitle(label))
      print(DimPlot(seuratObj, reduction = 'hto_tsne', group.by = globalClassificationField, label = FALSE) + ggtitle(label))
    }, error = function(e){
      print(e)
      print('Error generating tSNE, skipping')
    })
  }

  if (doHeatmap) {
    print(HTOHeatmap(seuratObj, assay = assay, classification = htoClassificationField, global.classification = globalClassificationField, ncells = min(3000, ncol(seuratObj)), singlet.names = NULL) + ggtitle(label))
  }
}

utils::globalVariables(
  names = c('HTO_classification.global.MultiSeq', 'HTO_classification.global.Seurat', 'n', 'HTO_classification.MultiSeq', 'HTO_classification.Seurat'),
  package = 'cellhashR',
  add = TRUE
)

#' @title ProcessEnsemblHtoCalls
#'
#' @description A description
#' @return A modified Seurat object.
#' @export
#' @import ggplot2
#' @param mc Multiseq calls dataframe
#' @param sc Seurat calls dataframe
#' @param barcodeMatrix The barcode count matrix
#' @param outFile The output TSV file
#' @param allCallsOutFile If provided, a more detailed output will be written here
#' @importFrom dplyr %>% group_by summarise
ProcessEnsemblHtoCalls <- function(mc, sc, barcodeMatrix,
                                   outFile = 'combinedHtoCalls.txt',
                                   allCallsOutFile = NA) {

  if (all(is.na(sc)) && all(is.na(mc))){
    print('MULTI-seq and Seurat failed to produce calls, aborting')
    return()
  }

  if (all(is.na(sc))){
    print('No calls for seurat found')
    dt <- data.frame(CellBarcode = mc$Barcode, HTO = mc$HTO_classification, HTO_Classification = mc$HTO_classification.global, key = 'CellBarcode', Seurat = c(F), MultiSeq = c(T))
    dt <- PrintFinalSummary(dt, barcodeMatrix)
    write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)

    return(dt)
  }

  if (all(is.na(mc))){
    print('No calls for MULTI-seq found')
    dt <- data.frame(CellBarcode = sc$Barcode, HTO = sc$HTO_classification, HTO_Classification = sc$HTO_classification.global, key = 'CellBarcode', Seurat = c(T), MultiSeq = c(F))
    dt <- PrintFinalSummary(dt, barcodeMatrix)
    write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)

    return(dt)
  }

  mc$Barcode <- as.character(mc$Barcode)
  sc$Barcode <- as.character(sc$Barcode)
  merged <- merge(mc, sc, all = T, by = 'Barcode', suffixes = c('.MultiSeq', '.Seurat'))

  merged$HTO_classification.MultiSeq <- as.character(merged$HTO_classification.MultiSeq)
  merged$HTO_classification.MultiSeq[is.na(merged$HTO_classification.MultiSeq)] <- 'Negative'
  merged$HTO_classification.MultiSeq <- naturalsort::naturalfactor(merged$HTO_classification.MultiSeq)

  merged$HTO_classification.Seurat <- as.character(merged$HTO_classification.Seurat)
  merged$HTO_classification.Seurat[is.na(merged$HTO_classification.Seurat)] <- 'Negative'
  merged$HTO_classification.Seurat <- naturalsort::naturalfactor(merged$HTO_classification.Seurat)

  merged$HTO_classification.global.MultiSeq <- as.character(merged$HTO_classification.global.MultiSeq)
  merged$HTO_classification.global.MultiSeq[is.na(merged$HTO_classification.global.MultiSeq)] <- 'Negative'
  merged$HTO_classification.global.MultiSeq <- naturalsort::naturalfactor(merged$HTO_classification.global.MultiSeq)

  merged$HTO_classification.global.Seurat <- as.character(merged$HTO_classification.global.Seurat)
  merged$HTO_classification.global.Seurat[is.na(merged$HTO_classification.global.Seurat)] <- 'Negative'
  merged$HTO_classification.global.Seurat <- naturalsort::naturalfactor(merged$HTO_classification.global.Seurat)

  merged$HasSeuratCall <- merged$HTO_classification.Seurat != 'Negative'
  merged$HasMultiSeqCall <- merged$HTO_classification.MultiSeq != 'Negative'

  #dont count situations where one side is negative as discordant
  merged$Concordant <- as.character(merged$HTO_classification.MultiSeq) == as.character(merged$HTO_classification.Seurat)
  merged$Concordant[!merged$Concordant & (merged$HTO_classification.MultiSeq == 'Negative' | merged$HTO_classification.Seurat == 'Negative')] <- TRUE

  merged$GlobalConcordant <- as.character(merged$HTO_classification.global.MultiSeq) == as.character(merged$HTO_classification.global.Seurat)
  merged$GlobalConcordant[!merged$GlobalConcordant & (merged$HTO_classification.global.MultiSeq == 'Negative' | merged$HTO_classification.global.Seurat == 'Negative')] <- TRUE

  print(paste0('Total concordant: ', sum(merged$Concordant)))
  print(paste0('Total discordant (HTO call): ', sum(!merged$Concordant)))
  print(paste0('Total discordant (HTO classification): ', sum(!merged$GlobalConcordant)))

  tbl <- data.frame(table(Concordant = merged$Concordant))

  print(ggplot(tbl, aes(x="", y=Freq, fill=Concordant)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
		coord_polar("y", start=0) +
    scale_fill_brewer(palette = "Set1") +
		theme_minimal() +
		theme(
			axis.text.x=element_blank(),
			axis.title = element_blank(),
			axis.ticks = element_blank(),
			panel.grid  = element_blank()
		) +
    ggtitle('Seurat/MultiSeq Concordance')
	)

  discord <- merged[!merged$GlobalConcordant,]
  discord <- discord %>% group_by(HTO_classification.global.MultiSeq, HTO_classification.global.Seurat) %>% summarise(Count = dplyr::n())

  print(qplot(x=HTO_classification.global.MultiSeq, y=HTO_classification.global.Seurat, data=discord, fill=Count, geom="tile") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
          ggtitle('Discordance By Global Call') + ylab('Seurat') + xlab('MULTI-seq')
  )

  discord <- merged[!merged$Concordant,]
  discord <- discord %>% group_by(HTO_classification.MultiSeq, HTO_classification.Seurat) %>% summarise(Count = dplyr::n())

  print(qplot(x=HTO_classification.MultiSeq, y=HTO_classification.Seurat, data=discord, fill=Count, geom="tile") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
          ggtitle('Discordance By HTO Call') + ylab('Seurat') + xlab('MULTI-seq')
  )

  discord <- merged[!merged$HasSeuratCall | !merged$HasMultiSeqCall,]
  discord <- discord[discord$HasSeuratCall | discord$HasMultiSeqCall,]
  discord <- discord %>% group_by(HTO_classification.MultiSeq, HTO_classification.Seurat) %>% summarise(Count = dplyr::n())
  print(qplot(x=HTO_classification.MultiSeq, y=HTO_classification.Seurat, data=discord, fill=Count, geom="tile") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    ggtitle('Calls Made By Single Caller') + ylab('Seurat') + xlab('MULTI-seq')
  )

  # These calls should be identical, except for possibly negatives from one method that are non-negative in the other
  # For the time being, accept those as correct.
  merged$FinalCall <- as.character(merged$HTO_classification.MultiSeq)
  merged$FinalCall[merged$HTO_classification.MultiSeq == 'Negative'] <- as.character(merged$HTO_classification.Seurat[merged$HTO_classification.MultiSeq == 'Negative'])
  merged$FinalCall[!merged$Concordant] <- 'Discordant'
  merged$FinalCall <- naturalsort::naturalfactor(merged$FinalCall)

  merged$FinalClassification <- as.character(merged$HTO_classification.global.MultiSeq)
  merged$FinalClassification[merged$HTO_classification.global.MultiSeq == 'Negative'] <- as.character(merged$HTO_classification.global.Seurat[merged$HTO_classification.global.MultiSeq == 'Negative'])
  merged$FinalClassification[!merged$GlobalConcordant] <- 'Discordant'
  merged$FinalClassification[!merged$Concordant] <- 'Discordant'
  merged$FinalClassification <- as.factor(merged$FinalClassification)

  df <- data.frame(
    TotalSinglet = c(sum(merged$HTO_classification.global.Seurat == 'Singlet'), sum(merged$HTO_classification.global.MultiSeq == 'Singlet'), sum(merged$FinalClassification == 'Singlet')),
    ConcordantSinglet = c(sum(merged$Concordant & merged$HTO_classification.global.Seurat == 'Singlet'), sum(merged$Concordant & merged$HTO_classification.global.MultiSeq == 'Singlet'), sum(merged$FinalClassification == 'Singlet'))
  )
  rownames(df) <- c('Seurat', 'MultiSeq', 'Final')
  df <- t(df)
  print(kableExtra::kbl(df) %>% kableExtra::kable_styling())

  if (!is.na(allCallsOutFile) && nrow(merged) > 0) {
    write.table(merged, file = allCallsOutFile, row.names = F, sep = '\t', quote = F)
  }

  if (nrow(merged) > 0){
    dt <- data.frame(CellBarcode = merged$Barcode, HTO = merged$FinalCall, HTO_Classification = merged$FinalClassification, key = 'CellBarcode', Seurat = merged$HasSeuratCall, MultiSeq = merged$HasMultiSeqCall)
    dt <- PrintFinalSummary(dt, barcodeMatrix)
    write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)

    return(dt)

  } else {
    print('No rows, not saving ')
  }
}


utils::globalVariables(
  names = c('HTO_Classification', 'TotalCounts'),
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

  htoNames <- simplifyHtoNames(as.character(merged$HTO))

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
    scale_fill_brewer(palette = "Set1") +
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

  print(ggplot(cellData, aes(x = TotalCounts, fill = HTO_Classification)) +
    geom_density() +
    scale_fill_brewer(palette = "Set1") +
    xlab('Total Counts/Cell (log)') +
    ylab('Density') +
    ggtitle('Total Counts By HTO') +
    facet_grid(HTO_Classification ~ ., scales = 'free')
  )

  print(ggplot(cellData[!(cellData$HTO %in% c('Negative', 'Doublet', 'Discordant')),], aes(x = Count, fill = HTO)) +
    geom_density() +
    scale_fill_brewer(palette = "Set1") +
    xlab('HTO Counts/Cell (log)') +
    ylab('Density') +
    ggtitle('Counts By HTO') +
    facet_grid(HTO ~ ., scales = 'free')
  )

  if (sum(merged$HTO == 'Negative') == 0) {
    print('There were no negative cells')
  } else {
    #Melt data:
    melted <- as.data.frame(merged)
    melted <- melted[melted$HTO == 'Negative', !(colnames(melted) %in% c('HTO_Classification', 'HTO', 'key', 'Seurat', 'MultiSeq', 'Count', 'TotalCounts')), drop = FALSE]
    print(str(melted))
    melted <- tidyr::gather(melted, key = 'HTO', value = 'Count', -CellBarcode)

    htoNames <- simplifyHtoNames(as.character(melted$HTO))
    melted$HTO <- naturalfactor(as.character(htoNames))
    melted$Count <- log10(melted$Count + 0.5)

    print(ggplot(melted, aes(x = Count, fill = HTO)) +
      geom_density() +
      scale_fill_brewer(palette = "Set1") +
      xlab('HTO Counts/Cell (log)') +
      ylab('Density') +
      ggtitle('Raw HTO Counts For Negative Cells (log10)') +
      facet_grid(HTO ~ ., scales = 'free')
    )
  }

  tbl <- table(HTO_Classification = merged$HTO_Classification)
  df <- data.frame(tbl)
  df$Pct <- round((df$Freq / sum(df$Freq)) * 100, 2)
  print(kableExtra::kbl(df) %>% kableExtra::kable_styling())

  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(max(3, min(9, length(names(tbl)))), "Set1"))
  colorValues <- getPalette(length(names(tbl)))

  print(ggplot(df, aes(x = '', y=Freq, fill=HTO_Classification)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start=0) +
    scale_fill_brewer(palette = "Set1") +
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



simplifyHtoNames <- function(v) {
  return(sapply(v, function(x){
    x <- unlist(strsplit(x, '-'))
    if (length(x) > 1) {
      x <- x[-(length(x))]
    }

    paste0(x, collapse = "-")
  }))
}



#' @title GenerateSummaryForExpectedBarcodes
#'
#' @description A description
#' @export
GenerateSummaryForExpectedBarcodes <- function(dt, whitelistFile, outputFile, barcodeMatrix) {
  categoryName <- "Cell Hashing Concordance"

  whitelist <- utils::read.table(whitelistFile, sep = '\t', header = F)
  names(whitelist) <- c('CellBarcode')
  df <- data.frame(Category = categoryName, MetricName = "InputBarcodes", Value = length(whitelist$CellBarcode))


  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalCounts", Value = sum(barcodeMatrix)))

  #Any called cell:
  calledCellBarcodes <- dt$CellBarcode[dt$HTO_Classification != 'Negative' & dt$HTO_Classification != 'Discordant']
  calledIntersect <- intersect(whitelist$CellBarcode, calledCellBarcodes)

  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalCalled", Value = length(calledCellBarcodes)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalCalledOverlapping", Value = length(calledIntersect)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionOfInputCalled", Value = (length(calledIntersect) / length(whitelist$CellBarcode))))

  totalCalledNotInInput <- length(calledCellBarcodes) - length(calledIntersect)
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalCalledNotInInput", Value = totalCalledNotInInput))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionCalledNotInInput", Value = (totalCalledNotInInput / length(calledCellBarcodes))))

  #Singlets:
  singletCellBarcodes <- dt$CellBarcode[dt$HTO_Classification == 'Singlet']
  singletIntersect <- intersect(whitelist$CellBarcode, singletCellBarcodes)
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalSinglet", Value = length(singletCellBarcodes)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalSingletOverlapping", Value = length(singletIntersect)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionOfInputSinglet", Value = (length(singletIntersect) / length(whitelist$CellBarcode))))

  totalSingletNotInInput <- length(singletCellBarcodes) - length(singletIntersect)
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalSingletNotInInput", Value = totalSingletNotInInput))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionSingletNotInInput", Value = (totalSingletNotInInput / length(singletCellBarcodes))))

  #Doublets:
  doubletCellBarcodes <- dt$CellBarcode[dt$HTO_Classification == 'Doublet']
  doubletIntersect <- intersect(whitelist$CellBarcode, doubletCellBarcodes)
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalDoublet", Value = length(doubletCellBarcodes)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalDoubletOverlapping", Value = length(doubletIntersect)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionOfInputDoublet", Value = (length(doubletIntersect) / length(whitelist$CellBarcode))))

  totalDoubletNotInInput <- length(doubletCellBarcodes) - length(doubletIntersect)
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalDoubletNotInInput", Value = totalDoubletNotInInput))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionDoubletNotInInput", Value = (totalDoubletNotInInput / length(doubletCellBarcodes))))

  #Discordant:
  discordantCellBarcodes <- dt$CellBarcode[dt$HTO_Classification == 'Discordant']
  discordantIntersect <- intersect(whitelist$CellBarcode, discordantCellBarcodes)
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalDiscordant", Value = length(doubletCellBarcodes)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionOfInputDiscordant", Value = (length(discordantIntersect) / length(whitelist$CellBarcode))))


  #By caller:
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "SeuratNonNegative", Value = sum(dt$HasSeuratCall & dt$HTO_Classification != 'Discordant')))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "MultiSeqNonNegative", Value = sum(dt$HasMultiSeqCall & dt$HTO_Classification != 'Discordant')))

  df$Value[is.na(df$Value)] <- 0

  write.table(df, file = outputFile, quote = F, row.names = F, sep = '\t')
}


#' @title DownloadAndAppendCellHashing
#'
#' @description Downloads matching Cell Hashing data using barcodePrefix on the seurat object and appends it to metadata
#' @param seuratObject, A Seurat object.
#' @return A modified Seurat object.
#' @export
DownloadAndAppendCellHashing <- function(seuratObject, outPath = '.'){
  if (is.null(seuratObject[['BarcodePrefix']])){
    stop('Seurat object lacks BarcodePrefix column')
  }

  for (barcodePrefix in unique(unique(unlist(seuratObject[['BarcodePrefix']])))) {
    print(paste0('Possibly adding cell hashing data for prefix: ', barcodePrefix))

    cellHashingId <- .FindMatchedCellHashing(barcodePrefix)
    if (is.null(cellHashingId)){
      print(paste0('Cell hashing not used for prefix: ', barcodePrefix, ', skipping'))
      next
    } else if (is.na(cellHashingId)){
      stop(paste0('Unable to find cellHashing calls table file for prefix: ', barcodePrefix))
    }

    callsFile <- file.path(outPath, paste0(barcodePrefix, '_cellHashingCalls.csv'))
    DownloadOutputFile(outputFileId = cellHashingId, outFile = callsFile, overwrite = T)
    if (!file.exists(callsFile)){
      stop(paste0('Unable to download calls table for prefix: ', barcodePrefix))
    }

    seuratObject <- AppendCellHashing(seuratObj = seuratObject, barcodeCallFile = callsFile, barcodePrefix = barcodePrefix)
  }

  return(seuratObject)
}

