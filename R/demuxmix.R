#' @include Utils.R
#' @include Visualization.R

utils::globalVariables(
  names = c('sortOrder'),
  package = 'cellhashR',
  add = TRUE
)

# rawFeatureMatrixH5 is the raw_feature_bc_matrix.h5 file from the gene expression data
GenerateCellHashCallsDemuxmix <- function(barcodeMatrix, rawFeatureMatrixH5, methodName = 'demuxmix', label = 'demuxmix', verbose= TRUE, metricsFile = NULL, doTSNE = TRUE) {
  if (verbose) {
    print('Starting demuxmix')
  }

  set.seed(GetSeed())

  tryCatch({
    # Calculate nFeature using the h5 file:
    gexCounts <- DropletUtils::read10xCounts(rawFeatureMatrixH5)
    gexData <- SingleCellExperiment::counts(gexCounts)
    colnames(gexData) <- gexCounts$Barcode
    rm(gexCounts)

    if (length(intersect(colnames(colnames(gexData)), colnames(barcodeMatrix))) == 0) {
      print('Removing cell barcode suffixes from input GEX matrix')
      colnames(gexData) <- sapply(colnames(gexData), function(x){
        return(unlist(strsplit(x, split = '-[0-9]'))[1])
      })
    }

    if (length(intersect(colnames(gexData), colnames(barcodeMatrix))) == 0) {
      stop('The cellbarcodes in the HTO and GEX matricies do not overlap. Please check inputs files')
    }

    if (!all(colnames(barcodeMatrix) %in% colnames(gexData))) {
      stop('There are cell barcodes present in the HTO matrix missing from the GEX matrix. Please check inputs files')
    }

    gexData <- as.matrix(gexData[,colnames(barcodeMatrix)])
    rna <- apply(gexData, 2, function(x){
      return(sum(!is.na(x) & x > 0))
    })

    dmm <- demuxmix::demuxmix(barcodeMatrix, rna = rna)

    demuxmix::plotDmmHistogram(dmm)
    demuxmix::plotDmmPosteriorP(dmm)

    classLabels <- demuxmix::dmmClassify(dmm)

    ind <- vapply(dmm@models, is, logical(1), "RegMixModel")
    if (sum(ind) > 0) {
      demuxmix::plotDmmScatter(dmm)
    }

    df <- data.frame(cellbarcode = rownames(classLabels), classification.global = classLabels$Type, classification = classLabels$HTO)

    df <- df[df$cellbarcode %in% colnames(barcodeMatrix),]
    df$classification[df$classification == ''] <- 'Negative'
    df$classification[is.na(df$classification) | df$classification == 'uncertain'] <- 'Negative'
    df$classification[df$classification == 'negative'] <- 'Negative'
    df$classification[df$classification == 'singlet'] <- 'Singlet'
    df$classification[df$classification == 'multiplet'] <- 'Doublet'
    df$classification[grepl(df$classification, pattern = ',')] <- 'Doublet'

    df$classification.global <- df$classification
    df$classification.global[!df$classification.global %in% c('Negative', 'Doublet')] <- 'Singlet'

    # Ensure order matches input:
    toMerge <- data.frame(cellbarcode = colnames(barcodeMatrix), sortOrder = 1:length(colnames(barcodeMatrix)))
    df <- merge(df, toMerge, by = 'cellbarcode', all.x = FALSE, all.y = FALSE)
    df <- dplyr::arrange(df, sortOrder)
    df <- df[c('cellbarcode', 'classification.global', 'classification')]
    df$classification[is.na(df$classification)] <- 'Negative'
    df$classification.global[is.na(df$classification.global)] <- 'Negative'

    ret <- data.frame(cellbarcode = df$cellbarcode, method = methodName, classification = df$classification, classification.global = df$classification.global, stringsAsFactors = FALSE)

    assay <- 'HTO'
    seuratObj <- Seurat::CreateSeuratObject(barcodeMatrix, assay = assay)

    toMerge <- ret$classification
    names(toMerge) <- ret$cellbarcode
    seuratObj$classification.demuxmix <- toMerge[colnames(seuratObj)]
    seuratObj$classification.demuxmix <- naturalsort::naturalfactor(seuratObj$classification.demuxmix)

    toMerge <- ret$classification.global
    names(toMerge) <- ret$cellbarcode
    seuratObj$classification.global.demuxmix <- toMerge[colnames(seuratObj)]
    seuratObj$classification.global.demuxmix <- naturalsort::naturalfactor(seuratObj$classification.global.demuxmix)
    SummarizeHashingCalls(seuratObj, label = label, columnSuffix = 'demuxmix', assay = assay, doTSNE = F, doHeatmap = F)

    return(ret)
  }, error = function(e){
    print('Error generating demuxmix calls, aborting')
    if (!is.null(e)) {
      print(conditionMessage(e))
      traceback()
    }

    return(NULL)
  })
}
