#' @include Utils.R
#' @include Preprocessing.R

utils::globalVariables(
  names = c('Marker', 'TotalCount'),
  package = 'cellhashR',
  add = TRUE
)

#' @title Load and Plot CiteSeq Count Data
#'
#' @description Load and create QC plots for CITE-seq data
#' @param rawCountData, The input barcode file or umi_count folder
#' @export
LoadCiteSeqCountData <- function(rawCountData = NA) {
  barcodeMatrix <- .LoadCountMatrix(rawCountData = rawCountData)

  featuresToPlot <- rownames(barcodeMatrix)
  seuratObj <- Seurat::CreateSeuratObject(barcodeMatrix, assay = 'ADT')
  setSize <- 2
  steps <- ceiling(length(featuresToPlot) / setSize) - 1

  for (i in 0:steps) {
    start <- (i * setSize) + 1
    end <- min((start + setSize - 1), length(featuresToPlot))
    features <- featuresToPlot[start:end]

    suppressMessages(print(RidgePlot(seuratObj, assay = 'ADT', features = features, ncol = 1)))
  }

  # Also total per ADT
  countsPerAdt <- rowSums(as.matrix(barcodeMatrix))
  countsPerAdt <- data.frame(Marker = names(countsPerAdt), TotalCount = countsPerAdt)

  P1 <- ggplot(countsPerAdt, aes(x = TotalCount)) +
    geom_density() +
    xlab('Total Count/ADT') +
    ylab('Density') +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    labs(title = 'Total Counts Per ADT')

  print(P1)

  P2 <- ggplot(countsPerAdt, aes(x = Marker, y = TotalCount)) +
    geom_bar(stat = 'identity') +
    xlab('Marker') +
    ylab('Total Count') +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    labs(title = 'Total Counts Per ADT')

  print(P2)
}
