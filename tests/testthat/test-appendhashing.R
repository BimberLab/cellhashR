context("appendhashing")

test_that("AppendCellHashing Works", {
	library(Seurat)
	library(dplyr)
	
	mat <- cellhashR:::.LoadCountMatrix(rawCountData = '../testdata/458-6-GEX.citeSeqCounts.2.citeseq/umi_count')
	mat <- mat[,1:100]

	datasetId1 <- '12345'
	countMat <- matrix(rep(0, length(mat)), nrow = nrow(mat), ncol = ncol(mat))
	rownames(countMat) <- rownames(mat)
	colnames(countMat) <- paste0(datasetId1, '_', colnames(mat))
	so1 <- Seurat::CreateSeuratObject(counts = Seurat::as.sparse(countMat), project = datasetId1)
	so1$BarcodePrefix <- datasetId1
	so1$DatasetId <- datasetId1
	
	datasetId2 <- 'ABCD'
	countMat <- matrix(rep(0, length(mat)), nrow = nrow(mat), ncol = ncol(mat))
	rownames(countMat) <- rownames(mat)
	colnames(countMat) <- paste0(datasetId2, '_', colnames(mat))

	# suppressWarnings added to avoid seurat 'missing row.names' warning
	so2 <- Seurat::CreateSeuratObject(counts = Seurat::as.sparse(countMat), project = datasetId1)
	so2$BarcodePrefix <- datasetId2
	so2$DatasetId <- datasetId2
	
	seuratObj <- merge(so1, so2)
	
	callFile <- 'test-calls.txt'
	
	cells <- colnames(mat)[1:50]
	htos <- sample(c('HTO-1', 'HTO-2'), size = length(cells), replace = TRUE)
	
	callsDf <- data.frame(cellbarcode = cells, consensuscall = htos, consensuscall.global = rep('Singlet', length(cells)))
	rownames(callsDf) <- colnames(so1)[1:50]
	write.table(callsDf, file = callFile, sep = '\t', row.names = FALSE, quote = FALSE)
	seuratObj <- AppendCellHashing(seuratObj, barcodeCallFile = callFile, barcodePrefix = datasetId1)
	
	expect_equal(sum(!is.na(seuratObj$HTO)), ncol(mat))
	expect_equal(sum(!is.na(seuratObj$HTO) & seuratObj$HTO == 'ND'), length(cells))
	
	meta <- seuratObj@meta.data[seuratObj@meta.data$DatasetId == datasetId1,]
	callsDf <- callsDf[rownames(meta),]
	callsDf$consensuscall[is.na(callsDf$cellbarcode)] <- 'ND'
	
	expect_equal(as.character(meta$HTO), as.character(callsDf$consensuscall))
	
	cells <- colnames(mat)[25:75]
	htos <- sample(c('HTO-3', 'HTO-4'), size = length(cells), replace = TRUE)
	callsDf <- data.frame(cellbarcode = cells, consensuscall = htos, consensuscall.global = rep('Doublet', length(cells)), prefixedCells = colnames(so2)[25:75])
	callsDf <- dplyr::arrange(callsDf, desc(htos)) #Skew order
	
	write.table(callsDf, file = callFile, sep = '\t', row.names = FALSE, quote = FALSE)
	seuratObj <- AppendCellHashing(seuratObj, barcodeCallFile = callFile, barcodePrefix = datasetId2)

	meta <- seuratObj@meta.data[seuratObj@meta.data$DatasetId == datasetId2,]
	expect_equal(sum(!is.na(meta$HTO)), ncol(mat))
	expect_equal(sum(!is.na(meta$HTO) & meta$HTO != 'ND'), length(cells))
	
	unlink(callFile)	
})
