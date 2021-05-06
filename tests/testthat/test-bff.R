context("scRNAseq")

source('testing-data.R')

test_that("BF Works", {
	barcodeData <- t(as.matrix(read.csv('../testdata/MS/cell_type_counts.csv', row.names = 1)))
	
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = c('bff_quantile'))
})