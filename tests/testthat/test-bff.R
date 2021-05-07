context("scRNAseq")

source('testing-data.R')

test_that("BFF Works", {
	barcodeData <- t(as.matrix(read.csv('../testdata/MS/cell_type_counts.csv', row.names = 1)))

	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = c('bff_quantile'))

	expectedCalls <- list(
		Bar1 = 197,
		Bar2 = 155,
		Bar3 = 440,
		Bar4 = 488,
		Bar5 = 151,
		Bar6 = 533,
		Bar7 = 602,
		Bar8 = 576,
		Bar9 = 489,
		Bar10 = 545,
		Bar11 = 520,
		Bar12 = 693,
		Doublet = 315,
		Negative = 92
	)

	for (hto in unique(df$consensuscall)) {
		#TODO: review why non-deterministic?
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto, tolerance = 5)
	}
})

test_that("563-3 Works", {
	barcodeData <- readRDS('../testdata/563-3-TCR.hashing.rawCounts.rds')
	
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = c('bff_quantile'))
})

test_that("568-1 Works", {
	barcodeData <- readRDS('../testdata/568-1-TCR.hashing.rawCounts.rds')

	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = c('bff_quantile'))
})

