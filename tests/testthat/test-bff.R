context("scRNAseq")

source('testing-data.R')

test_that("BFF Works", {
	barcodeData <- t(as.matrix(read.csv('../testdata/MS/cell_type_counts.csv', row.names = 1)))

	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = c('bff_cluster'))

	expectedCalls <- list(
		Bar1 = 193,
		Bar2 = 151,
		Bar3 = 429,
		Bar4 = 480,
		Bar5 = 146,
		Bar6 = 526,
		Bar7 = 598,
		Bar8 = 562,
		Bar9 = 482,
		Bar10 = 534,
		Bar11 = 506,
		Bar12 = 683,
		Doublet = 437,
		Negative = 70
	)

	for (hto in unique(df$consensuscall)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}
})

test_that("563-3 Works", {
	barcodeData <- readRDS('../testdata/563-3-TCR.hashing.rawCounts.rds')
	
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = c('bff_cluster'))
})

test_that("568-1 Works", {
	barcodeData <- readRDS('../testdata/568-1-TCR.hashing.rawCounts.rds')

	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = c('bff_cluster'))
})

