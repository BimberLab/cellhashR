context("scRNAseq")

source('testing-data.R')

test_that("BFF Works", {
	barcodeData <- t(as.matrix(read.csv('../testdata/MS/cell_type_counts.csv', row.names = 1)))

	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = c('bff_quantile'))

	expectedCalls <- list(
		Bar1 = 193,
		Bar2 = 151,
		Bar3 = 429,
		Bar4 = 483,
		Bar5 = 147,
		Bar6 = 526,
		Bar7 = 600,
		Bar8 = 564,
		Bar9 = 482,
		Bar10 = 535,
		Bar11 = 507,
		Bar12 = 683,
		Doublet = 428,
		Negative = 68
	)

	for (hto in unique(df$consensuscall)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
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

