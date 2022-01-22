context("scRNAseq")

source('testing-data.R')

test_that("BFF Works", {
	barcodeData <- t(as.matrix(read.csv('../testdata/MS/cell_type_counts.csv', row.names = 1)))

	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = c('bff_cluster'))

	expectedCalls <- list(
		Bar1 = 195,
		Bar2 = 154,
		Bar3 = 431,
		Bar4 = 483,
		Bar5 = 150,
		Bar6 = 535,
		Bar7 = 599,
		Bar8 = 569,
		Bar9 = 487,
		Bar10 = 541,
		Bar11 = 512,
		Bar12 = 687,
		Doublet = 385,
		Negative = 68
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

