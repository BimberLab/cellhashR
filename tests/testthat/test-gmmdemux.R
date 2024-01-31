context("scRNAseq")

source('testing-data.R')

test_that("GMM-Demux Works", {
	barcodeMatrix <- t(as.matrix(read.csv('../testdata/MS/cell_type_counts.csv', row.names = 1)))
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeMatrix, methods = c('gmm_demux'))
	print(table(df$consensuscall))
	
	expectedCalls <- list(
		Bar1 = 191,
		Bar2 = 146,
		Bar3 = 378,
		Bar4 = 468,
		Bar5 = 145,
		Bar6 = 498,
		Bar7 = 570,
		Bar8 = 529,
		Bar9 = 451,
		Bar10 = 498,
		Bar11 = 485,
		Bar12 = 638,
		Doublet = 546,
		Negative = 258
	)

	for (hto in unique(df$consensuscall)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}
})