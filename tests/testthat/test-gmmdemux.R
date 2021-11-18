context("scRNAseq")

source('testing-data.R')

test_that("GMM-Demux Works", {
	barcodeMatrix <- t(as.matrix(read.csv('../testdata/MS/cell_type_counts.csv', row.names = 1)))
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeMatrix, methods = c('gmm_demux'))
	print(table(df$consensuscall))
	
	expectedCalls <- list(
		Bar1 = 201,
		Bar2 = 144,
		Bar3 = 374,
		Bar4 = 454,
		Bar5 = 148,
		Bar6 = 481,
		Bar7 = 555,
		Bar8 = 520,
		Bar9 = 445,
		Bar10 = 493,
		Bar11 = 474,
		Bar12 = 631,
		Doublet = 663,
		Negative = 213
	)

	for (hto in unique(df$consensuscall)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}
})