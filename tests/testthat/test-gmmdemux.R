context("scRNAseq")

source('testing-data.R')

test_that("GMM-Demux Works", {
	barcodeMatrix <- t(as.matrix(read.csv('../testdata/MS/cell_type_counts.csv', row.names = 1)))
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeMatrix, methods = c('gmm_demux'))
	print(table(df$consensuscall))
	
	expectedCalls <- list(
		Bar1 = 192,
		Bar2 = 143,
		Bar3 = 377,
		Bar4 = 467,
		Bar5 = 144,
		Bar6 = 498,
		Bar7 = 564,
		Bar8 = 528,
		Bar9 = 448,
		Bar10 = 493,
		Bar11 = 482,
		Bar12 = 637,
		Doublet = 566,
		Negative = 257
	)

	for (hto in unique(df$consensuscall)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}
})