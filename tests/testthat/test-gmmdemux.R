context("scRNAseq")

source('testing-data.R')

test_that("GMM-Demux Works", {
	barcodeMatrix <- t(as.matrix(read.csv('../testdata/MS/cell_type_counts.csv', row.names = 1)))
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeMatrix, methods = c('gmm_demux'))
	print(table(df$consensuscall))
	
	expectedCalls <- list(
		Bar1 = 191,
		Bar2 = 143,
		Bar3 = 373,
		Bar4 = 467,
		Bar5 = 148,
		Bar6 = 487,
		Bar7 = 564,
		Bar8 = 524,
		Bar9 = 447,
		Bar10 = 496,
		Bar11 = 478,
		Bar12 = 634,
		Doublet = 588,
		Negative = 255
	)

	for (hto in unique(df$consensuscall)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}
})