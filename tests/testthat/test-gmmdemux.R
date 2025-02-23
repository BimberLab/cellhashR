context("scRNAseq")

source('testing-data.R')

test_that("GMM-Demux Works", {
	barcodeMatrix <- t(as.matrix(read.csv('../testdata/MS/cell_type_counts.csv', row.names = 1)))
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeMatrix, methods = c('gmm_demux'))
	print(table(df$consensuscall))
	
	expectedCalls <- list(
		Bar1 = 192,
		Bar2 = 140,
		Bar3 = 372,
		Bar4 = 466,
		Bar5 = 147,
		Bar6 = 488,
		Bar7 = 559,
		Bar8 = 524,
		Bar9 = 444,
		Bar10 = 496,
		Bar11 = 476,
		Bar12 = 633,
		Doublet = 605,
		Negative = 254
	)

	for (hto in unique(df$consensuscall)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto, tolerance = 3)
	}
})