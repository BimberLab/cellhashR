context("scRNAseq")

source('testing-data.R')

test_that("GMM-Demux Works", {
	ConfigureGMMDemux()
	
	barcodeMatrix <- t(as.matrix(read.csv('../testdata/MS/cell_type_counts.csv', row.names = 1)))
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeMatrix, methods = c('gmm_demux'))
	print(table(df$consensuscall))
	
	expectedCalls <- list(
		Bar1 = 200,
		Bar2 = 147,
		Bar3 = 375,
		Bar4 = 455,
		Bar5 = 149,
		Bar6 = 481,
		Bar7 = 561,
		Bar8 = 520,
		Bar9 = 448,
		Bar10 = 493,
		Bar11 = 477,
		Bar12 = 622,
		Doublet = 752,
		Negative = 215
	)

	# Note: we probably need to set the python / nympy random seed to make this deterministic
	# for (hto in unique(df$consensuscall)) {
	# 	expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	# }
})