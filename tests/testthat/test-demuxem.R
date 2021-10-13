context("scRNAseq")

source('testing-data.R')

test_that("demuxEM Works", {
	rawData <- '../testdata/438-21-GEX/umi_count'
	h5File <- '../testdata/438-21-GEX/438-21-raw_feature_bc_matrix.h5'
	barcodeMatrix <- ProcessCountMatrix(rawCountData = rawData, barcodeWhitelist = c('MS-11', 'MS-12'))
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeMatrix, methods = c('demuxem'), demuxem.rawFeatureMatrixH5 = h5File)
	print(table(df$consensuscall))
	
	expectedCalls <- list(
		Bar1 = 196,
		Bar2 = 139,
		Bar3 = 366,
		Bar4 = 449,
		Bar5 = 158,
		Bar6 = 467,
		Bar7 = 543,
		Bar8 = 505,
		Bar9 = 435,
		Bar10 = 487,
		Bar11 = 457,
		Bar12 = 621,
		Doublet = 752,
		Negative = 204
	)

	# Note: we probably need to set the python / nympy random seed to make this deterministic
	# for (hto in unique(df$consensuscall)) {
	# 	expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	# }
})