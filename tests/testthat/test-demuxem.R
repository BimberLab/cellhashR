context("scRNAseq")

source('testing-data.R')

test_that("demuxEM Works", {
	rawData <- '../testdata/438-21-GEX/umi_count'
	h5File <- '../testdata/438-21-GEX/438-21-raw_feature_bc_matrix.h5'
	barcodeMatrix <- ProcessCountMatrix(rawCountData = rawData, barcodeWhitelist = c('MS-11', 'MS-12'))
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeMatrix, methods = c('demuxem'), demuxem.rawFeatureMatrixH5 = h5File)
	print(table(df$consensuscall))
	print(unique(df$consensuscall))
	
	expectedCalls <- list(
		'MS-11' = 4380,
		'MS-11' = 6195,
		'Negative' = 514
	)

	if (sum(is.na(df$consensuscall)) > 0) {
		print('NA found in consensuscall')
		print(df[is.na(df$consensuscall)])
	}

	if (sum(is.null(df$consensuscall)) > 0) {
		print('NULL found in consensuscall')
		print(df[is.null(df$consensuscall)])
	}

	# Note: we probably need to set the python / nympy random seed to make this deterministic
	for (hto in unique(df$consensuscall)) {
		print(hto)
		print(sum(df$consensuscall == hto))
		print(expectedCalls[[hto]])
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}
})