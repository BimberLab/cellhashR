context("scRNAseq")

source('testing-data.R')

test_that("demuxEM Works", {
	rawData <- '../testdata/438-21-GEX/umi_count'
	h5File <- '../testdata/438-21-GEX/438-21-raw_feature_bc_matrix.h5'
	barcodeMatrix <- ProcessCountMatrix(rawCountData = rawData, barcodeWhitelist = c('MS-11', 'MS-12'))
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeMatrix, methods = c('demuxem'), demuxem.rawFeatureMatrixH5 = h5File)
	print(table(df$consensuscall))

	expectedCalls <- list(
		'MS-11' = 4380,
		'MS-12' = 6195,
		'Negative' = 514
	)

	for (hto in unique(df$consensuscall)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}
})