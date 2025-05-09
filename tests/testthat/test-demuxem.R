context("scRNAseq")

source('testing-data.R')

expectedCalls <- list(
	'MS-11' = 4388,
	'MS-12' = 6184,
	'Negative' = 517
)

test_that("demuxEM Works", {
	rawData <- '../testdata/438-21-GEX/umi_count'
	h5File <- '../testdata/438-21-GEX/438-21-raw_feature_bc_matrix.h5'
	barcodeMatrix <- ProcessCountMatrix(rawCountData = rawData, barcodeWhitelist = c('MS-11', 'MS-12'))
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeMatrix, methods = c('demuxem'), rawFeatureMatrixH5 = h5File)
	print(table(df$consensuscall))

	for (hto in unique(df$consensuscall)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}
})

test_that("demuxEM Works as a report", {
	html <- paste0(getwd(), '/test.html')
	output <- paste0(getwd(), '/test.txt')
	
	rawData <- '../testdata/438-21-GEX/umi_count'
	h5File <- '../testdata/438-21-GEX/438-21-raw_feature_bc_matrix.h5'

	cellhashR::CallAndGenerateReport(rawCountData = rawData, barcodeWhitelist = c('MS-11', 'MS-12'), reportFile = html, callFile = output, skipNormalizationQc = T, methods = c('demuxem', 'gmm_demux'), methodsForConsensus = c('demuxem'), rawFeatureMatrixH5 = h5File, minAllowableDoubletRateFilter = 0.5)
	
	df <- read.table(output, sep = '\t', header = T)
	print(table(df$consensuscall))
	
	for (hto in unique(df$consensuscall)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}
})