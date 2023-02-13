context("scRNAseq")

source('testing-data.R')

expectedCalls <- list(
	'MS-11' = 3729,
	'MS-12' = 4494,
	'Negative' = 2055,
	'Doublet' = 811
)

# https://bioconductor.org/packages/release/bioc/vignettes/demuxmix/inst/doc/demuxmix.html
test_that("demuxmix Works", {
	rawData <- '../testdata/438-21-GEX/umi_count'
	h5File <- '../testdata/438-21-GEX/438-21-raw_feature_bc_matrix.h5'
	barcodeMatrix <- ProcessCountMatrix(rawCountData = rawData, barcodeWhitelist = c('MS-11', 'MS-12'))
	
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeMatrix, methods = c('demuxmix'), rawFeatureMatrixH5 = h5File)
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

	cellhashR::CallAndGenerateReport(rawCountData = rawData, barcodeWhitelist = c('MS-11', 'MS-12'), reportFile = html, callFile = output, skipNormalizationQc = T, doTSNE = F, methods = c('demuxmix', 'gmm_demux'), methodsForConsensus = c('demuxmix'), rawFeatureMatrixH5 = h5File)
	
	df <- read.table(output, sep = '\t', header = T)
	print(table(df$consensuscall))
	
	for (hto in unique(df$consensuscall)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}
})