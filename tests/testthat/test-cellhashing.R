context("scRNAseq")

source('testing-data.R')

test_that("Cellbarcode Whitelist Works", {
	test <- tests[['438-21']]

	#Subset rows to run quicker:
	countData <- Seurat::Read10X(test$input, gene.column=1, strip.suffix = TRUE)
	countData <- countData[,1:2500]

	cellbarcodeWhitelist <- colnames(countData)[1:200]

	mat <- ProcessCountMatrix(rawCountData = test$input, cellbarcodeWhitelist = cellbarcodeWhitelist)
	expect_equal(colnames(mat), cellbarcodeWhitelist)

	fn <- 'whitelist.txt'
	write.table(cellbarcodeWhitelist, file = fn, row.names = FALSE, col.names = FALSE, quote = FALSE)
	mat <- ProcessCountMatrix(rawCountData = test$input, cellbarcodeWhitelist = fn)
	expect_equal(colnames(mat), cellbarcodeWhitelist)

	unlink(fn)
})

test_that("RMarkdown Copy works", {
	fn <- paste0(getwd(), '/foo.rmd')
	print(paste0('saving file to: ', fn))
	GetExampleMarkdown(dest = fn)
	expect_true(file.exists(fn))
	unlink(fn)
})

test_that("Workflow works", {
	html <- paste0(getwd(), '/test.html')
	output <- paste0(getwd(), '/test.txt')
	metricsFile <- paste0(getwd(), '/metrics.txt')
	rawCountsExport <- paste0(getwd(), '/rawCountsFile.rds')
	md  <- paste0(getwd(), '/test.md')
	
	test <- tests[['438-21']]
	
	#Subset rows to run quicker:
	countData <- Seurat::Read10X(test$input, gene.column=1, strip.suffix = TRUE)
	countData <- countData[,1:2500]
	
	subsetCountDir <- normalizePath('./subsetCounts/', mustWork = FALSE)
	DropletUtils::write10xCounts(path = subsetCountDir, countData, overwrite = TRUE)

	fn <- CallAndGenerateReport(rawCountData = subsetCountDir, reportFile = html, callFile = output, barcodeWhitelist = test$htos, title = 'Test 1', metricsFile = metricsFile, rawCountsExport = rawCountsExport)

	df <- read.table(output, sep = '\t', header = TRUE)
	expect_equal(nrow(df), 2500)
	expect_equal(sum(df$consensuscall == 'MS-12'), 1242)

	expect_true(file.exists(metricsFile))
	metrics <- read.table(metricsFile, sep = '\t', header = FALSE)
	expect_equal(nrow(metrics), 24)

	expect_true(file.exists(rawCountsExport))
	rawCountsMat <- readRDS(file = rawCountsExport)
	expect_equal(nrow(rawCountsMat), 2)
	expect_equal(ncol(rawCountsMat), 2500)
	expect_false(file.exists(md))

	unlink(html)
	unlink(output)
	unlink(metricsFile)
	unlink(rawCountsExport)

	# Repeat with skip normalization
	fn <- CallAndGenerateReport(rawCountData = subsetCountDir, reportFile = html, callFile = output, barcodeWhitelist = test$htos, title = 'Test 1', metricsFile = metricsFile, skipNormalizationQc = TRUE, keepMarkdown = TRUE)

	expect_true(file.exists(md))
	
	df <- read.table(output, sep = '\t', header = TRUE)
	expect_equal(nrow(df), 2500)
	expect_equal(sum(df$consensuscall == 'MS-12'), 1242)

	# Repeat with callerDisagreementThreshold, which should drop htodemux
	unlink(md)
	fn <- CallAndGenerateReport(rawCountData = subsetCountDir, reportFile = html, callFile = output, barcodeWhitelist = test$htos, title = 'Test 1', metricsFile = metricsFile, skipNormalizationQc = TRUE, keepMarkdown = TRUE, callerDisagreementThreshold = 0.05)

	expect_true(file.exists(md))
	
	df <- read.table(output, sep = '\t', header = TRUE)
	expect_equal(nrow(df), 2500)
	expect_equal(sum(df$consensuscall == 'MS-12'), 1242)
	expect_equal(sum(df$consensuscall == 'MS-11'), 980)
	expect_equal(sum(df$consensuscall.global == 'Singlet'), 2222)

	unlink(html)
	unlink(output)
	unlink(subsetCountDir, recursive = TRUE)
	unlink(metricsFile)
	unlink(md)
	unlink(paste0(getwd(), '/test_files'), recursive = TRUE)
})

test_that("Saturation plot works", {
	saturation <- PlotLibrarySaturation('../testdata/438-21-GEX/')
	expect_equal(0.36, saturation)
})

test_that("Cell hashing works", {
    for (testName in names(tests)) {
		print(paste0('Running test: ', testName))
		test <- tests[[testName]]

		l <- DoTest(test)
		barcodeData <- l$barcodeData
		df <- l$df
		metricsFile <- l$metricsFile

		expectedHtos <- sort(test$htos)
		actualHtosMatrix <- sort(unname(cellhashR:::SimplifyHtoNames(rownames(barcodeData))))
		expect_equal(expected = expectedHtos, object = actualHtosMatrix, info = testName)

		expect_true(file.exists(metricsFile), info = testName)
		metrics <- read.table(metricsFile, sep = '\t', header = FALSE, col.names = c('Category', 'MetricName', 'MetricValue'))
		expect_equal(object = nrow(metrics), expected = 18 + length(expectedHtos), info = testName)
		unlink(metricsFile)

		print(paste0('evaluating test: ', testName))
		expect_equal(expected = 0, object = sum(df$consensuscall == 'Not Called'), info = testName)
		expect_equal(expected = test[['CalledCells']], object = sum(df$consensuscall != 'Discordant'), info = testName)
		expect_equal(expected = test[['Singlet']], object = sum(df$consensuscall.global == 'Singlet'), info = testName)
		expect_equal(expected = test[['Doublet']], object = sum(df$consensuscall.global == 'Doublet'), info = testName)
		expect_equal(expected = test[['SeuratCalled']], object = sum(df$htodemux != 'Negative'), info = testName)
		expect_equal(expected = test[['MultiSeqCalled']], object = sum(df$multiseq != 'Negative'), info = testName)

		expect_equal(expected = sum(df$consensuscall.global == 'Discordant'), object = sum(df$consensuscall == 'Discordant'), info = testName)
		expect_equal(expected = test[['Discordant']], object = sum(df$consensuscall == 'Discordant'), info = testName)
		expect_equal(expected = test[['Discordant']], object = sum(df$consensuscall.global == 'Discordant'), info = testName)
    }
})


test_that("BFF calling works", {
	testName <- names(tests)[4]

	print(paste0('Running test: ', testName))
	test <- tests[[testName]]

	l <- DoTest(test, methods = c('dropletutils', 'bff_cluster', 'multiseq'), skipNormalizationQc = TRUE)
	barcodeData <- l$barcodeData
	df <- l$df
	metricsFile <- l$metricsFile
	unlink(metricsFile)
	
	print(unique(df$bff_cluster))
	expect_equal(expected = test[['BffQuantile']], object = sum(df$bff_cluster != 'Negative' & df$bff_cluster != 'ND'), info = testName)
})


test_that("Distinct methods and consensus work", {
	testName <- names(tests)[4]
	print(paste0('Running test: ', testName))
	test <- tests[[testName]]
	barcodeFile <- test$input
	barcodeData <- ProcessCountMatrix(rawCountData = barcodeFile, barcodeWhitelist = test$htos)
	if (ncol(barcodeData) > 5000) {
		print('Subsetting barcodeData to 5000 cells')
		barcodeData <- barcodeData[,1:5000]
	}
	
	df <- GenerateCellHashingCalls(barcodeMatrix = Seurat::as.sparse(barcodeData), methods = c('multiseq', 'gmm_demux'), methodsForConsensus = c('gmm_demux'))
	print(table(df$consensuscall))
	
	expectedCalls <- list(
		'MS-11' = 1486,
		'MS-12' = 2432,
		'Negative' = 501,
		'Doublet' = 581
	)
	
	for (hto in unique(df$consensuscall)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}
})


test_that("Consensus call works", {
	testName <- names(tests)[4]
	print(paste0('Running test: ', testName))
	test <- tests[[testName]]
	barcodeFile <- test$input
	barcodeData <- ProcessCountMatrix(rawCountData = barcodeFile, barcodeWhitelist = test$htos)
	if (ncol(barcodeData) > 5000) {
		print('Subsetting barcodeData to 5000 cells')
		barcodeData <- barcodeData[,1:5000]
	}

	print('Run using majorityConsensusThreshold=0.6, which will will allow 2/3 to win:')
	df <- GenerateCellHashingCalls(barcodeMatrix = Seurat::as.sparse(barcodeData), methods = c('multiseq', 'gmm_demux', 'bff_cluster'), majorityConsensusThreshold = 0.6)
	print(table(df$consensuscall))

	expectedCalls <- list(
		'MS-11' = 1612,
		'MS-12' = 2708,
		'Negative' = 59,
		'Doublet' = 621,
		'Discordant' = 0
	)

	for (hto in names(expectedCalls)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}

	print('Run using majorityConsensusThreshold = 0.75. With three callers, this effectively requires all of them to agree')
	df <- GenerateCellHashingCalls(barcodeMatrix = Seurat::as.sparse(barcodeData), methods = c('multiseq', 'gmm_demux', 'bff_cluster'), majorityConsensusThreshold = 0.75)
	print(table(df$consensuscall))

	expectedCalls <- list(
		'MS-11' = 1579,
		'MS-12' = 2692,
		'Negative' = 59,
		'Doublet' = 209,
		'Discordant' = 461
	)

	for (hto in names(expectedCalls)) {
		expect_equal(sum(df$consensuscall == hto), expectedCalls[[hto]], info = hto)
	}
})

test_that("HTO Reporting works", {
	testName <- names(tests)[4]

	print(paste0('Running test: ', testName))
	test <- tests[[testName]]

	barcodeFile <- test$input
	mf <- 'metrics.txt'
	barcodeData <- ProcessCountMatrix(rawCountData = barcodeFile, barcodeWhitelist = c('MS-11', 'MS-9'), metricsFile = mf)

	dat <- read.table(mf, header = FALSE, sep = '\t')
	expect_equal(dat$V3[dat$V2 == 'HTOsDroppedAboveMinRetained'], 'MS-6;MS-10;MS-12')
	unlink(mf)
})

