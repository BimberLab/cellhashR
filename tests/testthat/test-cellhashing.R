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

	test <- tests[['438-21']]
	
	#Subset rows to run quicker:
	countData <- Seurat::Read10X(test$input, gene.column=1, strip.suffix = TRUE)
	countData <- countData[,1:2500]
	
	subsetCountDir = normalizePath('./subsetCounts/', mustWork = FALSE)
	DropletUtils::write10xCounts(path = subsetCountDir, countData, overwrite = TRUE)

	fn <- CallAndGenerateReport(rawCountData = subsetCountDir, reportFile = html, callFile = output, citeSeqCountDir = test$citeSeqCountDir, barcodeWhitelist = test$htos, title = 'Test 1', metricsFile = metricsFile)

	df <- read.table(output, sep = '\t', header = TRUE)
	expect_equal(nrow(df), 2500)
	expect_equal(sum(df$consensuscall == 'MS-12'), 1124)

	expect_true(file.exists(metricsFile))
	metrics <- read.table(metricsFile, sep = '\t', header = FALSE)
	expect_equal(nrow(metrics), 23)

	unlink(html)
	unlink(output)
	unlink(metricsFile)

	# Repeat with skip normalization
	fn <- CallAndGenerateReport(rawCountData = subsetCountDir, reportFile = html, callFile = output, citeSeqCountDir = test$citeSeqCountDir, barcodeWhitelist = test$htos, title = 'Test 1', metricsFile = metricsFile, skipNormalizationQc = TRUE)

	df <- read.table(output, sep = '\t', header = TRUE)
	expect_equal(nrow(df), 2500)
	expect_equal(sum(df$consensuscall == 'MS-12'), 1124)

	unlink(html)
	unlink(output)
	unlink(subsetCountDir, recursive = TRUE)
	unlink(metricsFile)
})

test_that("Saturation plot works", {
  saturation <- PlotLibrarySaturation('../testdata/438-21-GEX/')
	expect_equal(0.36, saturation)
})

test_that("Cell hashing works", {
    for (testName in names(tests)) {
				print(paste0('Running test: ', testName))
				test <- tests[[testName]]

				callsFile <- paste0(testName, '-calls.txt')
				summaryFile <- NULL
				if (!is.null(test$gexBarcodeFile)) {
					summaryFile <- paste0(testName, '-summary.txt')
				}

				l <- DoTest(test, callsFile=callsFile, summaryFile=summaryFile)
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

	callsFile <- paste0(testName, '-calls.txt')
	summaryFile <- NULL
	if (!is.null(test$gexBarcodeFile)) {
		summaryFile <- paste0(testName, '-summary.txt')
	}

	l <- DoTest(test, callsFile=callsFile, summaryFile=summaryFile, methods = c('dropletutils', 'bff_quantile', 'multiseq'), skipNormalizationQc = TRUE)
	barcodeData <- l$barcodeData
	df <- l$df
	metricsFile <- l$metricsFile
	unlink(metricsFile)
	
	print(unique(df$bff_quantile))
	expect_equal(expected = test[['BffQuantile']], object = sum(df$bff_quantile != 'Negative' & df$bff_quantile != 'ND'), info = testName)
})
