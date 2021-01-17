context("scRNAseq")

tests <- list(
		'282-1' = list(
        input = '../testdata/cellHashing/282-1-HTO_cellHashingRawCounts.txt',
        htos = paste0('HTO-', c(2:3, 8, 10, 12)),
        gexBarcodeFile = '../testdata/cellHashing/282-1-whitelist.txt',
        CalledCells = 3491,
        Singlet = 2213,
				Doublet = 947,
        MultiSeqCalled = 4010,
        Discordant = 1509,
        SeuratCalled = 3179
    ),
		'283' = list(
        input = '../testdata/cellHashing/283-cellbarcodeToHTO.calls.citeSeqCounts.txt',
        htos = paste0('HTO-', c(2:6)),
        gexBarcodeFile = '../testdata/cellHashing/283-validCellIndexes.csv',
        CalledCells = 3715,
        Singlet = 2362,
				Doublet = 723,
        MultiSeqCalled = 3223,
        Discordant = 1285,
        SeuratCalled = 4116
    ),
    '438-21' = list(
      input = '../testdata/438-21-GEX/umi_count',
      citeSeqCountDir = '../testdata/438-21-GEX/',
      htos = paste0('MS-', c(11, 12)),
      CalledCells = 4833,
      Singlet = 3971,
			Doublet = 366,
      MultiSeqCalled = 4504,
      Discordant = 167,
      SeuratCalled = 4096
    ),
    '438-24' = list(
      input = '../testdata/438-24-GEX/umi_count',
      htos = paste0('MS-', c(11, 12)),
      CalledCells = 4624,
      Singlet = 3920,
			Doublet = 251,
      MultiSeqCalled = 4547,
      Discordant = 376,
      SeuratCalled = 3038,
			SeqNDCalled = 3959
    ),
    '449-1' = list(
      input = '../testdata/449-1-GEX/umi_count',
      htos = paste0('MS-', c(2:16)),
      CalledCells = 400,
      Singlet = 295,
			Doublet = 102,
      MultiSeqCalled = 1081,
      Discordant = 706,
      SeuratCalled = 1103
    ),
    '457-1' = list(
      input = '../testdata/457-1-GEX/umi_count',
      htos = paste0('MS-', c(1:3, 5:8)),
      CalledCells = 1816,
      Singlet = 1481,
			Doublet = 330,
      MultiSeqCalled = 2383,
      Discordant = 621,
      SeuratCalled = 2432
    ),
    '471-1' = list(
      input = '../testdata/471-1-GEX/umi_count',
      htos = paste0('MS-', c(1, 2)),
      CalledCells = 3948,
      Singlet = 3368,
			Doublet = 580,
      MultiSeqCalled = 5000,
      Discordant = 1052,
      SeuratCalled = 3948
    ),
    '471-2' = list(
      input = '../testdata/471-2-GEX/umi_count',
      htos = paste0('MS-', c(3, 4)),
      CalledCells = 4966,
      Singlet = 2085,
			Doublet = 553,
      MultiSeqCalled = 38,
      Discordant = 34,
      SeuratCalled = 2672
    ),
    '483-3' = list(
      input = '../testdata/483-3-GEX/umi_count',
      htos = paste0('MS-', c(2:4, 6:8, 10:13)),
      CalledCells = 56,
      Singlet = 45,
			Doublet = 10,
      MultiSeqCalled = 156,
      Discordant = 110,
      SeuratCalled = 165
    )
)

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
	expect_equal(nrow(metrics), 16)

	unlink(html)
	unlink(output)
	unlink(subsetCountDir, recursive = TRUE)
	unlink(metricsFile)
})

test_that("Saturation plot works", {
  saturation <- PlotLibrarySaturation('../testdata/438-21-GEX/')
	expect_equal(0.36, saturation)
})

DoTest <- function(test, callsFile, summaryFile, methods = c('multiseq', 'htodemux'), skipNormalizationQc = FALSE) {
	barcodeFile <- test$input
	barcodeData <- ProcessCountMatrix(rawCountData = barcodeFile, barcodeWhitelist = test$htos)

	if (nrow(barcodeData) == 0) {
		stop('No passing HTOs')
	}

	if (ncol(barcodeData) == 0) {
		stop('No passing cells')
	}

	#Subset to keep reasonable
	if (ncol(barcodeData) > 5000) {
		print('Subsetting barcodeData to 5000 cells')
		barcodeData <- barcodeData[,1:5000]
	}

	# This is giving memory errors on github runners:
	if (!skipNormalizationQc) {
		PlotNormalizationQC(barcodeData)
	}
	
	metricsFile <- 'metrics.txt'
	if (file.exists(metricsFile)) {
		unlink(metricsFile)
	}
	df <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = methods, metricsFile = metricsFile)

	return(list(barcodeData = barcodeData, df = df, metricsFile = metricsFile))
}

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
				expect_equal(expected = expectedHtos, object = actualHtosMatrix)

				expect_true(file.exists(metricsFile))
				metrics <- read.table(metricsFile, sep = '\t', header = FALSE, col.names = c('Category', 'MetricName', 'MetricValue'))
				expect_equal(object = nrow(metrics), expected = 13)
				unlink(metricsFile)

				print(paste0('evaluating test: ', testName))
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


test_that("SeqND calling works", {
	testName <- names(tests)[4]

	print(paste0('Running test: ', testName))
	test <- tests[[testName]]

	callsFile <- paste0(testName, '-calls.txt')
	summaryFile <- NULL
	if (!is.null(test$gexBarcodeFile)) {
		summaryFile <- paste0(testName, '-summary.txt')
	}

	l <- DoTest(test, callsFile=callsFile, summaryFile=summaryFile, methods = c('seqnd', 'multiseq', 'htodemux'), skipNormalizationQc = TRUE)
	barcodeData <- l$barcodeData
	df <- l$df
	metricsFile <- l$metricsFile
	unlink(metricsFile)
	
	expect_equal(expected = test[['SeqNDCalled']], object = sum(df$seqnd != 'Negative'), info = testName)
})
