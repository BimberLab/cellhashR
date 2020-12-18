context("scRNAseq")

tests <- list(
		'282-1' = list(
        input = '../testdata/cellHashing/282-1-HTO_cellHashingRawCounts.txt',
        htos = paste0('HTO-', c(2:3, 8, 10, 12)),
        gexBarcodeFile = '../testdata/cellHashing/282-1-whitelist.txt',
        CalledCells = 3491,
        Singlet = 2943,
				Doublet = 947,
        MultiSeqCalled = 4010,
        Discordant = 779,
        SeuratCalled = 3179,
        TotalRows = 8000,
        DoRowFilter = T
    ),
		'283' = list(
        input = '../testdata/cellHashing/283-cellbarcodeToHTO.calls.citeSeqCounts.txt',
        htos = paste0('HTO-', c(2:6)),
        gexBarcodeFile = '../testdata/cellHashing/283-validCellIndexes.csv',
        CalledCells = 3715,
        Singlet = 2468,
				Doublet = 723,
        MultiSeqCalled = 3223,
        Discordant = 1285,
        SeuratCalled = 4116,
        TotalRows = 6027,
        DoRowFilter = T
    ),
    '438-21' = list(
      input = '../testdata/438-21-GEX/umi_count',
      citeSeqCountDir = '../testdata/438-21-GEX/',
      htos = paste0('MS-', c(11, 12)),
      CalledCells = 6296,
      Singlet = 4207,
			Doublet = 3,
      MultiSeqCalled = 5860,
      Discordant = 1704,
      SeuratCalled = 3623,
      TotalRows = 8000,
      DoRowFilter = T
    ),
    '438-24' = list(
      input = '../testdata/438-24-GEX/umi_count',
      htos = paste0('MS-', c(11, 12)),
      CalledCells = 6296,
      Singlet = 4207,
			Doublet = 4,
      MultiSeqCalled = 5860,
      Discordant = 1704,
      SeuratCalled = 3623,
      TotalRows = 8000,
      DoRowFilter = T
    ),
    '449-1' = list(
      input = '../testdata/449-1-GEX/umi_count',
      htos = paste0('MS-', c(2:16)),
      CalledCells = 6296,
      Singlet = 4207,
			Doublet = 5,
      MultiSeqCalled = 5860,
      Discordant = 1704,
      SeuratCalled = 3623,
      TotalRows = 8000,
      DoRowFilter = T
    ),
    '457-1' = list(
      input = '../testdata/457-1-GEX/umi_count',
      htos = paste0('MS-', c(1:3, 5:8)),
      CalledCells = 6296,
      Singlet = 4207,
			Doublet = 6,
      MultiSeqCalled = 5860,
      Discordant = 1704,
      SeuratCalled = 3623,
      TotalRows = 8000,
      DoRowFilter = T
    ),
    '471-1' = list(
      input = '../testdata/471-1-GEX/umi_count',
      htos = paste0('MS-', c(1, 2)),
      CalledCells = 6296,
      Singlet = 4207,
			Doublet = 7,
      MultiSeqCalled = 5860,
      Discordant = 1704,
      SeuratCalled = 3623,
      TotalRows = 8000,
      DoRowFilter = T
    ),
    '471-2' = list(
      input = '../testdata/471-2-GEX/umi_count',
      htos = paste0('MS-', c(3, 4)),
      CalledCells = 6296,
      Singlet = 4207,
			Doublet = 8,
      MultiSeqCalled = 5860,
      Discordant = 1704,
      SeuratCalled = 3623,
      TotalRows = 8000,
      DoRowFilter = T
    ),
    '483-3' = list(
      input = '../testdata/483-3-GEX/umi_count',
      htos = paste0('MS-', c(2:4, 6:8, 10:13)),
      CalledCells = 6296,
      Singlet = 4207,
			Doublet = 9,
      MultiSeqCalled = 5860,
      Discordant = 1704,
      SeuratCalled = 3623,
      TotalRows = 8000,
      DoRowFilter = T
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

test_that("Cell hashing works", {
    for (testName in names(tests)) {
        DoTest <- function(test, callsFile, summaryFile) {
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
          PlotNormalizationQC(barcodeData)

					metricsFile <- 'metrics.txt'
					if (file.exists(metricsFile)) {
						unlink(metricsFile)
					}
          df <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = c('multiseq', 'htodemux'), metricsFile = metricsFile)

					return(list(barcodeData = barcodeData, df = df, metricsFile = metricsFile))
				}

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
				metrics <- read.table(metricsFile, sep = '\t', header = FALSE, col.names = c('MetricName', 'MetricValue'))
				expect_equal(object = nrow(metrics), expected = 13)
				unlink(metricsFile)

				print(paste0('evaluating test: ', testName))
				expect_equal(expected = test[['CalledCells']], object = sum(df$consensuscall != 'Discordant'))
				expect_equal(expected = test[['Singlet']], object = sum(df$consensuscall.global == 'Singlet'))
				expect_equal(expected = test[['Doublet']], object = sum(df$consensuscall.global == 'Doublet'))
				expect_equal(expected = test[['SeuratCalled']], object = sum(df$htodemux != 'Negative'))
				expect_equal(expected = test[['MultiSeqCalled']], object = sum(df$multiseq != 'Negative'))
				expect_equal(expected = test[['Discordant']], object = sum(df$consensuscall == 'Discordant'))
				expect_equal(expected = test[['Discordant']], object = sum(df$consensuscall.global == 'Discordant'))
    }
})

