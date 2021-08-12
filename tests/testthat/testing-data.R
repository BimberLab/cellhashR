tests <- list(
		'282-1' = list(
        input = '../testdata/cellHashing/282-1-HTO_cellHashingRawCounts.txt',
        htos = paste0('HTO-', c(2:3, 8, 10, 12)),
        gexBarcodeFile = '../testdata/cellHashing/282-1-whitelist.txt',
        CalledCells = 3476,
        Singlet = 2198,
				Doublet = 947,
        MultiSeqCalled = 4010,
        Discordant = 1524,
        SeuratCalled = 3179
    ),
		'283' = list(
        input = '../testdata/cellHashing/283-cellbarcodeToHTO.calls.citeSeqCounts.txt',
        htos = paste0('HTO-', c(2:6)),
        gexBarcodeFile = '../testdata/cellHashing/283-validCellIndexes.csv',
        CalledCells = 3600,
        Singlet = 2247,
				Doublet = 723,
        MultiSeqCalled = 3223,
        Discordant = 1400,
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
      input = '../testdata/438-24-GEX/umi_count/',
      htos = paste0('MS-', c(11, 12)),
      CalledCells = 4624,
      Singlet = 3920,
			Doublet = 251,
      MultiSeqCalled = 4547,
      Discordant = 376,
      SeuratCalled = 3038,
			BffQuantile = 4942
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
      input = '../testdata/457-1-GEX/umi_count/',
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
      CalledCells = 5000,
      Singlet = 3368,
			Doublet = 580,
      MultiSeqCalled = 5000,
      Discordant = 0,
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

DoTest <- function(test, methods = c('multiseq', 'htodemux'), skipNormalizationQc = FALSE) {
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
