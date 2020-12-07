context("scRNAseq")

tests <- list(
		'282-1' = list(
        input = '../testdata/cellHashing/282-1-HTO_cellHashingRawCounts.txt',
        htos = c(2:3, 8, 10, 12),
        gexBarcodeFile = '../testdata/cellHashing/282-1-whitelist.txt',
        CalledCells = 6296,
        Singlet = 4207,
        MultiSeq = 5860,
        Discordant = 1704,
        Seurat = 3623,
        TotalRows = 8000,
        DoRowFilter = T
    ),
		'283' = list(
        input = '../testdata/cellHashing/283-cellbarcodeToHTO.calls.citeSeqCounts.txt', htos = c(2:6),
        gexBarcodeFile = '../testdata/cellHashing/283-validCellIndexes.csv',
        CalledCells = 4759,
        Singlet = 3365,
        MultiSeq = 4786,
        Discordant = 1268,
        Seurat = 3459,
        TotalRows = 6027,
        DoRowFilter = T
    )
    # 'NewFormat' = list(
    #     input = '../testdata/cellHashing/umi_count',
    #     htos = c(1:4),
    #     gexBarcodeFile = NULL,
    #     CalledCells = 5346,
    #     Singlet = 2697,
    #     MultiSeq = 853,
    #     Seurat = 2900,
    #     TotalRows = 100,
    #     DoRowFilter = F
    # )
)

test_that("Cell hashing works", {
    for (testName in names(tests)) {
        DoTest <- function(barcodeFile, callsFile, summaryFile, whitelistFile, doRowFilter = F) {

            barcodeData <- ProcessCiteSeqCount(bFile=barcodeFile, doRowFilter = doRowFilter)

            if (nrow(barcodeData) == 0) {
                stop('No passing HTOs')
            }

            if (ncol(barcodeData) == 0) {
                stop('No passing cells')
            }

            GenerateQcPlots(barcodeData)

            #Subset to keep reasonable
            if (ncol(barcodeData) > 8000) {
                print('Subsetting barcodeData')
                barcodeData <- barcodeData[1:8000]
            }

            sc <- cellhashR:::GenerateCellHashCallsSeurat(barcodeData)

            mc <- cellhashR:::GenerateCellHashCallsMultiSeq(barcodeData)

            dt <- cellhashR:::ProcessEnsemblHtoCalls(mc, sc, barcodeData, outFile = callsFile)

            if (!is.null(whitelistFile) && !is.null(summaryFile)){
              cellhashR:::GenerateSummaryForExpectedBarcodes(dt, whitelistFile = whitelistFile, outputFile = summaryFile, barcodeData = barcodeData)
            }

            return(list(barcodeData = barcodeData, dt = dt))
        }

        print(paste0('Running test: ', testName))
        test <- tests[[testName]]

        callsFile <- paste0(testName, '-calls.txt')
        summaryFile <- NULL
        if (!is.null(test$gexBarcodeFile)) {
            summaryFile <- paste0(testName, '-summary.txt')
        }

        l <- DoTest(test$input, callsFile=callsFile, summaryFile=summaryFile, whitelistFile=test$gexBarcodeFile, doRowFilter = test$DoRowFilter)
        barcodeData <- l$barcodeData

        expectedHtos <- sort(paste0('HTO-', test$htos))
        actualHtosMatrix <- sort(unname(cellhashR:::simplifyHtoNames(rownames(barcodeData))))

        expect_equal(expectedHtos, actualHtosMatrix)

        dt <- l$dt

        expect_equal(test[['CalledCells']], sum(dt$HTO_Classification != 'Discordant'))
        expect_equal(test[['Singlet']], sum(dt$HTO_Classification == 'Singlet'))
        expect_equal(test[['Seurat']], sum(dt$Seurat))
        expect_equal(test[['MultiSeq']], sum(dt$MultiSeq))
        expect_equal(test[['Discordant']], sum(dt$HTO == 'Discordant'))
        expect_equal(test[['Discordant']], sum(dt$HTO_Classification == 'Discordant'))

        d <- read.table(callsFile, header = T, sep = '\t')
        expect_equal(test[['TotalRows']], nrow(d))
        unlink(callsFile)

        if (!is.null(summaryFile)) {
            d <- read.table(summaryFile, header = T, sep = '\t')
            expect_equal(21, nrow(d))
            unlink(summaryFile)
        }
    }
})

