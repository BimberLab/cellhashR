context("citeseq")

test_that("CITE-seq plots work", {

	LoadCiteSeqCountData(rawCountData = '../testdata/458-6-GEX.citeSeqCounts.2.citeseq/umi_count')

	PlotLibrarySaturationByMarker('../testdata/458-6-GEX.citeSeqCounts.2.citeseq/')

})
