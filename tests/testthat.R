library(testthat)
library(cellhashR)

options(testthat.progress.max_fails = 10000)

test_check("cellhashR", reporter = "progress")