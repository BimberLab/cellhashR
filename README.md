![R Build and Checks](https://github.com/BimberLab/cellhashR/workflows/R%20Build%20and%20Checks/badge.svg)

# cellhashR
An R package designed to demultiplex cell hashing data.


## Table of Contents
* [Example Usage](#example)
* [Installation](#installation)
* [Development Guidelines](#developers)


### <a name="example">Example Usage</a>

[Click here to download an example R Markdown report](https://htmlpreview.github.io/?https://github.com/BimberLab/cellhashR/blob/master/examples/cellhashR.html)

Below are the primary functions of cellhashR needed to QC and score hashing data:
```r
# Example 1: parse CITE-seq-Count output, printing QC
barcodeData <- ProcessCountMatrix(rawCountData = 'myCountDir/umi_count', minCountPerCell = 5)

# Example 2: parse CITE-seq-Count output, providing a barcode whitelist. 
barcodeData <- ProcessCountMatrix(rawCountData = 'myCountDir/umi_count', minCountPerCell = 5, barcodeWhitelist = c('HTO-1', 'HTO-2', 'HTO-3', 'HTO-4', 'HTO-6'))

# Create QC plots of barcode normalization
PlotNormalizationQC(barcodeData)

# Generate the final cell hashing calls
calls <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = c('multiseq', 'htodemux'))

# Inspect negative cells:
SummarizeCellsByClassification(calls = calls, barcodeMatrix = barcodeData)


```

Or export/save a template RMarkdown file outlining the default workflow, which can be run interactively or headlessly as part of a pipeline:
 
```r
GetExampleMarkdown(dest = 'cellhashR_template.rmd')
```

Finally, the workflow can be executed using this wrapper around the Rmarkdown, producing a TSV of calls and HTML QC report:
 
```r
CallAndGenerateReport(rawCountData = 'myCountDir/umi_count', reportFile = 'report.html', callFile = 'calls.txt', barcodeWhitelist = c('HTO-1', 'HTO-2', 'HTO-3'), title = 'Cell Hashing For Experiment 1')
```
### <a name="installation">Installation</a>

```{r}

# Make sure to update your Rprofile to include Bioconductor repos, such as adding this line to ~/.Rprofile:
local({options(repos = BiocManager::repositories())})

#Latest version:
devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'master', dependencies = TRUE, upgrade = 'always')

```
    
Pre-packaged Docker images with all needed dependencies installed can be found on our [github repository](https://hub.docker.com/r/bimberlab/oosap): 

```

docker pull ghcr.io/bimberlab/cellhashR:latest

```

### <a name="developers">Development Guidelines</a>

* New development should occur on a branch, and go through a Pull Request before merging into Dev.  [See here for information on the pull request workflow](https://guides.github.com/introduction/flow/).  Ideally PRs would be reviewed by another person.  For the PR, please review the set of changed files carefully to make sure you are only merging the changes you intend.   

* New functions should have [Roxygen2 documentation](https://kbroman.org/pkg_primer/pages/docs.html).

* As part of each PR, you should run 'devtools::document()' to update documentation and include these changes with your commits.

* It is a good idea to run 'R CMD check' locally to make sure your changes will pass.  [See here for more information](http://r-pkgs.had.co.nz/check.html)

* Code should only be merged after the build and tests pass.  The Dev branch should always be stable.

* New features should ideally have at least a basic test (see [R testthat](http://r-pkgs.had.co.nz/tests.html)).  There is existing test data in ./tests/testdata.  This can be expanded, but please be conscious about file size and try to reuse data across tests if appropriate.


