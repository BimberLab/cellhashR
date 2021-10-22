![R Build and Checks](https://github.com/BimberLab/cellhashR/workflows/R%20Build%20and%20Checks/badge.svg)

# cellhashR
An R package designed to demultiplex cell hashing data. [Please see our documentation for more detail](https://bimberlab.github.io/cellhashR/).

## Table of Contents
* [Overview](#overview)
* [Example Usage](#example)
* [Installation](#installation)
* [Known Issues](#issues)
* [Development Guidelines](#developers)


### <a name = "overview">Overview</a>

Cell hashing is a method that allows sample multiplexing or super-loading within single-cell RNA-seq platforms, such as 10x genomics, originally developed at New York Genome Center in collaboration with the Satija lab. [See here for more detail on the technique](https://cite-seq.com/cell-hashing/). The general idea is that cells are labeled with a staining reagent (such as an antibody) tagged with a short nucleotide barcode. Other staining methods have been published, such as the lipid-based Multi-Seq ([https://www.ncbi.nlm.nih.gov/pubmed/31209384](https://www.ncbi.nlm.nih.gov/pubmed/31209384)).  In all methods, the hashtag oligo/barcode is sequenced in parallel with cellular mRNA, creating a separate cell hashing library. After sequencing, the cell barcode and hashing index are parsed using tools like Cite-seq-Count ([https://github.com/Hoohm/CITE-seq-Count](https://github.com/Hoohm/CITE-seq-Count)), creating a count matrix with the total hash tag counts per cell. 

Once the count matrix is created, an algorithm must be used to demultiplex cells and assign them to hash tags (i.e. sample). This is where cellhashR comes in. This package provides several functions:
- Quality control reports for the cell hashing library, covering read counts and normalization. Think [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), except for cell hashing data.
- A single interface to run one or more demutiplexing algorithms, including the novel demultiplexing algorithms BFF_raw and BFF_cluster.  Each algorithm has pros and cons, and will perform better or worse under certain conditions (though in our experience, of the algorithms we have tested, the BFF algorithms work most consistently and under the widest variety of conditions). If you select multiple algorithms (our default workflow), cellhashR will score cells using the consensus call from the set. Various QC summaries are produced during this process as well, if debugging is needed.  In addition to the BFF demultiplexing algorithms, other algorithms that can be run from cellhashR include:
    - [GMM-Demux](https://github.com/CHPGenetics/GMM-Demux)
    - [demuxEM](https://github.com/klarman-cell-observatory/demuxEM) [(see extra requirements below)](#demuxEM)
    - [deMULTIplex](https://github.com/chris-mcginnis-ucsf/MULTI-seq)
    - [HTODemux from Seurat](https://satijalab.org/seurat/v3.1/hashing_vignette.html)
    - [hashedDrops from DropletUtils](https://github.com/MarioniLab/DropletUtils)
- The workflow produces a unified table with the results of each caller and the consensus call. Final QC plots and summaries are created. 

Each step of the workflow can either be run interactively in R (through the terminal or RStudio), or it can be executed as a pipeline that runs all commands and creates the call table and an HTML report. 

[Click here to view an example QC report](https://bimberlab.github.io/cellhashR/articles/V01-QC-example.html)

### <a name="example">Example Usage</a>

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
    
Pre-packaged Docker images with all needed dependencies installed can be found on our [GitHub Packages page](https://github.com/orgs/BimberLab/packages/container/package/cellhashr). We recommend using a specific release, which you can do using tags:    

```
docker pull ghcr.io/bimberlab/cellhashr:latest
```

### <a name="issues">Known Issues</a>

If you receive an error along the lines of:
```
"ERROR; return code from pthread_create() is 22\n"
```
Please manually install preprocessCore with threading disabled:
```
devtools::install_github('bmbolstad/preprocessCore', dependencies = T, upgrade = 'always', configure.args = '--disable-threading')
```


### <a name="demuxEM">demuxEM</a>

Unlike the other algorithms, which just require the HTO count matrix, demuxEM also requires the 10x h5 gene expression counts. This can be supplied as follows. This example runs BFF and demuxEM:
```
  rawData <- '../testdata/438-21-GEX/umi_count'
  h5File <- '../testdata/438-21-GEX/438-21-raw_feature_bc_matrix.h5'
  barcodeMatrix <- ProcessCountMatrix(rawCountData = rawData, barcodeWhitelist = c('MS-11', 'MS-12'))
  df <- GenerateCellHashingCalls(barcodeMatrix = barcodeMatrix, methods = c('bff_cluster', 'demuxem'), demuxem.rawFeatureMatrixH5 = h5File)
```

### <a name="developers">Development Guidelines</a>

* New development should occur on a branch, and go through a Pull Request before merging into the master branch.  [See here for information on the pull request workflow](https://guides.github.com/introduction/flow/).  Ideally PRs would be reviewed by another person.  For the PR, please review the set of changed files carefully to make sure you are only merging the changes you intend.   

* New functions should have [Roxygen2 documentation](https://kbroman.org/pkg_primer/pages/docs.html).

* As part of each PR, you should run 'devtools::document()' to update documentation and include these changes with your commits.

* It is a good idea to run 'R CMD check' locally to make sure your changes will pass.  [See here for more information](http://r-pkgs.had.co.nz/check.html)

* Code should only be merged after the build and tests pass.  The master branch should always be stable.

* New features should ideally have at least a basic test (see [R testthat](http://r-pkgs.had.co.nz/tests.html)).  There is existing test data in ./tests/testdata.  This can be expanded, but please be conscious about file size and try to reuse data across tests if appropriate.


