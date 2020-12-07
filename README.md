# cellhashR
An R package designed to demultiplex cell hashing data.



```

# Example 1: parse CITE-seq-Count output, printing QC
barcodeData <- ProcessCountMatrix(rawCountData = 'myCountDir', minCountPerCell = 5)


# Example 2: parse CITE-seq-Count output, providing a barcode whitelist. 
barcodeData <- ProcessCountMatrix(rawCountData = 'myCountDir', minCountPerCell = 5, barcodeWhitelist = c('HTO-1', 'HTO-2', 'HTO-3', 'HTO-4', 'HTO-6'))

```