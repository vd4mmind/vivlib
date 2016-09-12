# vivlib
Vivek Bhardwaj  

## Useful functions and wrappers around RNAseq analysis

vivlib is a set of useful functions built around RNA-Seq analysis focussed on featurecounts and edgeR/DESeq2. It includes exploring featurecounts output, performing edgeR/DESeq2, annotating the outputs and some downstream analysis. Some useful operations on BAM files and bed files are also included.

### Installation

You can install vivlib using either `biocLite` from bioconductor or using `install_github` command from devtools package.

Using biocLite :


```r
source("https://bioconductor.org/biocLite.R")
biocLite("vivekbhr/vivlib")
```

### Functions

**DISCLAIMER : I keep updating/adding the functions in this package without making any new releases.**

For the updated set of functions, install vivlib and do :


```r
library(vivlib)
ls(pos = "package:vivlib")
```

### Final Disclaimer

Internal testing and version increments for this package is not being done. So use it at your own risk!
