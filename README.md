==========================================================
vivlib
==========================================================

*by @vivekbhr*

[![Build Status](https://travis-ci.org/vivekbhr/vivlib.svg?branch=master)](https://travis-ci.org/vivekbhr/vivlib) [![DOI](https://zenodo.org/badge/51520252.svg)](https://zenodo.org/badge/latestdoi/51520252)

## Useful functions and wrappers around RNAseq analysis

*vivlib* is a set of useful functions for RNA-Seq analysis built around featureCounts and edgeR/DESeq2. It includes exploring featureCounts output, performing edgeR/DESeq2, annotating the outputs and some downstream analysis. Some useful operations on BAM files and bed files are also included.

### Installation

You can install vivlib using either `biocLite` from bioconductor or using `install_github` command from devtools package.

Using biocLite :


```r
source("https://bioconductor.org/biocLite.R")
biocLite("vivekbhr/vivlib")
```

### Functions

To get a list of available functions, install vivlib and do :

```r
library(vivlib)
ls(pos = "package:vivlib")
```

### Disclaimer

Although most functions in the latest release have been tested, the package is always under active (but infrequent) development.
