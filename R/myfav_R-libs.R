

## Load a bunch of libraries that I mostly use (NOT WORKING STILL)

## default pkg options
pkglist <- list(
        usual = c("plyr","dplyr","reshape2","ggplot2"),
        genomic = c("GenomicRanges","DESeq2","DEXSeq","devtools"),
        plotlibs = c("UpSetR","RColorBrewer","pheatmaps","gridExtra"),
        efficient = c("parallel","doParallel","readr","data.table")
)

## load pkg list
loadMyLibs <- function(groupName = "usual",silent=TRUE,pkglist = pkglist){
        if(groupName %in% names(pkglist)){
                pkgs = pkglist[[groupName]]
                ## Load packages
                lapply(pkgs,require,character.only=TRUE)
        } else {
                print("Package list is not predefined")
        }

}


