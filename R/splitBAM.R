

#' Split a bam file by given criteria
#'
#' @param bamFile Bam file path
#' @param splitby How to split the file ("strand","chr","flag","mapq")
#' @param chrnameList list of Chromosomes to split by
#' @param flagList List of samflags to split by
#' @param mapqList List of mapq cutoffs to split by (reads with mapq greater than cutoff will be kept)
#' @param outfile_prefix prefix for output filename
#' @param nthreads Number of threads to use
#'
#' @return Splitted, sorted and indexed bam files
#' @export
#'
#' @examples
#' bam <- system.file("extdata", "test_2L2R.bam", package="vivlib")
#' splitBAM(bam = bam, splitby = "chr", chrnameList = NULL, flagList = NULL,
#'              mapqList = NULL, outfile_prefix = "test", nthreads = 20)
#'

splitBAM <- function(bamFile, splitby = "strand", chrnameList = NULL, flagList = NULL,
                     mapqList = NULL, outfile_prefix, nthreads = 20) {

        ## make filterRules (split by chr)
        # splitby could be "chr", "flag", "strand", "mapq"

        ##  ---- Make filterfuncs for every kind of filter ---- ##

        message(paste0("Making rules to filter by : ", splitby))

        # chrname
        make_FilterFunc_chr <- function(name){
                function(df){
                        df_pos = data.frame(data = df$rname , stringsAsFactors = FALSE)
                        return(grepl(name, df_pos$data))
                }
        }
        # flag
        make_FilterFunc_flag <- function(name){
                function(df){
                        df_pos = data.frame(data = df$flag , stringsAsFactors = FALSE)
                        return(grepl(name, df_pos$data))
                }
        }

        # strand
        make_FilterFunc_strand <- function(name){
                function(df){
                        df_pos = data.frame(data = df$strand , stringsAsFactors = FALSE)
                        return(grepl(name, df_pos$data))
                }
        }

        # mapq
        make_FilterFunc_mapq <- function(name){
                function(df){
                        df_pos = data.frame(data = df$mapq , stringsAsFactors = FALSE)
                        return(df_pos$data > name)
                }
        }


        ## -------------  Apply filter func as per the filter type --------- ##

        if(splitby == "chr"){
                if(is.null(chrnameList)){
                        print("Chrname not given for filtering. Splitting file by all chromosomes. Keeping only mapped reads")
                        bam.dat = Rsamtools::scanBam(bamFile, param = Rsamtools::ScanBamParam(what = "rname"),
                                                     flag = scanBamFlag(isUnmappedQuery = FALSE) )
                        chrnameList = unique(as.character(bam.dat[[1]]$rname))
                }
                filtfuncs <- lapply(chrnameList,make_FilterFunc_chr)
                nameSuffix <- chrnameList

        } else if(splitby == "flag") {
                if(is.null(flagList)){
                        stop("Provide flags to split by!")
                } else {
                        filtfuncs <- lapply(flagList,make_FilterFunc_flag)
                        nameSuffix <- flagList
                }

        } else if(splitby == "strand") {
                strandList <- c("+","-")
                filtfuncs <- lapply(strandList,make_FilterFunc_strand)
                nameSuffix <- c("plus","minus")

        } else if(splitby == "mapq") {
                if(is.null(mapqList)){
                        stop("Provide at least one mapq cutoff!")
                } else {
                        filtfuncs <- lapply(mapqList,make_FilterFunc_mapq)
                        nameSuffix <- mapqList
                }

        }

        ##     ---------------------------------------          ##

        ## Make filter rules
        make_FilterRules <- function(FilterFunc){
                return(S4Vectors::FilterRules(list(FilterFunc)) )
        }
        filtrules <- lapply(filtfuncs, make_FilterRules)

        ## Now filter

        message("Filtering the BAM file")
        destinations <- paste0(outfile_prefix,"_", nameSuffix,".bam")

        param = BiocParallel::MulticoreParam(workers = nthreads)
        BiocParallel::bplapply(seq_along(destinations), function(i, file, destinations, filtrules) {
                Rsamtools::filterBam(file, destinations[i], filter = filtrules[[i]])
        }, bamFile, destinations, filtrules,
        BPPARAM = param)

        message("Done!")
}
