
#' annotate DESeq Output file
#'
#' @param DESeqOutput tab seperated DESeq2 output file
#' @param Output Annotated output file
#' @param remote Whether use biomart to annotate file
#' @param genome When remote = TRUE, which genome to use? (available = "hg38","mm10","dm6")
#' @param map_file If remote = FALSE, provide a map file (with ENS id and Gene id in column 1 and 2) respectively
#'
#' @return annotated output file
#' @export
#'
#' @examples
#' annotate_DESeqOutput(DESeqOutput = "test.out", Output = "test_annotated.out", remote = TRUE, genome = "hg38")
#'

annotate_DESeqOutput <- function(DESeqOutput, Output, remote = TRUE, genome = "hg38", map_file ){

        seqout <- read.delim(DESeqOutput,header = TRUE)
        if(remote == TRUE) {
                genomes = data.frame(id = c("hg38","mm10","dm6"),
                                     path = c("hsapiens_gene_ensembl","mmusculus_gene_ensembl",
                                              "dmelanogaster_gene_ensembl"),
                                     stringsAsFactors = FALSE)

                assertthat::assert_that(genome %in% genomes$id)
                message("fetching annotation")
                tofetch <- dplyr::filter(genomes,id == genome)$path
                ensembl = biomaRt::useMart("ensembl",tofetch)
                ext.data <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","description"),
                                mart = ensembl)
                message("merging and writing")
                outfile <- merge(seqout,ext.data,by.x = 0,by.y = 1,all.x = TRUE)
                write.table(outfile,Output,sep="\t", row.names = FALSE,quote=FALSE)

        } else {
                assertthat::assert_that(assertthat::is.readable(map_file))
                message("merging and writing")
                bed <- read.delim(map_file,header = TRUE)
                outfile <- merge(seqout,bed,by.x = 0,by.y = 1,all.x = TRUE)
                write.table(outfile,Output,sep="\t", row.names = FALSE,quote=FALSE)
        }
        message("Done!")
}
