
#' Fetch annotation from ENSEMBL ID to gene ID from Biomart
#'
#' @param genome genome name ("hg38","mm10","dm6")
#'
#' @return annotation
#'

fetch_annotation <- function(genome) {
        genomes = data.frame(id = c("hg38","mm10","dm6"),
                             path = c("hsapiens_gene_ensembl",
                                      "mmusculus_gene_ensembl",
                                      "dmelanogaster_gene_ensembl"),
                             stringsAsFactors = FALSE)

        assertthat::assert_that(genome %in% genomes$id)
        message("fetching annotation")
        tofetch <- dplyr::filter(genomes,id == genome)$path
        ensembl <- biomaRt::useMart("ensembl",tofetch)
        ext.data <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","description"),
                                   mart = ensembl)
        return(ext.data)
}
