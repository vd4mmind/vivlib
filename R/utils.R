
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



#' A function to get normalized mean of replicates by sample
#'
#' @param name name of sample (from the given df)
#' @param df A data frame of counts
#'
#' @return data frame
#'

mean_bysample <- function(name, df){
        df2 <- dplyr::select(df, dplyr::contains(name))
        # get library-norm counts for the df using DESeq2
        coldata <- data.frame(row.names = colnames(df2), sample = rep(name, ncol(df2)))
        dds <- DESeq2::DESeqDataSetFromMatrix(df2, colData = coldata,design = ~1)
        dds <- DESeq2::estimateSizeFactors(dds)
        df2 <- DESeq2::counts(dds,normalized=TRUE)

        # return rowmeans
        rmeans <- rowMeans(df2)
        return(rmeans)
}
