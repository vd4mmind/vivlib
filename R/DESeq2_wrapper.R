
#' A Wrapper for DESeq2 over featurecounts output
#'
#' @param fcountOutput featurecounts output (with control and test columns)
#' @param numReplicates Number of replicates (could be an integer if the number is same for control and test,
#'                      or a vector with number of replicates for control and for test seperately)
#' @param fdr fdr Cutoff
#' @param Output Output tab seperated file
#' @param pdfReport A report of processing
#' @param col_order How are the columns arranged? (1 = control->test, 2 = test->control)
#' @param heatmap_topN Number of top DE genes to plot in a heatmap (will make rlog trasformed and clustered heatmap)
#'                      type "all" for all DE genes.
#'
#' @return Output file and Pdf Report of DESeq2 analysis
#' @export
#'
#' @examples
#' DESeq_wrapper(fcountOutput = "test.out",numReplicates = 4, fdr = 0.01, Output = "deseq_output.tab",
#'              pdfReport = "deseq_report.pdf")


DESeq_wrapper <- function(fcountOutput,numReplicates = 4, fdr = 0.01, Output = "deseq_output.tab",
                          col_order = 1, heatmap_topN = 20, pdfReport = "deseq_report.pdf"){

        # Prepare data
        message("reading data")
        data <- read.table(fcountOutput,header = T)
        message("preparing data")
        input <- data.frame(row.names = data[,1], data[,c(7:ncol(data))])
        if (length(numReplicates) != 1) {
                numReplicates_cont <- as.integer(numReplicates[1])
                numReplicates_test <- as.integer(numReplicates[2])
        } else {
                numReplicates_cont <- as.integer(numReplicates)
                numReplicates_test <- as.integer(numReplicates)
        }

        if(col_order == 1) {
                colnames(input) <- c(paste0("Control_",1:numReplicates_cont),paste0("KD_",1:numReplicates_test))
        } else {
                colnames(input) <- c(paste0("KD_",1:numReplicates_test),paste0("Control_",1:numReplicates_cont))
        }

        # DESeq
        samples <- data.frame(row.names = colnames(input),
                              condition = rep(c("Cnt","KD"),each = numReplicates))
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = input,
                                              colData = samples, design = ~condition)
        dds <- DESeq2::DESeq(dds)
        ddr <- DESeq2::results(dds,alpha = fdr)
        ddr.df <- as.data.frame(ddr)
        df.filt <- ddr.df[which(ddr.df$padj < fdr),]
        df.plot <- data.frame(Status = c("Up","Down"),
                              Genes = c(length(which(df.filt$log2FoldChange > 0)),
                                        length(which(df.filt$log2FoldChange < 0))
                              )
        )

        ## select top DE genes for heatmap
        rld <- DESeq2::rlog(dds)

        # order by fold change (by abs foldch if only few top genes requested)
        select <- order(abs(df.filt$log2FoldChange), decreasing=TRUE)
        if(heatmap_topN != "all"){
                select <- select[1:heatmap_topN]
        }
        names <- rownames(df.filt)[select]
        data <- SummarizedExperiment::assay(rld)[names,]

        ## take out cooks distance statistics for outlier detection
        W <- ddr$stat
        maxCooks <- apply(assays(dds)[["cooks"]],1,max)
        idx <- !is.na(W)
        m <- ncol(dds)
        p <- 3

        ## Write output
        message("writing results")
        pdf(pdfReport)

        DESeq2::plotSparsity(dds)
        DESeq2::plotDispEsts(dds)
        print(DESeq2::plotPCA(rld))
        # plot cooks
        plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic",
             ylab="maximum Cook's distance per gene",
             ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
        abline(h=qf(.99, p, m - p))
        #plot heatmap
        pheatmap::pheatmap(data,cluster_rows = TRUE, clustering_method = "average",
                           show_rownames=TRUE,
                           cluster_cols=FALSE,main = sprintf("Heatmap : Top %d DE genes (by p-value)",heatmap_topN)
                        )

        DESeq2::plotMA(ddr)
        print(ggplot2::ggplot(df.plot,ggplot2::aes(Status,Genes,fill=Status)) +
                      ggplot2::geom_bar(stat = "identity", position = "dodge")
        )

        dev.off()

        write.table(ddr.df,file = Output,sep = "\t",quote = FALSE)
        save(dds,ddr, file = paste0(Output,"_DESeq.Rdata"))
}


#' annotate DESeq2 Output file
#'
#' @param DESeqOutput Tab-seperated DESeq2 output file
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
                                     path = c("hsapiens_gene_ensembl",
                                              "mmusculus_gene_ensembl",
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
