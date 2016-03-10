

#' A Wrapper fr DESeq2 over featurecounts output
#'
#' @param fcountOutput featurecounts output (with control and test columns)
#' @param numReplicates Number of replicates (could be an integer if the number is same for control and test,
#'                      or a vector with number of replicates for control and for test seperately)
#' @param fdr fdr Cutoff
#' @param Output Output tab seperated file
#' @param pdfReport A report of processing
#' @param col_order How are the columns arranged? (1 = control->test, 2 = test->control)
#'
#' @return Output file and Pdf Report of DESeq2 analysis
#' @export
#'
#' @examples
#' DESeq_wrapper(fcountOutput = "test.out",numReplicates = 4, fdr = 0.01, Output = "deseq_output.tab",
#'              pdfReport = "deseq_report.pdf")


DESeq_wrapper <- function(fcountOutput,numReplicates = 4, fdr = 0.01, Output = "deseq_output.tab",
                          col_order = 1, pdfReport = "deseq_report.pdf"){

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

        colnames(input) <- ifelse(col_order == 1,
                                  c(paste0("Control_",1:numReplicates_cont),paste0("KD_",1:numReplicates_test)),
                                  c(paste0("KD_",1:numReplicates_test),paste0("Control_",1:numReplicates_cont))
                                )
        # DESeq
        samples <- data.frame(row.names = colnames(input),
                              condition = rep(c("Cnt","KD"),each = numReplicates))
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = input,
                                              colData = samples, design = ~condition)
        dds <- DESeq2::DESeq(dds)
        ddr <- DESeq2::results(dds,alpha = fdr)
        ddr.df <- as.data.frame(ddr)
        df.filt <- dplyr::filter(ddr.df,padj < 0.01)
        df.plot <- data.frame(Status = c("Up","Down"),
                              Genes = c(length(which(df.filt$log2FoldChange > 0)),
                                        length(which(df.filt$log2FoldChange < 0))
                              )
        )
        # Write output
        rld <- DESeq2::rlog(dds)
        select <- order(rowMeans(DESeq2::counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]

        message("writing results")
        pdf(pdfReport)
        DESeq2::plotSparsity(dds)
        DESeq2::plotDispEsts(dds)
        print(DESeq2::plotPCA(rld))
        pheatmap::pheatmap(SummarizedExperiment::assay(rld)[select,],
                           cluster_rows=TRUE, show_rownames=TRUE,
                           cluster_cols=FALSE,main = "Heatmap : Top 20 expressed genes")
        DESeq2::plotMA(ddr)
        print(ggplot2::ggplot(df.plot,ggplot2::aes(Status,Genes,fill=Status)) +
                      ggplot2::geom_bar(stat = "identity", position = "dodge")
        )
        dev.off()

        write.table(ddr.df,file = Output,sep = "\t",quote = FALSE)
}
