

#' A Wrapper fr DESeq2 over featurecounts output
#'
#' @param fcountOutput featurecounts output (with control and test columns)
#' @param numReplicates Number of replicates
#' @param fdr fdr Cutoff
#' @param Output Output tab seperated file
#' @param pdfReport A report of processing
#'
#' @return Output file and Pdf Report of DESeq2 analysis
#' @export
#'
#' @examples
#' DESeq_wrapper(fcountOutput = "test.out",numReplicates = 4, fdr = 0.01, Output = "deseq_output.tab",
#'              pdfReport = "deseq_report.pdf")


DESeq_wrapper <- function(fcountOutput,numReplicates = 4, fdr = 0.01, Output = "deseq_output.tab",
                          pdfReport = "deseq_report.pdf"){
        message("reading data")
        data <- read.table(fcountOutput,header = T)
        message("preparing data")
        input <- data.frame(row.names = data[,1], data[,c(7:ncol(data))])
        colnames(input) <- c(paste0("Control_",1:numReplicates),paste0("KD_",1:numReplicates))
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
                           cluster_rows=FALSE, show_rownames=FALSE,
                           cluster_cols=FALSE,main = "Heatmap : Top 20 expressed genes")
        DESeq2::plotMA(ddr)
        print(ggplot2::ggplot(df.plot,ggplot2::aes(Status,Genes,fill=Status)) +
                      ggplot2::geom_bar(stat = "identity", position = "dodge")
        )
        dev.off()

        write.table(ddr.df,file = Output,sep = "\t",quote = FALSE)
}
