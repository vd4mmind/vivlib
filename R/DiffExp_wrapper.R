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
#' fc <- system.file("extdata", "fcount_mouse.tsv", package="vivlib")
#' DESeq_wrapper(fcountOutput = fc,numReplicates = 3, fdr = 0.01,
#' Output = "deseq_output.tsv", pdfReport = "deseq_report.pdf")


DESeq_wrapper <- function(fcountOutput,numReplicates = 4, fdr = 0.01, Output = "deseq_output.tsv",
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
        maxCooks <- apply(SummarizedExperiment::assays(dds)[["cooks"]],1,max)
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

#' A Wrapper for EdgeR over featurecounts output
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
#' fc <- system.file("extdata", "fcount_mouse.tsv", package="vivlib")
#' EdgeR_wrapper(fcountOutput = fc,numReplicates = 3, fdr = 0.01,
#' Output = "edger_output.tsv", pdfReport = "edger_report.pdf")


EdgeR_wrapper <- function(fcountOutput,numReplicates = 4, fdr = 0.01, Output = "edger_output.tsv",
                          col_order = 1, heatmap_topN = 20, pdfReport = "edger_report.pdf"){

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
                condition <- c(rep("Cnt", numReplicates_cont), rep("KD", numReplicates_test) )
        } else {
                colnames(input) <- c(paste0("KD_",1:numReplicates_test),paste0("Control_",1:numReplicates_cont))
                condition <- c(rep("KD", numReplicates_test), rep("Cnt", numReplicates_cont) )
        }

        # edgeR
        samples <- data.frame(row.names = colnames(input),
                              condition = condition)

        design <- model.matrix(~samples$condition)
        y <- edgeR::DGEList(counts = input, group = samples$condition, remove.zeros = TRUE)
        y <- edgeR::calcNormFactors(y)
        y <- edgeR::estimateDisp(y, design)
        print("Using QLF test")
        fit <- edgeR::glmQLFit(y, design)
        qlf.2vs1 <- edgeR::glmQLFTest(fit,coef = 2)
        # result
        ddr.df <- as.data.frame(edgeR::topTags(qlf.2vs1,n = Inf) )
        df.filt <- ddr.df[which(ddr.df$FDR < fdr),]
        df.plot <- data.frame(Status = c("Up","Down"),
                              Genes = c(length(which(df.filt$logFC > 0)),
                                        length(which(df.filt$logFC < 0))
                              )
        )

        ## select top DE genes for heatmap
        logcpm <- edgeR::cpm(y, prior.count=2, log=TRUE)

        decounts <- edgeR::getCounts(y)
        decounts <- merge(decounts,df.filt,by=0)

        # order by fold change (by abs foldch if only few top genes requested)
        select <- order(abs(decounts$logFC), decreasing=TRUE)
        if(heatmap_topN != "all"){
                select <- select[1:heatmap_topN]
        }
        names <- decounts$Row.names[select]
        data <- logcpm[names,]

        ## Write output
        message("writing results")
        pdf(pdfReport)

        edgeR::plotMeanVar(y)
        edgeR::plotBCV(y)
        edgeR::plotMDS.DGEList(y)

        #plot heatmap
        pheatmap::pheatmap(data,cluster_rows = TRUE, clustering_method = "average",
                           show_rownames=TRUE,
                           cluster_cols=FALSE,main = sprintf("Heatmap : Top %d DE genes (by p-value)",heatmap_topN)
        )

        edgeR::plotMD.DGELRT(qlf.2vs1,main = "MA plot")
        print(ggplot2::ggplot(df.plot,ggplot2::aes(Status,Genes,fill=Status)) +
                      ggplot2::geom_bar(stat = "identity", position = "dodge")
        )

        dev.off()
        # Changing the column names so that other functions built for DESeq output can use them
        colnames(ddr.df) <- c("log2FoldChange","logCPM","F","pvalue","padj")
        write.table(ddr.df,file = Output,sep = "\t",quote = FALSE)
        save(y, fit , file = paste0(Output,"_edgeR.Rdata"))
}


#' annotate the output file from Differential Expression wrapper
#'
#' @param DEoutput Tab-seperated output file from \code{\link{EdgeR_wrapper}} or \code{\link{DESeq_wrapper}}
#' @param Output Annotated output file name
#' @param remote Whether use biomart to annotate file
#' @param genome When remote = TRUE, which genome to use (available = "hg38","mm10","dm6")
#' @param map_file If remote = FALSE, provide a map file (with ENS id and Gene id in column 1 and 2) respectively
#'
#' @return annotated output file
#' @export
#'
#' @examples
#' ed <- system.file("extdata", "edger_output.tsv", package="vivlib")
#' annotate_DEoutput(DEoutput = ed, Output = "edger_output_annotated.tsv", remote = TRUE, genome = "hg38")
#'

annotate_DEoutput <- function(DEoutput, Output, remote = TRUE, genome = "hg38", map_file ){

        seqout <- read.delim(DEoutput,header = TRUE)
        if(remote == TRUE) {
                ext.data <- fetch_annotation(genome)
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
