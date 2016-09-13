## Wrappers for stuff downstream of featurecounts
# plot featurecount summary, boxplot of expression catgories, correlation between replicates. etc.

#' Plot the output of featureCounts summary
#'
#' @param summaryFile featureCounts output.summary file
#' @param CutFromHeader unnecessory text to remove from the header
#' @param outFile Output file name
#'
#' @return A plot of featurecounts summary in a pdf file
#' @export
#'
#' @examples
#' fcsum <- system.file("extdata", "test_fcountsummary.tsv", package="vivlib")
#' plot_fCountSummary(fcsum,"/long/path/")
#'


plot_fCountSummary <- function(summaryFile, CutFromHeader, outFile = NULL){
        f <- read.table(summaryFile,header=T)
        cuttxt <- gsub("/|-",".",CutFromHeader)
        colnames(f) <- gsub(cuttxt,"",colnames(f))
        f <- f[which(rowMeans(f[,2:ncol(f)]) > 0),]
        fp <- reshape2::melt(f)
        fp$million_reads <- fp$value/1000000
        # plot
        p <- ggplot2::ggplot(fp,ggplot2::aes(variable,million_reads,fill=Status)) +
                ggplot2::geom_bar(stat = "identity",position = "dodge") +
                ggplot2::theme_gray(base_size = 14) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,vjust = 0.6))

        if(!is.null(outFile)) {
                pdf(outFile)
                print(p)
                dev.off()
        } else {
                return(p)
        }
}


#' Transform a featurecount output to a data frame of library-normalized mean counts
#'
#' @param fcoutOutput Featurecount output in a data frame
#' @param samplenames The name of samples (exclusing replicate ID etc). Replicates of the samples
#'                      will be grouped together for calculating mean..
#' @param filterByCount A value, if given it will first filter genes by raw mean expression <= given threshold.
#' @param boxplot Logical. Draw a boxplot or normalized counts?
#'
#' @return A data frame of library-normalized mean counts
#' @export
#'
#' @examples
#' fc <- system.file("extdata", "fcount_mouse.out", package="vivlib")
#' samples <- rep(c("cnt","KD"), each = 3)
#' fcount_meantransform(fcountOutput = fc, samplenames = samples, filterByCount = 1000, boxplot = TRUE)
#'

fcount_meantransform <- function(fcountOutput, samplenames, filterByCount = NULL, boxplot = TRUE){
        # make rownames = geneIDs
        rownames(fcountOutput) <- fcountOutput[,1]
        fcountOutput <- fcountOutput[7:ncol(fcountOutput)]

        # Filter genes by mean count over all samples
        if(!is.null(filterByCount)){
                meanexp <- rowMeans(fcountOutput,na.rm = TRUE)
                oldnum <- nrow(fcountOutput)
                fcountOutput <- fcountOutput[meanexp <= filterByCount,]
                filteredNum <- oldnum - nrow(fcountOutput)
                print(paste0("Filtered ", filteredNum, " entries!"))
        }

        # A function to get normalized mean of replicates by sample
        mean_bysample <- function(name, df){
                df2 <- dplyr::select(df, contains(name))
                # get library-norm counts for the df using DESeq2
                coldata <- data.frame(row.names = colnames(df2), sample = rep(name, ncol(df2)))
                dds <- DESeq2::DESeqDataSetFromMatrix(df2, colData = coldata,design = ~1)
                dds <- DESeq2::estimateSizeFactors(dds)
                df2 <- DESeq2::counts(dds,normalized=TRUE)

                # return rowmeans
                rmeans <- rowMeans(df2)
                return(rmeans)
        }

        output <- sapply(samplenames, mean_bysample, fcountOutput)

        ## make samplewise boxplots if asked
        if(boxplot){
                plotdat <- reshape2::melt(output)
                print(ggplot(plotdat, aes(Var2, value, fill = Var2)) +
                              geom_boxplot(outlier.color = "red",notch = TRUE) +
                              scale_y_log10() +
                              labs(x = "Sample", y = "Counts", fill = "Sample", title = "Mean expression by sample")
                )
        }

        return(output)
}


