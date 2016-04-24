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
#' plot_fCountSummary(test.summary,"/long/path/to/cut")
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

