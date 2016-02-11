



#' plot the output of featureCounts summary
#'
#' @param summaryFile featureCounts output.summary
#' @param CutFromHeader Enter the long path that's present in featurecounts output to cut
#'
#' @return a plot
#' @export
#'
#' @examples
#' plotFeatureCounts(test.summary,"/long/path/to/cut")
#'


plotFeatureCounts <- function(summaryFile,CutFromHeader){
        f <- read.table(summaryFile,header=T)
        cuttxt <- gsub("/|-",".",CutFromHeader)
        colnames(f) <- gsub(cuttxt,"",colnames(f))
        f <- f[which(rowMeans(f[,2:ncol(f)]) > 0),]
        fp <- reshape2::melt(f)
        fp$million_reads <- fp$value/1000000
        print(ggplot2::ggplot(fp,ggplot2::aes(variable,million_reads,fill=Status)) +
                ggplot2::geom_bar(stat = "identity",position = "dodge") +
                ggplot2::theme_gray(base_size = 14) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,vjust = 0.6))
        )

}

