
#' Make an informative volcano plot using DESeq2 output
#'
#' @param DESeqOutput Tab-seperated DESeq2 output file
#' @param fdr FDR cutoff for plotting
#' @param foldChangeLine Where to place a line for a fold change cutoff. Needs two values
#'                      (for positive andnegative fold change). Default is to plot no line
#' @param markGenes geneIDs for the genes to mark (a circle with gene name will be added)
#' @param colorGenes geneIDs for the genes to color
#' @param useGeneNames Use gene names instead of default geneIDs. only works if the output
#'                      has been annotated by \code{\link{annotate_DESeqOutput}} function.
#' @param outFile Name of output pdf file
#'
#' @return A pdf file with volcano plot
#' @export
#'
#' @examples
#' plotVolcano("DESeqOutput.tab", fdr = 0.05, foldChangeLine = NULL, markGenes = NULL,
#'              colorGenes = NULL, useGeneNames = TRUE, outFile = "volcano.pdf")
#'

plotVolcano <- function(DESeqOutput, fdr = 0.05, foldChangeLine = NULL, markGenes = NULL,
                        colorGenes = NULL, useGeneNames = TRUE, outFile) {

        # get data and format
        res_obj <- read.delim(DESeqOutput,header = TRUE, stringsAsFactors = TRUE)

        if(isTRUE(useGeneNames)){
                geneID <- res_obj$external_gene_name
        } else {
                geneID <- res_obj$Row.names
        }

        plotdata <- data.frame(geneID = geneID, log2FoldChange = res_obj$log2FoldChange,
                               padj = res_obj$padj )
        # color genes
        if (!(is.null(colorGenes))) {
                plotdata$marked <- ifelse(plotdata$geneID %in% colorGenes, TRUE, FALSE)
        }
        plotdata = plotdata[!is.na(plotdata$padj),]

        # set colors
        xlim = c(-3,3)
        ylim = c(0,30)
        cex = c(0.3,0.5)

        plotdata$cex = cex[1]
        plotdata$pch = 16
        plotdata$col = "#525252"
        plotdata$col[plotdata$padj < fdr] = "#cd0000"

        if (!(is.null(colorGenes))) {
                plotdata$col[plotdata$marked] = "#006400"
        }

        # set limits
        plotdata$pch[plotdata$log2FoldChange < xlim[1] ] = 5
        plotdata$cex[plotdata$log2FoldChange < xlim[1] ] = cex[2]
        plotdata$log2FoldChange[plotdata$log2FoldChange < xlim[1] ] = xlim[1]

        plotdata$pch[plotdata$log2FoldChange > xlim[2] ] = 5
        plotdata$cex[plotdata$log2FoldChange > xlim[2] ] = cex[2]
        plotdata$log2FoldChange[plotdata$log2FoldChange > xlim[2] ] = xlim[2]

        plotdata$pch[-log10(plotdata$padj) > ylim[2] ] = 2
        plotdata$cex[-log10(plotdata$padj) > ylim[2] ] = cex[2]
        plotdata$padj[-log10(plotdata$padj) > ylim[2] ] = 10^ - ylim[2]

        de_up <- length(which(res_obj$log2FoldChange > 0 & res_obj$padj < fdr ))
        de_down <- length(which(res_obj$log2FoldChange < 0 & res_obj$padj < fdr ))

        # plot and save
        pdf(outFile)
        plot(plotdata$log2FoldChange, -log10(plotdata$padj),
             main=sprintf("Volcano plot\n(FDR: %.2f, up: %d, down: %d)",fdr,de_up,de_down),
             xlab="log2-fold change",
             ylab="-log10 adjusted-pvalue",
             xlim=xlim,
             ylim=ylim,
             cex=plotdata$cex, pch=plotdata$pch,
             col=plotdata$col)

        abline(h=-log10(fdr), col=rgb(0,0,1,0.5), lwd=4)
        #abline(v=0, col=rgb(0,0,1,0.5), lwd=4)
        abline(v=0, col= "grey40", lwd=4)

        if(!(is.null(foldChangeLine))) {
                abline(v = log2(foldChangeLine[1]), col=rgb(1,0,1,0.5), lwd=3)
                abline(v = -log2(abs(foldChangeLine[2])), col=rgb(1,0,1,0.5), lwd=3)
        }

        # finally, mark genes
        if (!(is.null(markGenes))) {
                markedata <- plotdata[which(plotdata$geneID %in% markGenes),]
        }
        points(markedata$log2FoldChange,-log10(markedata$padj), pch = 1, lwd = 3, col = "black")
        calibrate::textxy(markedata$log2FoldChange, -log10(markedata$padj),markedata$geneID, cex = 0.6)

        dev.off()
}
