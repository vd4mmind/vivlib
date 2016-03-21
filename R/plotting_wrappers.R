
#' Make an informative volcano plot using DESeq2 output
#'
#' @param DESeqOutput Tab-seperated DESeq2 output file
#' @param fdr FDR cutoff for plotting
#' @param foldChangeLine Where to place a line for a fold change cutoff. Needs two values
#'                      (for positive andnegative fold change). Default is to plot no line
#' @param markGenes geneIDs for the genes to mark (a circle with gene name will be added)
#' @param colorGenes geneIDs for the genes to color
#' @param useGeneNames Use gene names instead of default geneIDs. only works if the input
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
                        colorGenes = NULL, useGeneNames = TRUE, outFile = NULL) {

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
        if(!is.null(outFile)) {
                pdf(outFile)
        }

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

        if(!is.null(outFile)) {
                dev.off()
        }
}



#' Plot heatmap of raw counts for top DEgenes using DESeq2 output
#'
#' @param DESeqOutput A tab-seperated DESeq2 output file
#' @param fdr FDR cutoff for DE genes
#' @param fcountOutput Featurecounts output (for raw counts)
#' @param sampleNames samplenames for heatmap column label (must be the same order as in featurecounts file)
#' @param topNgenes How many genes to plot. Type NULL for all genes (might take long time)
#' @param clusterbyCorr Whether to cluster row/columns by correlations. Default is eucledian distances.
#' @param useGeneNames Use gene names to plot instead of default geneIDs. only works if the input
#'                      has been annotated by \code{\link{annotate_DESeqOutput}} function.
#' @param markGenes Genes to mark in bold on heatmap.
#' @param outFile File name to save the output. NULL prints heatmap on screen.
#'
#' @return A heatmap.
#' @export
#'
#' @examples
#'
#'

plotHeatmap <- function(DESeqOutput, fdr = 0.05, fcountOutput, sampleNames, topNgenes = 100,
                        clusterbyCorr = FALSE, useGeneNames = TRUE, markGenes = NULL,
                        outFile = NULL) {

        ## read the data and subset it
        deseqRes <- read.delim(DESeqOutput,header = TRUE)
        deseqRes <- deseqRes[which(deseqRes$padj < fdr),]
        deseqRes <- deseqRes[order(deseqRes$padj),]
        deseqRes <- deseqRes[!duplicated(deseqRes$Row.names),]
        if(!is.null(topNgenes)) {
                deseqRes <- deseqRes[1:topNgenes,]
        }

        ## read and subset fcountoutput
        fcout <- read.delim(fcountOutput,skip = 1, header = TRUE, row.names = 1)
        fcout <- fcout[which(rownames(fcout) %in% deseqRes$Row.names), 6:ncol(fcout)]
        fcout <- fcout[order(deseqRes$padj),]
        colnames(fcout) <- sampleNames

        ## make heatmap
        # cluster rows by pearson and cols by spearman
        if(clusterbyCorr){
                hr <- as.dist(1-cor(t(as.matrix(fcout)), method="pearson"))
                hc <- as.dist(1-cor(as.matrix(fcout), method="spearman"))
        } else {
                hr <- NULL
                hc <- NULL
        }
        # use gene names
        if(useGeneNames){
                rowlab = deseqRes$external_gene_name
        } else rowlab = NULL

        pheatmap::pheatmap(fcout, clustering_distance_rows = hr, cutree_rows = splitClusters,
                           clustering_distance_cols = hc, labels_row = rowlab,
                           main = sprintf("Raw counts: Top %d DE genes",topNgenes), filename = outFile)
        # splitting clusters (not implemented)
        if(markGenes) print("marking genes not implemented yet!")
}



#' Plot Stacked barchart of DE genes using DESeq2 output
#'
#' @param DESeqOutput A tab-seperated DESeq2 output file
#' @param fdr FDR cutoff for DE genes
#' @param foldCh Which scale of fold-change to plot. Choose from "abs" (absolute)
#'                      and "log2" (log2).
#' @param sampleName Samplename to print on the plot.
#' @param outFile Output pdf file name. If not, plot will be printed on screen.
#'
#' @return Stacked bar chart of DE gene counts.
#' @export
#'
#' @examples
#'
#'

plotStackedBars <- function(DESeqOutput, fdr = 0.05, foldCh = "abs", sampleName = NULL, outFile = NULL) {

        deseqRes <- read.delim(DESeqOutput,header = TRUE)

        # subset by fdr
        deseqRes <- deseqRes[which(deseqRes$padj < fdr),]
        deseqRes <- na.omit(deseqRes[!duplicated(deseqRes$Row.names),])
        deseqRes$Status <- ifelse(deseqRes$log2FoldChange < 0, "Down", "Up")
        deseqRes$absfoldch <- 2^abs(deseqRes$log2FoldChange)

        # make fold ch category
        if(foldCh == "abs") {
                deseqRes$category <- ifelse(deseqRes$absfoldch < 2, "< 2 fold",
                                            ifelse(deseqRes$absfoldch < 6, "2 to 6 fold",
                                                   ifelse(deseqRes$absfoldch < 10, "6 to 10 fold", "> 10 fold")
                                            ))
                deseqRes$category <- factor(deseqRes$category,
                                            levels = c("< 2 fold", "2 to 6 fold","6 to 10 fold", "> 10 fold"))
        } else {
                deseqRes$category <- ifelse(abs(deseqRes$log2FoldChange) < 2, "< 2 fold",
                                            ifelse(abs(deseqRes$log2FoldChange) < 4, "2 to 4 fold",
                                                   ifelse(abs(deseqRes$log2FoldChange) < 8, "4 to 8 fold", "> 8 fold")
                                            ))
                deseqRes$category <- factor(deseqRes$category,
                                            levels = c("< 2 fold", "2 to 4 fold", "4 to 8 fold", "> 8 fold"))
        }

        # plot
        p <- ggplot(deseqRes,aes(Status, fill = Status, alpha = category)) +
                geom_bar(colour = "black") + theme_grey(base_size = 16) +
                scale_y_continuous(breaks = round(seq(0, nrow(deseqRes) + 500, by = 500)) ) +
                scale_fill_manual(values = c("darkred","darkgreen")) +
                labs(y = "No. of Genes", fill = "Level",
                     title = paste0(sampleName, " DE genes (divided by fold change)" ))

        if(!is.null(outFile)) {
                pdf(outFile)
                print(p)
                dev.off()
        } else {
                return(p)
        }
}
