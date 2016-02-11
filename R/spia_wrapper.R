

#' A wrapper over SPIA pathway analysis for mouse and human
#'
#' @param DESeqOut Tab-seperated DESeq2 output
#' @param organism select from "mmu" (mouse) or "hsa" (human)
#' @param padjCutoff FDR cutoff for the genes to include in the list
#' @param outFile output file to write affected pathways
#' @param outPlot output SPIA two-way evidence plot + output summary (pdf)
#'
#' @return a txt file and a pdf file
#' @export
#'
#' @examples
#'
#'
spia_wrapper <- function(DESeqOut, organism = "mmu", padjCutoff, outFile, outPlot){
        df <- read.table(DESeqOut, sep= "\t", header = TRUE, quote = "" )
        # load corresponding lib
        orgmap <- data.frame(org = c("mmu","hsa"), db = c("org.Mm.eg.db","org.Hs.eg.db"),stringsAsFactors = FALSE)
        assertthat::assert_that(organism %in% orgmap$org)
        db <- orgmap[grep(organism,orgmap$org),2]
        library(db,character.only = TRUE)

        # get ENTREZ id
        ensToEntrez <- AnnotationDbi::select(org.Mm.eg.db,as.character(df$Row.names),"ENTREZID", keytype = "ENSEMBL" )
        df <- merge(df,ensToEntrez,by = 1)
        df.dg <- dplyr::filter(df, padj < padjCutoff,!(is.na(ENTREZID)),!(duplicated(ENTREZID)))
        df.map <- df.dg$log2FoldChange
        names(df.map) <- df.dg$ENTREZID
        allgenes <- na.omit(as.character(df$ENTREZID)) # can take from any df. it's the universe

        # SPIA
        spia.degenes <- SPIA::spia(df.map,allgenes,organism = organism, nB = 2000)
        spia.degenes$Name <- substr(spia.degenes$Name,1,20)
        spia.degenes = spia.degenes[order(spia.degenes$pGFWER),]
        top <- head(spia.degenes[c(1:5,7,9,11)])
        colnames(top) <- c("TOP PATHWAYS : NAME",colnames(top[2:7]))

        # Write output
        pdf(outPlot)
        SPIA::plotP(spia.degenes)
        gplots::textplot(top)
        dev.off()
        write.table(spia.degenes,outFile, sep = "\t", quote = FALSE, row.names = FALSE)
}


#' Make a bubblePlot for SPIA pathway output
#'
#' @param SPIAout output from spia_wrapper
#' @param outfileName output pdf file name for plot
#' @param top How many top pathways to plot (by pGFdr value)
#' @param title Title of the plot
#'
#' @return plot A bubbleplot in pdf format
#' @export
#'
#' @examples
#' spia_plotBubble(spia_wrapper_output, outfileName = "test.out, top = 20, title = "test plot)
#'

spia_plotBubble <- function(SPIAout,outfileName, top = 20, title = NULL){
        path <- read.delim(pipe(paste0("cut -f1,3,9,11 ",SPIAout)), header = TRUE)
        path$pGFdr <- -log10(path$pGFdr)
        path <- path[order(path$pGFdr,decreasing = TRUE),]
        path <- path[1:top,]

        pdf(outfileName)
        print(ggplot2::ggplot(path,aes(Name,pGFdr,fill = Status, size = pSize)) +
                      ggplot2::geom_point(alpha = 0.7, shape = 21) +
                      ggplot2::scale_size_area(max_size = 15) + theme_bw(base_size = 15) +
                      ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
                      ggplot2::scale_fill_manual(values = c("forestgreen","Red")) +
                      ggplot2::labs(x = "Pathway size", y = "-log10(p-value)",fill = "Status", title = title)
        )
        dev.off()
}
