#' A wrapper over CAMERA test.
#'
#' @param counts A tab-seperated file containing gene names followed by counts
#' @param design  A tab-seperated file containing design information (colnames, condition).
#'                      colnames should correspond to columns in count file and condition could be
#'                      control/test or any set of factors.
#' @param bmGeneNames Optionally provide alternative gene symbols downloaded from biomart as a
#'                      tab-seperated file. The columns should be ("ensembl_gene_id","external_gene_id")
#' @param name File name to output filtering plots.
#' @param moduleFile A file with modules to test. Download from
#'                      \code{\link{http://software.broadinstitute.org/gsea/msigdb/index.jsp}}
#' @param outfileName Output pdf file name.
#'
#' @return CAMERA result object.
#' @export
#'
#' @examples
#'runCamera(counts,design,bmGeneNames,name="myVoomInput",moduleFile="msigdb.v5.0.symbols.gmt")
#'

runCamera <- function(counts,design,bmGeneNames, moduleFile="msigdb.v5.0.symbols.gmt",  outfileName){

        design <- read.table(design, header=T)
        design <- model.matrix(~ condition, design)

        # get count data and filter by rowmeans
        counts <- read.table(counts, header = T, row.names = 1)
        rownames(counts) = gsub('(ENS.*)\\.[0-9]*','\\1',rownames(counts))
        means = rowMeans(counts)
        counts = counts[which(means > 1),]

        # add gene names from biomart file
        bmGeneNames = read.table(bmGeneNames,sep="\t", header=TRUE, row.names=1)
        matchingIds = merge(counts, bmGeneNames,
                            by.x = 0,
                            by.y = "ensembl_gene_id",
                            all.x = TRUE)
        matchingIds = matchingIds[c("Row.names","external_gene_id")]

        # voom transform
        y <- voom(counts,design)
        y$Gene = tolower(as.character(matchingIds$external_gene_id))
        fit <- lmFit(y, design = design)

        fit$df.residual <- fit$df.residual - 1 # for CAMERA
        fit <- eBayes(fit,trend = T) # for CAMERA
        fit$gene = tolower(as.character(matchingIds$external_gene_id))
        print(summary(decideTests(fit)))

        # read and prepare the ModuleFile
        Mods <- read.table(moduleFile,fill = TRUE, sep="\t")
        Mods <- Mods[,c(1,3:length(Mods))]
        nams <- Mods[,1]
        dat <- Mods[,-1]
        ldat <- split(dat,seq_len(nrow(dat)))
        ldat <- lapply(ldat,function(x) x[x != ""])
        names(ldat) <- nams
        ldat <- lapply(ldat,tolower)

        # Run CAMERA and print output
        index <- symbols2indices(ldat,fit$gene,remove.empty = TRUE)
        gst <- camera(index = index,y = y$E,design = design,allow.neg.cor=FALSE)
        gst <- gst[order(gst[,5],decreasing = FALSE),]

        print("No of Significant modules : ")
        print(table(gst[,5] < 0.05))
        print("List of significant modules : ")
        print(gst[which(gst$FDR < 0.05),])

        # write output
        write.table(gst,outFileName,sep="\t")

}

#' Make a bubble plot for CAMERA output
#'
#' @param CAMERAoutput A tab-seperated file with CAMERA results.
#' @param outfileName Output pdf file name.
#' @param top Number of top gene-sets to plot.
#' @param title Title of the plot
#'
#' @return A pdf file with bubble plot
#' @export
#'
#' @examples
#'
#' Camera_plotbubble(CAMERAoutput, outfileName, top = 20, title = NULL)
#'

Camera_plotbubble <- function(CAMERAoutput, outfileName = NULL, top = 20, title = NULL){
        cam.dat <- read.table(CAMERAoutput, sep="\t", row.names = NULL)
        cam.dat <- cam.dat[order(cam.dat$FDR),]
        cam.dat <- cam.dat[1:top,]
        cam.dat$FDR <- -log10(cam.dat$FDR)
        colnames(cam.dat) <- c("Pname",colnames(cam.dat[2:ncol(cam.dat)]))

        # plot and save
        p <- ggplot(camtop,aes(Correlation, FDR, fill = Direction, size = NGenes, label = Pname)) +
                ggplot2::geom_point(alpha=0.7,shape=21) +
                ggplot2::geom_text(size=4) +
                ggplot2::scale_size_area(max_size = 15) +
                ggplot2::theme_bw(base_size = 15) +
                ggplot2::labs(x = "Inter-gene Correlation", y = "-log10(p-value)",fill="Category")
        if (!is.null(title)) {
                p <- p + ggplot2::ggtitle(title)
        }

        if(!is.null(outfileName)) pdf(outfileName)
        print(p)
        if(!is.null(outfileName)) dev.off()
}
