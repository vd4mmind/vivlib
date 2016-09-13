#' A wrapper over CAMERA test.
#'
#' @param fcountOutput A tab-seperated file containing gene names followed by fcountOutput
#' @param design  A tab-seperated file containing design information (colnames, condition).
#'                      colnames should correspond to columns in count file and condition could be
#'                      control/test or any set of factors.
#' @param bmGeneNames Optionally provide alternative gene symbols downloaded from biomart as a
#'                      tab-seperated file. The columns should be ("ensembl_gene_id","external_gene_id")
#' @param name File name to output filtering plots.
#' @param moduleFile A file with modules to test. Download from
#'                      \link{http://software.broadinstitute.org/gsea/msigdb/index.jsp}
#' @param outfileName Output pdf file name.
#'
#' @return CAMERA result object.
#' @export
#'
#' @examples
#'
#' \dontrun{
#' fc <- system.file("extdata", "fcount_mouse.out", package="vivlib")
#' design <- data.frame(rownames = c(paste0("cnt_", 1:3), paste0("KD_",1:3)),
#'                      condition = rep(c("cnt","KD"), each = 3) )
#'
#' runCamera(fcountOutput,design,bmGeneNames,name="myVoomInput",outfileName = "test")
#' }
#'
#'

runCamera <- function(fcountOutput,design,bmGeneNames, moduleFile="msigdb.v5.0.symbols.gmt", outfileName){


        # get count data and filter by rowmeans
        fcountOutput <- read.delim(fcountOutput, header = T, row.names = 1)
        fcountOutput <-  fcountOutput[,c(6:ncol(fcountOutput))]
        means = rowMeans(fcountOutput)
        fcountOutput = fcountOutput[which(means > 1),]
        # design matrix
        colnames(design) <- "condition"
        design <- model.matrix(~ condition, design)

        # add gene names from biomart file
        ext.data <- fetch_annotation(genome)
        fcountOutputMerged <- merge(fcountOutput, ext.data, by.x = 0, by.y = 1, all.x = TRUE, sort = FALSE)

        # voom transform
        y <- limma::voom(fcountOutput,design)
        y$Gene = tolower(as.character(fcountOutputMerged$external_gene_name))
        fit <- limma::lmFit(y, design = design)

        fit$df.residual <- limma::fit$df.residual - 1 # for CAMERA
        fit <- limma::eBayes(fit,trend = T) # for CAMERA
        fit$gene = tolower(as.character(fcountOutputMerged$external_gene_name))
        print(summary(limma::decideTests(fit)))

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
        index <- limma::symbols2indices(ldat,fit$gene,remove.empty = TRUE)
        gst <- limma::camera(index = index,y = y$E,design = design,allow.neg.cor=FALSE)
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
#' \dontrun{
#' Camera_plotbubble(CAMERAoutput, outfileName = "test", top = 20, title = NULL)
#' }
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
