
#' Cluster DE genes by fold change from multiple files
#'
#' @param DEoutList A vector, with names of DESeq2 output files to use.
#' @param sampleNames Name of samples corresponding to the DESeq output file list.
#' @param FDRcutoff FDR cutoff to select DE genes from the list
#' @param method Which clustering method to use. "correlation" or "biclustering",
#' will allow output of gene names per cluster. Other methods are also supported
#' (all methods of hclust + kmeans), but won't output genes per cluster
#' @param cut_cluster A number to which the cluster dendrogram will be cut into
#' (NA means do not cut clusters)
#'
#' @param row_annotation A data frame for the annotation of genes, with rownames
#' corresponding to the Row.names column of the DESeq2 output. The columns can have
#' annotations like chromosome, gene type etc..
#'
#' @param keepNAs Many genes will have fold change = NA in some samples after merging,
#' since they are undetected in some samples. Select this if you want to still keep
#' those genes (NAs will be converted to zeros for clustering).
#'
#' @param outFile_prefix A prefix for output files.
#'
#' @return A pdf with clustered heatmap, an .Rdata file with the hclust objects and the
#' genes sorted by clustered output, and text file with genes divided by clusters
#' (if cut_clusters is selected).
#'
#' @export
#'
#' @examples
#'
#' clusterDEgenes(DEoutList, sampleNames, FDRcutoff = 0.05, method = "correlation")
#'

clusterDEgenes <- function(DEoutList, sampleNames, FDRcutoff = 0.05, method = "correlation",
                           cut_cluster = NA, row_annotation = NULL, keepNAs = TRUE, outFile_prefix = NULL) {
        # Read files
        dedata <- lapply(DEoutList, function(x){
                read.delim(pipe(paste0("cut -f1,3,7 ",x)),stringsAsFactors = FALSE)
        })
        names(dedata) <- sampleNames
        # take cutoff and rename columns by sample
        dedata <- mapply(function(x,y){
                x <- x[which(x[,3] < FDRcutoff),1:2]
                colnames(x) <- c("GeneID",y)
                return(x)
        },dedata, sampleNames,SIMPLIFY = FALSE)
        # merge all DE genes
        dedata <- Reduce(function(x,y) merge(x, y, by = "GeneID", all.x = TRUE, all.y = TRUE), dedata)
        dedata <- unique(dedata)
        # many DEgenes will have FoldCh = NA in other samples. How many?
        nagenes <- nrow(dedata) - nrow(na.omit(dedata))
        print(paste0("No. of unique DE genes after merging : ", nrow(dedata)))
        print(paste0("No. of genes with FoldChange = NA in atleast one sample : ", nagenes))
        # get in shape
        rownames(dedata) <- dedata$GeneID
        dedata <- dedata[,2:ncol(dedata)]
        if(keepNAs){
                dedata[is.na(dedata)] <- 0
                message("NAs converted to zeros for clustering!")
        } else {
                dedata <- na.omit(dedata)
                message("Genes with FoldChange = NA removed for clustering!")
        }

        ## compute distance

        if(method == "correlation") {
                hr <- as.dist(1 - cor(t(as.matrix(dedata)), method = "spearman"))
                hc <- as.dist(1 - cor(as.matrix(dedata), method = "spearman"))
                k <- NA
        } else if(method == "bicluster") {
                print("Biclustering not implemented yet")
        } else if(method == "kmeans") {
        	    warning("Choosing kmeans! You won't get genes splitted by clusters..")
        	    hr <- "euclidean"
                hc <- "euclidean"
                k <- cut_cluster

        } else {

                hr <- dist(as.matrix(dedata), method = method)
                hc <- dist(t(as.matrix(dedata)), method = method)
                k <- NA
        }

        # make row annotation
        if(!is.null(row_annotation)){
        	rowannot <- row_annotation[which(rownames(row_annotation) %in% rownames(dedata)),]
        } else {
        	rowannot <- NA
        }

        ## plot clusters
        paletteLength <- 100
        colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))(paletteLength)
        # set breaks to get white to zero
        myBreaks <- c(seq(min(dedata), 0,
                          length.out=ceiling(paletteLength/2) + 1),
                      seq(max(dedata)/paletteLength, max(dedata),
                          length.out=floor(paletteLength/2)))

        if(!is.null(outFile_prefix)) pdf(paste0(outFile_prefix,"_heatmap.pdf"))
        lapply(c("row","none"), function(x){
                pheatmap::pheatmap(dedata,color = colors, breaks = myBreaks, kmeans_k = k,
                                   clustering_distance_rows = hr,clustering_distance_cols = hc,
                                   annotation_row = rowannot,
                                   show_rownames = FALSE, scale = x,cutree_rows = cut_cluster,
                                   main = paste0("DE genes clustered. Scale : ",x)
                )
        })
        if(!is.null(outFile_prefix)) dev.off()

        ## sort genes as per cluster
        if(class(hr) == "dist"){
        	message("Sorting genes as per clusters.")
        	hr.h <- hclust(hr)
        	hc.h <- hclust(hc)
        	sortedRes <- dedata[hr.h$labels[hr.h$order],hc.h$labels[hc.h$order]]
        }

        ## cut clusters

        if(!is.na(cut_cluster)) {
        	message(paste0("Cutting tree into ",cut_cluster," clusters"))
                cluscut <- dendextend::cutree(hclust(hr), k = cut_cluster, order_clusters_as_data = FALSE)

                declusters <- lapply(seq(1:cut_cluster), function(i){
                        clust_ids <- names(cluscut[cluscut == i])
                        return(rownames(dedata[clust_ids,]))
                })
                names(declusters) <- paste0("Cluster_",seq(1:cut_cluster))
                declusters <- plyr::ldply(declusters, data.frame)
                colnames(declusters) <- c("Cluster_ID","GeneID")
                write.table(declusters, file = paste0(outFile_prefix,"_clusters.txt"),
                            sep = "\t", quote = FALSE, row.names = FALSE )
        } else {
                declusters <- list()
        }
        save(dedata,sortedRes, hr, hc, file = paste0(outFile_prefix,"_clustering.Rdata") )
}



#' plot intersections of list of DE genes from multiple DESeq2 outputs
#'
#' @param DEoutList A vector, with names of DESeq2 output files to use.
#' @param sampleNames Name of samples corresponding to the DESeq output file list.
#' @param FDRcutoff FDR cutoff to select DE genes from the list
#' @param outFile Name of output pdf file. If NULL, prints output on the screen
#'
#' @return A plot of Intersection of gene sets.
#' @export
#'
#' @examples
#' plotDEgeneOverlap(DEoutList, sampleNames, FDRcutoff = 0.05, outFile = NULL)
#'

plotDEgeneOverlap <- function(DEoutList, sampleNames, FDRcutoff = 0.05, outFile = NULL){
        dedata <- lapply(DEoutList, function(x){
                read.delim(pipe(paste0("cut -f1,3,7 ",x)),stringsAsFactors = FALSE)
        })
        names(dedata) <- sampleNames
        # take cutoff and rename columns by sample
        dedata <- mapply(function(x,y){
                x <- x[which(x[,3] < FDRcutoff),c(1,3)]
                colnames(x) <- c("GeneID",y)
                return(x)
        },dedata, sampleNames,SIMPLIFY = FALSE)

        # merge all DE genes
        dedata <- Reduce(function(x,y) merge(x, y, by = "GeneID", all.x = TRUE, all.y = TRUE), dedata)
        dedata <- unique(dedata)
        # get in format for plot
        rownames(dedata) <- dedata$GeneID
        dedata <- dedata[2:length(dedata)]
        dedata[dedata < FDRcutoff] <- 1
        dedata[dedata < 1] <- 0
        dedata[is.na(dedata)] <- 0

        if(!is.null(outFile)) pdf(outFile,width=1200,height = 1000)
        UpSetR::upset(dedata,sets = colnames(dedata),number.angles = 30, point.size = 5,
              name.size = 16, line.size = 2)

        if(!is.null(outFile)) dev.off()


}
