
#' Cluster DE genes by fold change from multiple files
#'
#' @param DEoutList A vector, with names of DESeq2 output files to use.
#' @param sampleNames Name of samples corresponding to the DESeq output file list.
#' @param FDRcutoff FDR cutoff to select DE genes from the list
#' @param method Which clustering method to use. Either "correlation" or "biclustering"
#' @param cut_cluster A number to which the cluster dendrogram will be cut into (NA means do not cut clusters)
#' @param row_annotation A data frame for the annotation of genes, with rownames corresponding to the Row.names column of the DESeq2 output. The columns can have annotations like chromosome, gene type etc..
#' @param outFile_prefix A prefix for output files.
#'
#' @return A pdf with clustered heatmap, an .Rdata file with the hclust objects and the genes sorted by clustered output, and text file with genes divided by clusters (if cut_clusters is selected).
#' @export
#'
#' @examples
#'
#' clusterDEgenes(DEoutList, sampleNames, FDRcutoff = 0.05, method = "correlation")
#'

clusterDEgenes <- function(DEoutList, sampleNames, FDRcutoff = 0.05, method = "correlation",
                           cut_cluster = NA, row_annotation = NA, outFile_prefix = NULL) {
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
        dedata[is.na(dedata)] <- 0
        message("NAs converted to zeros for clustering!")

        ## compute distance
        if(method == "correlation") {
                hr <- as.dist(1 - cor(t(as.matrix(dedata)), method = "spearman"))
                hc <- as.dist(1 - cor(as.matrix(dedata), method = "spearman"))
        } else if(method == "bicluster") {
                print("Biclustering not implemented yet")
        } else {
                warning("Method neither correlation or bicluster! Picking eucledian..")
                hr <- NA
                hc <- NA
        }

        # make row annotation
        rowannot <- row_annotation[which(rownames(row_annotation) %in% rownames(dedata)),]
        ## plot clusters
        paletteLength <- 100
        colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))(pallength)
        # set breaks to get white to zero
        myBreaks <- c(seq(min(dedata), 0,
                          length.out=ceiling(paletteLength/2) + 1),
                      seq(max(dedata)/paletteLength, max(dedata),
                          length.out=floor(paletteLength/2)))

        if(!is.null(outFile)) pdf(paste0(outFile_prefix,"_heatmap.pdf"))
        lapply(c("row","none"), function(x){
                pheatmap::pheatmap(dedata,color = colors, breaks = myBreaks,
                                   clustering_distance_rows = hr,clustering_distance_cols = hc,
                                   clustering_method = "complete", annotation_row = rowannot,
                                   show_rownames = FALSE, scale = x,cutree_rows = cut_cluster,
                                   main = paste0("DE genes clustered. Scale : ",x)
                )
        })
        if(!is.null(outFile_prefix)) dev.off()

        ## sort genes as per cluster
        hr.h <- hclust(hr)
        hc.h <- hclust(hc)
        sortedRes <- dedata[rev(hr.h$labels[hr.h$order]),hc.h$labels[hc.h$order]]

        ## cut clusters
        message(paste0("Cutting tree into ",cut_cluster," clusters"))
        if(!is.na(cut_cluster)) {
                cluscut <- cutree(hclust(hr), k = cut_cluster)

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
