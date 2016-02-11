#' Get introns by genes
#'
#' @param txdb A TxDb object
#' @param seqout The annotated DESeq Output
#' @param padjFilter p-value cutoff
#' @param genome_mart the name of mart object to annotate the input with (eg : hsapiens_gene_ensembl)
#' @param outfile output file with introns, to write back
#'
#' @return annotated DESeq output with introns for each gene
#' @export
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'
#' getIntrons_byGenes(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, seqout = "DESeq_outputs/XR1_annotated.out",
#' padjFilter = 0.05, genome_mart = "hsapiens_gene_ensembl", outfile = "annotated_introns.out")
#'

getIntrons_byGenes <- function(txdb, seqout = "DESeq_outputs/XR1_annotated.out",
                               padjFilter = 0.05, genome_mart = "hsapiens_gene_ensembl", outfile){

        ## Get introns for all genes in Txdb
        message("Creating introns from Txdb")
        introns <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE) # use ucsc txnames
        ulst <- unlist(introns)
        intronsNoDups <- ulst[!duplicated(ulst)]
        txnames <- names(intronsNoDups)
        map <- select(txdb, keys=txnames, keytype='TXNAME', cols='GENEID')
        idx <- map$GENEID[!is.na(map$GENEID)]
        intronsByGene <- split(intronsNoDups[!is.na(map$GENEID)], idx)
        names(intronsByGene) <- unique(idx)

        # convert GRL to dataframe
        gr <- unlist(intronsByGene)
        df <- data.frame(seqnames=seqnames(gr),
                         starts=start(gr)-1,
                         ends=end(gr),names = names(gr),
                         scores=c(rep(".", length(gr))), strands = strand(gr)
        )
        df <- tidyr::separate(df,names,c("refgeneID","ucscID","num"),by = "//.")

        ## Read deseq output and get entrez ids for genes
        message("fetching data for the input genes")
        seqdat <- read.table(seqout,sep = "\t",header = T,quote = "",stringsAsFactors = FALSE)
        seqdat <- as.data.frame(dplyr::filter(seqdat, padj < padjFilter))
        mart <- biomaRt::useMart("ensembl",dataset = genome_mart)
        martdata <- biomaRt::getBM(attributes = c("ensembl_gene_id","entrezgene"),mart = mart)
        martgenes <- martdata[which(martdata$ensembl_gene_id %in% seqdat$Row.names),2]
        if(length(na.omit(martgenes) < length(martgenes))){
                warning("Sorry! you lost some genes due to ID conversion issues")
        } else {
                message("Hurrey! No genes lost due to ID conversion issues")
        }

        ## merge the dfs and write back
        message("Merging and writing")
        seqdat <- merge(seqdat,martdata, by = 1, all.x = TRUE)
        seqdat <- merge(seqdat,df, by.x = "entrezgene", by.y = "refgeneID", all.x = TRUE)
        write.table(seqdat, outfile, quote = FALSE,sep = "\t",row.names = FALSE)
}
