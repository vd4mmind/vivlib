
#' Make Voom transformed input for BarCodePlot and ROAST test
#'
#' @param counts A tab-seperated file containing gene names followed by counts
#' @param design A tab-seperated file containing design information (colnames, condition).
#'                      colnames should correspond to columns in count file and condition could be
#'                      control/test or any set of factors.
#' @param bmGeneNames Optionally provide alternative gene symbols downloaded from biomart as a
#'                      tab-seperated file. The columns should be ("ensembl_gene_id","external_gene_id")
#' @param name File name to output filtering plots.
#'
#' @return A plot (as pdf) and a voom transformed output as list.
#' @export
#'
#' @examples
#'makeVoomInput(counts,design,bmGeneNames,name="name")
#'

makeVoomInput <- function(counts,design,bmGeneNames,name="name"){

  # get gene count data
  counts <- read.table(counts, header = T, row.names = 1)
  rownames(counts) = gsub('(ENS.*)\\.[0-9]*','\\1',rownames(counts))
  # filter by rowmean
  means = rowMeans(counts)
  counts = counts[which(means > 1),]

  # add gene names from biomart file
  
  bmGeneNames = read.table(bmGeneNames,sep="\t", header=TRUE, row.names=1)
  matchingIds = merge(counts, bmGeneNames,
                      by.x = 0,
                      by.y = "ensembl_gene_id",
                      all.x = TRUE)
  matchingIds = matchingIds[c("Row.names","external_gene_id")]

  # get design matrix for voom
  design <- read.table(design, header=T)
  design <- model.matrix(~ condition, design)

  # voom
  y <- limma::voom(counts,design)
  y$Gene = tolower(as.character(matchingIds$external_gene_id))
  fit <- limma::lmFit(y, design = design)
  fit <- limma::eBayes(fit)
  voomOutput <- list(y = y, fit = fit)

  # make plots of filtering
  pdf(paste0(name,"filtering_plots.pdf"))
  boxplot(log2(counts),notch = TRUE,col="steelblue",cex.names=0.5,main="Log Counts after filtering")
  textplot(paste0("Number of genes after filtering : ",nrow(counts)))
  hist(fit$t[,"conditiontreatment"],col="steelblue",xlab="Range of test-statistic",
       main="Distribution of test statistics post-filtering")
  dev.off()

  return(voomOutput)
}


#' Perform ROAST and make BarCodePlot using a voom transformed output and a microarry output
#'
#' @description This function was originally written for comparing a  voom-transformed mouse input
#' and a microarry (affy) input from human or mouse. That's why the option to provide a Human-Mouse ID map.
#' But if you don't want that, you can simply use it to compare any set of microarry outputs to any voom-transformed input.
#'
#' @param GSEfile A microarry data file. Filename should start from 'GSE'.
#'              The columns should have 'gname' and 'SPOT_ID'.
#' @param ourVoomFile Voom transformed output from \code{\link{makeVoomInput}}
#' @param batchAnalyse If there are multiple microarry files in the folder, give the folder path.
#' @param VoomInputName Provide a name for our Voom-transformed input data.
#' @param humanMouseNameMap A tab-seperated file with column ("Sym_Human") which provides alternative
#'                      human symbol for our mouse gene IDs
#' @param dfCutoff Which cutoff to use for filtering genes from GSEfile.
#'                      Default is by pvalue. Otherwise a log Fold Change cutoff can be used.
#' @param logFoldCh logFoldChange cutoff for filtering GSE file (if used)
#' @param padj adjusted p-value cutoff for filtering GSE file (if used)
#' @param outFolder Folder to write the output pdf file.
#'
#' @return Pdf file with barcode plots
#' @export
#'
#' @examples
#'plotBarCodes(GSEfile,ourVoomFile)
#'

plotBarCodes <- function(GSEfile,ourVoomFile, batchAnalyse=TRUE, VoomInputName=NULL,
                         humanMouseNameMap = "/data/akhtar/bhardwaj/my_annotations/Human_mouse_Gene_orthologous.txt",
                         dfCutoff="pvalue",logFoldCh = 0, padj = 0.05,
                         outFolder="/data/akhtar/bhardwaj/2015_OtherAnalysis/Bilal_Oncogene_paper/BarCodePlots/MPC5") {

        # one file or many?
        if(batchAnalyse == FALSE) {
                GSE_files <- GSEfile
                } else {
                        GSE_files <- list.files(path=GSEfile,pattern="GSE")
                }

  # seperate matrix and fit inputs
  y = ourVoomFile$y
  fit = ourVoomFile$fit

  for (d in GSE_files){
    GSEname = gsub('.*(GSE.*).txt','\\1',d)
        print(paste(VoomInputName,"vs",GSEname))

        ## parse GSE data files
        GSE_df <- read.table(d, header=T)
        colnames(GSE_df) <- gsub('SPOT_ID','gname',colnames(GSE_df))
        GSE_df$gname <-tolower(GSE_df$gname)

        if(dfCutoff== "pvalue"){
          ## take only UPs and Downs from Online Data
          print(paste0("Filtering by P-value:",padj))
          GSE_df <- subset(GSE_df, adj.P.Val < padj)

          GSE_df <- aggregate(logFC ~ gname, data=GSE_df, FUN=mean)
          GSE_df[which(GSE_df$logFC < 0),2] <- -1
          GSE_df[which(GSE_df$logFC > 0),2] <- 1

        } else {
          print(paste0("Filtering by logFC : +- ", logFoldCh))
          GSE_df <- aggregate(logFC ~ gname, data=GSE_df, FUN=mean)
          GSE_df <- subset(GSE_df,abs(logFC) > logFoldCh)

          GSE_df[which(GSE_df$logFC < 0),2] <- -1
          GSE_df[which(GSE_df$logFC > 0),2] <- 1
        }

        print(paste0("No. of selected genes in sample:",d," = ",nrow(GSE_df)))

        if(!is.null(humanMouseNameMap)){ # convert human gene names to mouse
          print("Converting Human IDs to mouse")
          geneMap <- read.table(humanMouseNameMap,sep="\t",header=T)

          geneMap$Sym_Human = tolower(geneMap$Sym_Human)
          GSE_df = merge(GSE_df,geneMap,by.x = "gname",by.y = "Sym_Human",all.x = TRUE)
        }

        # match our Gene names to the gene names in the GSE datasets
        GSE_df$gname = tolower(as.character(GSE_df$gname))
        GSE_df$idx <- match(GSE_df$gname, y$Gene)
        GSE_df <- na.omit(GSE_df) # remove NAs

        print("Running Test and making plots")
        pdf(file = paste0(outFolder,"/",GSEname,"_vs_",VoomInputName, ".pdf"), width = 6, height = 6)

        ## the stats[index] is based on the gene symbol
        limma::barcodeplot(statistics = fit$t[,"conditiontreatment"],index = GSE_df$idx[GSE_df$logFC > 0],
                                 index2 = GSE_df$idx[GSE_df$logFC < 0], main = (paste(GSEname," vs ",VoomInputName))
                    )

        ## Do roast and add testScores
        testRes <- roast(index = GSE_df$idx, y = y, design = y$design,contrast = 2,
                                   nrot = 9999, gene.weights = GSE_df$logFC)
        gplots::textplot(as.data.frame(testRes),cex=0.5,col.colnames="steelblue",col.rownames = "red")

        dev.off()
      }
  }
