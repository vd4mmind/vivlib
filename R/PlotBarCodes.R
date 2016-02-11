
# The script can be used to compare a count dataset to online microarrary datasets, 
# with columns "gname/SPOTID" (for gene name) & "adj.P.Val". 
# Author: Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de) (11/06/2015)

library(limma)
library(gplots)

makeVoomInput <- function(counts,design,bmGeneNames,name="name"){
  design <- read.table(design, header=T)
  design <- model.matrix(~ condition, design)
  counts <- read.table(counts, header = T)
  rownames(counts) = gsub('(ENS.*)\\.[0-9]*','\\1',rownames(counts))
  means = rowMeans(counts)
  counts = counts[which(means > 1),]
  bmGeneNames = read.table(bmGeneNames,sep="\t", header=TRUE, row.names=1)
  matchingIds = merge(counts, bmGeneNames,
                      by.x = 0,
                      by.y = "ensembl_gene_id",
                      all.x = TRUE)
  matchingIds = matchingIds[c("Row.names","external_gene_id")]
  y <- voom(counts,design)
  y$Gene = tolower(as.character(matchingIds$external_gene_id))
  fit <- lmFit(y, design = design)
  fit <- eBayes(fit)
  voomInput <- list(y = y, fit = fit)
  
  # make plots of filtering
  pdf(paste0(name,"filtering_plots.pdf"))
  boxplot(log2(counts),notch = TRUE,col="steelblue",cex.names=0.5,main="Log Counts after filtering")
  textplot(paste0("Number of genes after filtering : ",nrow(counts)))
  hist(fit$t[,"conditiontreatment"],col="steelblue",xlab="Range of test-statistic",
       main="Distribution of test statistics post-filtering")
  dev.off()
  
  return(voomInput)
}


plotBarCodes <- function(GSEfile,ourVoomFile,batchAnalyse=TRUE,isHuman=FALSE,name=NULL,
                         humanMouseNameMap = "/data/akhtar/bhardwaj/my_annotations/Human_mouse_Gene_orthologous.txt",
                         dfCutoff="pvalue",logFoldCh = 0, padj = 0.05,
                         outFolder="/data/akhtar/bhardwaj/2015_OtherAnalysis/Bilal_Oncogene_paper/BarCodePlots/MPC5"){
                         
  if(batchAnalyse == FALSE){
    GSE_files <- GSEfile
  } else {
    GSE_files <- list.files(path=GSEfile,pattern="GSE")
  }
  # seperate matrix and fit inputs
  y = ourVoomFile$y
  fit = ourVoomFile$fit
    
  for (d in GSE_files){
    GSEname = gsub('.*(GSE.*).txt','\\1',d)
        print(paste(name,"vs",GSEname))
        ## parse GSE data files
        GSE_df <- read.table(d, header=T)
        colnames(GSE_df) <- gsub('SPOT_ID','gname',colnames(GSE_df))
        GSE_df$gname <-tolower(GSE_df$gname)
        
        if(dfCutoff== "pvalue"){
          ## take only UPs and Downs from Online Data
          print(paste0("Filtering by P-value:",padj," and logFC : ", logFoldCh))
          GSE_df <- subset(GSE_df, adj.P.Val < padj)
          GSE_df <- aggregate(logFC ~ gname, data=GSE_df, FUN=mean)
          GSE_df[which(GSE_df$logFC < 0),2] <- -1
          GSE_df[which(GSE_df$logFC > 0),2] <- 1
          
        } else {
          print(paste0("Filtering only by logFC : +- ", logFoldCh))
          GSE_df <- aggregate(logFC ~ gname, data=GSE_df, FUN=mean)
          GSE_df <- subset(GSE_df,abs(logFC) > logFoldCh)
          GSE_df[which(GSE_df$logFC < 0),2] <- -1
          GSE_df[which(GSE_df$logFC > 0),2] <- 1
        }
        
        print(paste0("No. of diffExp genes in sample:",d," = ",nrow(GSE_df)))
           
        if(isHuman == TRUE){ # convert human gene names to mouse
          print("Converting Human IDs to mouse")
          geneMap <- read.table(humanMouseNameMap,sep="\t",header=T)
          geneMap$Sym_Human = tolower(geneMap$Sym_Human)
          GSE_df = merge(GSE_df,geneMap,by.x = "gname",by.y = "Sym_Human",all.x = TRUE)
        }
        
        GSE_df$gname = tolower(as.character(GSE_df$gname))
        GSE_df$idx <- match(GSE_df$gname, y$Gene)## match our Gene names to the gene names in the GSE datasets
        GSE_df <- na.omit(GSE_df)## remove NAs
        
        print("Running Test and making plots")
        pdf(file = paste0(outFolder,"/",GSEname,"_vs_",name, ".pdf"), width = 6, height = 6)
        ## the stats[index] is based on the gene symbol
        barcodeplot(statistics = fit$t[,"conditiontreatment"],index = GSE_df$idx[GSE_df$logFC > 0],
                                 index2 = GSE_df$idx[GSE_df$logFC < 0],main = (paste(GSEname," vs ",name)))
        ## Add testScores
        testRes <- roast(index = GSE_df$idx, y = y, design = y$design,contrast = 2,
                                   nrot = 9999, gene.weights = GSE_df$logFC)
        textplot(as.data.frame(testRes),cex=0.5,col.colnames="steelblue",col.rownames = "red")
        dev.off()
      }
  }  



runCamera <- function(counts,design,bmGeneNames,name="name",moduleFile="msigdb.v5.0.symbols.gmt"){
        design <- read.table(design, header=T)
        design <- model.matrix(~ condition, design)
        counts <- read.table(counts, header = T)
        rownames(counts) = gsub('(ENS.*)\\.[0-9]*','\\1',rownames(counts))
        means = rowMeans(counts)
        counts = counts[which(means > 1),]
        bmGeneNames = read.table(bmGeneNames,sep="\t", header=TRUE, row.names=1)
        matchingIds = merge(counts, bmGeneNames,
                            by.x = 0,
                            by.y = "ensembl_gene_id",
                            all.x = TRUE)
        matchingIds = matchingIds[c("Row.names","external_gene_id")]
        y <- voom(counts,design)
        y$Gene = tolower(as.character(matchingIds$external_gene_id))
        fit <- lmFit(y, design = design)
        fit$df.residual <- fit$df.residual - 1 # for CAMERA
        fit <- eBayes(fit,trend = T) # for CAMERA
        fit$gene = tolower(as.character(matchingIds$external_gene_id))
        print(summary(decideTests(fit)))
        # read and prepare ModuleFile
        Mods <- read.table(moduleFile,fill = TRUE, sep="\t")
        Mods <- Mods[,c(1,3:length(Mods))]
        nams <- Mods[,1]
        dat <- Mods[,-1]
        ldat <- split(dat,seq_len(nrow(dat)))
        ldat <- lapply(ldat,function(x) x[x != ""])
        names(ldat) <- nams
        ldat <- lapply(ldat,tolower)
        index <- symbols2indices(ldat,fit$gene,remove.empty = TRUE)
        # Run CAMERA
        gst<-camera(index = index,y = y$E,design = design,allow.neg.cor=FALSE)
        gst<-gst[order(gst[,5],decreasing = FALSE),]
        
        print("No of Significant modules : ")
        print(table(gst[,5]<0.05))
        print("List of significant modules : ")
        print(gst[which(gst$FDR < 0.05),])
        
        return(gst)
}
