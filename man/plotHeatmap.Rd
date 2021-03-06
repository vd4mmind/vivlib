% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_wrappers.R
\name{plotHeatmap}
\alias{plotHeatmap}
\title{Plot heatmap of raw counts for top DEgenes using edgeR/DESeq2 output}
\usage{
plotHeatmap(DEoutput, fdr = 0.05, fcountOutput, sampleNames,
  topNgenes = 100, clusterbyCorr = FALSE, useGeneNames = TRUE,
  markGenes = NULL, outFile = NA)
}
\arguments{
\item{DEoutput}{A tab-seperated edgeR/DESeq2 output file, using \code{\link{EdgeR_wrapper}} or
\code{\link{DESeq_wrapper}}}

\item{fdr}{FDR cutoff for DE genes}

\item{fcountOutput}{Featurecounts output (for raw counts)}

\item{sampleNames}{samplenames for heatmap column label (must be the same order as in featurecounts file)}

\item{topNgenes}{How many genes to plot. Type NULL for all genes (might take long time)}

\item{clusterbyCorr}{Whether to cluster row/columns by correlations. Default is eucledian distances.}

\item{useGeneNames}{Use gene names to plot instead of default geneIDs. only works if the input
has been annotated by \code{\link{annotate_DEoutput}} function.}

\item{markGenes}{Genes to mark in bold on heatmap.}

\item{outFile}{File name to save the output. NULL prints heatmap on screen.}
}
\value{
A heatmap.
}
\description{
Plot heatmap of raw counts for top DEgenes using edgeR/DESeq2 output
}
\examples{
deout <- system.file("extdata", "edgeR_output_annotated.tsv", package="vivlib")
fc <- system.file("extdata", "fcount_mouse.tsv", package="vivlib")
plotHeatmap(DEoutput = deout, fcountOutput = fc, sampleNames = rep(c("cnt","KD"), each = 3))

}

