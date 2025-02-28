#!/usr/bin/env Rscript

######################################################################################
# Functions
######################################################################################

# Perform analysis for a given combination
perform_analysis <- function(dds, metadata, combination, min_padj, min_lfc, dir_name) {
  log_info("Performing analysis for: {combination}")
  res <- results(dds, name=combination)

  # Order the results based on padj and LFC
  res_ordered <- res[order(res$padj, -res$log2FoldChange), ]
  # Create file name to store the results
  file_name <- paste(combination, "_all_results.csv", sep="")
  file_name <- file.path(dir_name, file_name)
  write.csv(as.data.frame(res_ordered), 
            file=file_name,
            row.names=TRUE)

  log_info("Obtained results for {combination}. Here is the summary:")
  summary(res)

  # Create plots and save to pdf file
  # Create file name for the pdf
  pdf_file_name <- paste(combination, "_plots.pdf", sep="")
  pdf_file_name <- file.path(dir_name, pdf_file_name)
  pdf(pdf_file_name)

  log_info("Creating heatmap of sample-to-sample distances using the rLog transformed data.")
  # Create heatmap of sample-to-sample distances using the rlog transformed data
  condition <- strsplit(combination, "_")[[1]][1]
  (mycols <- 1:length(unique(metadata[[condition]])))
  (mycols <- brewer.pal(8, "Dark2")[1:length(unique(metadata[condition]))])
  sampleDists <- as.matrix(dist(t(assay(rld))))
  heatmap.2(as.matrix(sampleDists), key=T, trace="none",
            col=colorpanel(100, "purple", "white"),
            ColSideColors = mycols[metadata[[condition]]], RowSideColors = mycols[metadata[[condition]]],
            margin = c(10,10), main = "Sample Distance Matrix")
  

  log_info("Performing log fold shrinkage using apeglm.")
  # plot log fold shrinkage
  res_lfcs <- lfcShrink(dds, coef=combination, type = "apeglm")
  plotMA(res_lfcs)

  # Close pdf
  dev.off()

  # write csv file with only the significant results
  resdata <- merge(as.data.frame(res_lfcs), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "miRNA"

  res_table_thresh <- resdata %>%
    mutate(threshold=padj<min_padj & abs(log2FoldChange)>min_lfc)
  
  resdata_file_name = paste(combination, "_only_significant_results.csv", sep="")
  resdata_file_name = file.path(dir_name, resdata_file_name)
  write.csv(res_table_thresh, resdata_file_name, row.names=FALSE)

  if (length(res_table_thresh) == 0) {
    log_info("No significant results")
    dev.off()
    return
  } else if (length(res_table_thresh) < 10) {
    top_mirnas <- res_table_thresh
    log_info("Less than 10 significant results")
  } else {
    top_mirnas <- res_table_thresh %>%
      arrange(desc(abs(log2FoldChange))) %>%
      head(10)
      log_info("At least 10 significant results")
  }

  log_info("Making enchanced volcano plot")
  
  # Create Enhanced Volcano Plot
  volcano_plot <- EnhancedVolcano(res_lfcs,
                                  lab = row.names(res_lfcs),
                                  pCutoff = min_padj,
                                  FCcutoff = min_lfc,
                                  x = "log2FoldChange",
                                  y = "padj",
                                  selectLab = top_mirnas$miRNA,
                                  encircleCol = 'black',
                                  encircleSize = 3,
                                  encircleFill = 'green',
                                  encircleAlpha = 1/5,
                                  gridlines.major = FALSE,
                                  gridlines.minor = FALSE,
                                  labSize = 3,
                                  drawConnectors = TRUE,
                                  boxedLabels = TRUE)
  # Save volcano plot to png
  ggsave(paste(dir_name, "/", combination, "_volcano_plot.png", sep=""), plot=volcano_plot, width=8, height=6, dpi=600)
  log_info("Completed analysis of {combination}")


}

# Function to make a dataframe with only the samples in the metadata
get_metadata_samples <- function(counts, metadata){
  # Get the sample ids from the metadata
  sample_ids <- rownames(metadata)
  log_info("Sample ids from metadata: {sample_ids}")

  log_info("Sample ids from counts: {colnames(counts)}")
  
  # Get the read counts for samples within the metadata
  counts <- counts[, sample_ids]
  
  return(counts)
}

######################################################################################
# Get inputs
######################################################################################
library(argparse)
parser <- ArgumentParser()

parser$add_argument("--all_raw_counts", help="Path to the raw counts file for all samples", type="character")
parser$add_argument("--metadata", help="Path to the metadata file", type="character")
parser$add_argument("--min_padj", help="Minimum adjusted p-value for differential expression", type="numeric", default=0.05)
parser$add_argument("--min_log2fc", help="Minimum log2 fold change for differential expression", type="numeric", default=1.5)
parser$add_argument("--meta2_condition", help="Essentially the symbol for the metadata file", type="character", default="metadata")

args <- parser$parse_args()

all_raw_counts_file <- args$all_raw_counts
metadata_file <- args$metadata
min_padj <- args$min_padj
min_lfc <- args$min_log2fc
meta2_condition <- args$meta2_condition

######################################################################################
# Load libraries
######################################################################################

library(DESeq2)
library(ggplot2)
library(logger)
library(dplyr)
library(png)
library(RColorBrewer)
library(apeglm)
library(gplots)
library(ggrepel)
library(EnhancedVolcano)

######################################################################################
# Start preprocessing of data
######################################################################################

log_threshold(INFO)
log_formatter(formatter_paste)
log_appender(appender_file(paste(meta2_condition, "_DESeq2.log")))


log_info("Starting DESeq2 analysis")

# Get the read counts for the miRNAs in the samples
counts <- as.matrix(read.csv(all_raw_counts_file, row.names=1, sep='\t'))

log_info("Read counts for all samples.")
log_info("raw_counts_columns: {colnames(counts)}")

# Remove the rows from HTSeq that start with __ (just some extra information)
counts <- counts[!grepl("^__", rownames(counts)), ]

# Get the metadata
metadata <- read.csv(metadata_file, row.names=1)

log_info("Metadata for all samples. {colnames(metadata)}")

# Get the read counts for samples within the metadata
counts <- get_metadata_samples(counts, metadata)

# Verify that the columns in the counts data are the same as the rows in the metadata
all(colnames(counts) %in% rownames(metadata))
all(colnames(counts) == rownames(metadata))

log_info("Read counts for samples within the metadata:")
head(counts)

######################################################################################
# Start differential expression analysis
######################################################################################

# Create the design_string
for_design <- colnames(metadata)
if (length(for_design) == 1) {
  design_string <- paste("~", for_design[1])
} else {
  design_string <- paste("~", paste(for_design, collapse = " + "))
} 

log_info("Design string: {design_string}")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = metadata, 
                              design = as.formula(design_string))

log_info("colData(dds): {colData(dds)}")

# Run DESeq2 analysis
dds <- DESeq(dds)

# Print some diagnostic plots: dispersion plot and PCA plot
basic_plots_file_name = paste(meta2_condition, "_basic_plots_design_", paste(unlist(colnames(metadata)), collapse = "_"), ".pdf", sep="")
pdf(basic_plots_file_name)

# Dispersion plot
dispersion_plot_title = paste("Dispersion Plot for Design: ", design_string)
plotDispEsts(dds, main=dispersion_plot_title)

# PCA plot
# perform regularized log transformation
rld <- rlogTransformation(dds)
pca_plot_title = paste("PCA Plot for Design: ", design_string)
plotPCA(rld, intgroup=unlist(colnames(metadata)))
title(pca_plot_title)

dev.off()

# Iterate over all possible combinations available
for (combination in resultsNames(dds)) {
  if (combination == "Intercept") {
    next
  }

  # Create a folder for this combination
  main_dir <- paste(meta2_condition, "_deseq2_results")
  dir_name <- file.path(main_dir, combination)
  dir.create(dir_name, recursive=TRUE, showWarnings = FALSE)

  # Perform the analysis
  perform_analysis(dds, metadata, combination, min_padj, min_lfc, dir_name)
}


######################################################################################
# Save R sesions info
######################################################################################

log_info("Saving R session info")

sessionInfo_file_name = paste(meta2_condition, "_R_sessionInfo.log")
sink(sessionInfo_file_name)
sessionInfo()
sink()

log_info("Completed DESeq2 analysis")