#!/usr/bin/env Rscript

library(Seurat)
library(stringr)
library(tools)
library(optparse)
# source("/gpfs/oe-scrna/pipeline/scRNA-seq_further_analysis/function/seuratFc.r")
source("/home/liuchenglong/script/lclFunc.r")

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input file path (RDS file)"),
  make_option(c("-o", "--outdir"), type = "character", default = "./", help = "Output directory to save the CSV file"),
  make_option(c("-c", "--subnew_celltype"), type = "character", default = "all", help = "Subnew cell type"),
  make_option(c("-s", "--subsampleid"), type = "character", default = "all", help = "Subsample ID"),
  make_option(c("-g", "--subgroup"), type = "character", default = "all", help = "Subgroup"),
  make_option(c("-p", "--predicate"), type = "character", default = NULL, help = "Filtering predicate (e.g., '')"),
  make_option(c("-l", "--subcluster"), type = "character", default = "all", help = "Subcluster")
)

# Parse the command line arguments
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

rds <- readSeurat(args$input)
rds = subsetRDS(rds, args$subsampleid, args$subgroup, args$subcluster, args$subnew_celltype, args$predicate)

output_file <- paste0(args$outdir, "/counts.csv")
write.csv(t(as.matrix(rds@assays$RNA@counts)), file = output_file)
