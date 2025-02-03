library(stringr)
library(optparse)
source("/gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/scripts/script/netplot/utils/func.R")

# 参数解析配置
option_list <- list(
  make_option(c("--regulons"), type="character", help="scenic_regulons_importance.tsv file path"),
  make_option(c("--TF"), type="character", help="interested TF (comma-separated)"),
  make_option(c("--highGene"), type="character", default="F", help="interested Gene (comma-separated)"),
  make_option(c("--outdir"), type="character", help="the output dir"),
  make_option(c("--top"), type="integer", default=30, help="top N genes [default %default]"),
  make_option(c("--showlabel"), type="double", default=0.3, help="label display threshold [default %default]"),
  make_option(c("--highGeneSize"), type="double", default=3, help="highlight gene size [default %default]"),
  make_option(c("--highGeneColor"), type="character", default="00AA55", help="highlight gene color [default %default]"),
  make_option(c("--TFColor"), type="character", default="1E90FF", help="TF color [default %default]"),
  make_option(c("--geneColor"), type="character", default="E87D72", help="gene color [default %default]")
)

# 解析参数
opt <- parse_args(OptionParser(option_list=option_list))

# 参数处理
int_tfs <- unlist(strsplit(opt$TF, ","))
high_genes <- if(opt$highGene == "F") NULL else unlist(strsplit(opt$highGene, ","))

# 创建目录
dir.create(opt$outdir, showWarnings = FALSE)
setwd(opt$outdir)

# 核心处理函数
process_tf <- function(tf, df, top) {
  subset_df <- df[df$TF == tf, ]
  ordered_df <- subset_df[order(-subset_df$Weight), ]
  ordered_df$topNub <- seq_len(nrow(ordered_df))
  head(ordered_df, top)
}

# 主数据处理
regulon_data <- read.delim(opt$regulons, stringsAsFactors = FALSE)
edge_list <- lapply(int_tfs, process_tf, df = regulon_data, top = opt$top)
edge_df <- do.call(rbind, edge_list)

# 高亮基因处理
if(!is.null(high_genes)) {
  find_highlight <- function(tf_df) {
    matches <- tf_df[tf_df$Gene %in% high_genes, ]
    if(nrow(matches) > 0) matches[which.max(matches$Weight), ]
  }
  highlight_edges <- lapply(split(edge_df, edge_df$TF), find_highlight)
  edge_df <- unique(rbind(edge_df, do.call(rbind, highlight_edges)))
}

# 顶点数据构建
all_genes <- unique(c(edge_df$Gene, edge_df$TF))
show_labels <- c(int_tfs, 
                if(!is.null(high_genes)) high_genes, 
                edge_df$Gene[seq_len(round(opt$top * opt$showlabel))])

vertex_df <- data.frame(
  name = all_genes,
  label = ifelse(all_genes %in% show_labels, all_genes, ""),
  size = ifelse(all_genes %in% int_tfs, 10,
               ifelse(all_genes %in% high_genes, opt$highGeneSize,
                     ifelse(all_genes %in% show_labels, 3, 1))),
  color = ifelse(all_genes %in% int_tfs, paste0("#", opt$TFColor),
                ifelse(all_genes %in% high_genes, paste0("#", opt$highGeneColor),
                      ifelse(all_genes %in% show_labels, paste0("#", opt$geneColor), "#FFFFFF")))
)

# 结果输出
write.csv(edge_df, "edge.csv", row.names = FALSE)
write.csv(vertex_df, "vertex.csv", row.names = FALSE)