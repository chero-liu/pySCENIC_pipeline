#!/usr/bin/env Rscript
source("/home/liuchenglong/script/lclFunc.r")
source("/home/liuchenglong/script/Get_colors.R")
# source('/gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/scripts/tools/utils.r')

library(Seurat)
library(tidyr)
library(dplyr)
library(pheatmap)
library(optparse)
library(SCENIC,lib.loc = '/gpfs/oe-software/scrna_software_bak/conda_env/Scenic/lib/R/library')
library(AUCell,lib.loc = '/gpfs/oe-software/scrna_software_bak/conda_env/Scenic/lib/R/library')
library(SCopeLoomR,lib.loc = '/gpfs/oe-software/scrna_software_bak/conda_env/Scenic/lib/R/library')

# 自定义函数模块 ---------------------------------------------------------------
handle_parameters <- function() {
    option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Path to the pyscenic sce_SCENIC.loom file path", metavar = "FILE"),
    make_option(c("-r", "--rds_filepath"), type = "character", help = "Path to input RDS file", metavar = "FILE"),
    make_option(c("-o", "--outdir"), type = "character", default = 'auto', help = "Output directory to save the CSV file"),
    make_option(c("-c", "--subnew_celltype"), type = "character", default = "all", help = "Subnew cell type"),
    make_option(c("-s", "--subsampleid"), type = "character", default = "all", help = "Subsample ID"),
    make_option(c("-g", "--subgroup"), type = "character", default = "all", help = "Subgroup"),
    make_option(c("-l", "--subcluster"), type = "character", default = "all", help = "Subcluster"),
    make_option(c("-p", "--predicate"), type = "character", default = NULL, help = "Filtering predicate (e.g.,"),
    make_option(c("-b", "--binmethod"), type = "character", default = "aucell", help = "Binning method"),
    make_option(c("--groupby"), type = "character", default = "group_new_celltype", help = "Grouping variable"),
    make_option(c("-t", "--threshold"), type = "numeric", default = 0, help = "Threshold value"),
    make_option(c("--topGenes"), type = "integer", default = 4, help = "Number of top genes to consider"),
    make_option(c("-e", "--extended"), action = "store_true", default = FALSE, help = "Whether to use extended mode"),
    make_option(c("--nclust"), type = "integer", default = 4, help = "the number of csi modules to use"),
    make_option(c("--utils_path"), type = "character", help = "")
    )
  parse_args(OptionParser(option_list = option_list))
}

analyze_regulon_activity <- function(seurat_obj, regulonAUC, args) {

  # 注释处理
  cell_info <- seurat_obj@meta.data %>% arrange(!!sym(args$groupby))

  # 热图生成
  color_palette <- colorRampPalette(c("#406AA8", "white", "#D91216"))
  annotation_colors <- list(
    setNames(getColorOrder(seurat_obj, args$groupby), 
            levels(seurat_obj@meta.data[[args$groupby]]))
  )
  names(annotation_colors) <- args$groupby

  heatmap_plot <- pheatmap(
    mat = regulonAUC[, rownames(cell_info)],
    scale = "row",
    cluster_cols = FALSE,
    show_colnames = FALSE,
    cluster_rows = FALSE,
    color = color_palette(200),
    annotation_col = cell_info[, args$groupby, drop = FALSE],
    annotation_colors = annotation_colors,
    treeheight_col = 10,
    border_color = NA,
    breaks = unique(c(seq(-2.5,0, length=100),  seq(0,2.5, length=100))),
    fontsize_row = 6
  )

  save_plot(heatmap_plot, 
              file.path(args$outdir, "1.1.regulon_activity_heatmap_groupby_cells"),
              width = 8, height = 6)

  # activity处理
  regulon_activity <- process_regulon_activity(
    cell_info, regulonAUC, args$groupby, nmads = 3
  )

  write.table(
    regulon_activity,
    file.path(args$outdir, "1.2.centered_regulon_activity_groupby_design.xls"),
    sep = "\t", row.names = TRUE, quote = FALSE
  )

  if (dim(regulon_activity)[1]<11) {
    hig <- 5
  } else {
    hig <- dim(regulon_activity)[1]*0.45
  }

  p_mean_pheatmap = pheatmap(regulon_activity,
                    cellwidth = 18,
                    cellheight = 18,
                    color=color_palette(200),
                    angle_col = 45,
                    treeheight_col=20, 
                    treeheight_row=20,
                    border_color=NA)

  save_plot(p_mean_pheatmap, 
              file.path(args$outdir, "1.3.regulon_activity_heatmap"),
              width = dim(regulon_activity)[2]/2+3,height = hig)

}

#' RSS分析模块
perform_rss_analysis <- function(seurat_obj, regulonAUC, args) {
  rss_matrix <- calcRSS(regulonAUC, seurat_obj@meta.data[[args$groupby]])
  rss_df <- na.omit(rss_matrix) %>%
    as.data.frame()
  rss_df$regulon = rownames(rss_df)
  rss_df = gather(rss_df,key = !!sym(args$groupby), value = "RSS", -regulon)

  # 可视化
  rss_plot <- RSSRanking(
    rss_df, 
    group.by = args$groupby,
    top_genes = args$topGenes,
    plot_extended = args$extended
  )

  write.table(rss_df, file.path(args$outdir, "2.1.regulon_RSS_annotation.xls"), sep = "\t", col.names = T, row.names =F, quote =F)
  save_plot(rss_plot, 
              file.path(args$outdir, "2.2.RSS_ranking_plot"),
              height = ceiling(length(unique(rss_df[[args$groupby]]))/2) * 4)

  #筛选出至少有一个值大于 rss_threshold 的所有行
  if ( is.null(args$threshold) ){
    rss_threshold = 0
  }else{
    rss_threshold = as.numeric(args$threshold)
  }
  rss_specific <- na.omit(rss_matrix)[apply(na.omit(rss_matrix),MARGIN = 1 ,FUN = function(x) any(x > rss_threshold)),]
  p = pheatmap(rss_specific,
    cellwidth = 18,
    cellheight = 18,
    color=colorRampPalette(c("#406AA8","white","#D91216"))(299),
    angle_col = 45,
    treeheight_col=20,
    treeheight_row=20, 
    border_color=NA)
  write.table(rss_df, file.path(args$outdir, "2.3.RSS_heatmap.xls"), sep = "\t", col.names = T, row.names =F, quote =F)
  save_plot(p,file.path(args$outdir, "2.3.RSS_heatmap"),	width = dim(rss_specific)[2]/2+3,height = dim(rss_specific)[1]*0.4)
}

#' CSI模块分析
perform_csi_analysis <- function(regulonAUC, args) {
  # CSI计算与聚类
  regulons_csi <- calculate_csi(regulonAUC)
  
  # 矩阵转换
  csi_matrix <- regulons_csi %>% tidyr::spread(regulon_2,CSI)
  future_rownames <- csi_matrix$regulon_1
  csi_matrix <- as.matrix(csi_matrix[,2:ncol(csi_matrix)])
  rownames(csi_matrix) <- future_rownames
  
  # 层次聚类
  hclust_result <- hclust(dist(csi_matrix, method = "euclidean"))
  clusters <- cutree(hclust_result, k = as.numeric(args$nclust))
  
  # 结果保存
  clusters_df <- data.frame(
    regulon = names(clusters),
    csi_module = clusters,
    row.names = NULL
  )
  write.table(
    clusters_df,
    file.path(args$outdir, "3.1.csi_module_annotation.xls"),
    sep = "\t", 
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
  
  # 热图注释数据准备
  row_annotation <- clusters_df
  rownames(row_annotation) <- NULL
  row_annotation$csi_module <- as.character(row_annotation$csi_module)
  rownames(row_annotation) <- row_annotation$regulon
  row_annotation <- row_annotation[, -which(names(row_annotation) == "regulon"),drop = F]

  # 可视化
  csi_plot <- plot_csi_modules(
    regulons_csi,
    args$groupby,
    nclust = as.numeric(args$nclust),
    row_anno = row_annotation,
    font_size_regulons = 8
  )
  
  save_plot(
    csi_plot,
    file.path(args$outdir, "3.2.regulons_csi_correlation_heatmap")
  )

  return(clusters_df)
}

#' CSI模块活动分析
analyze_csi_activity <- function(clusters_df, seurat_obj, regulonAUC, args) {
  # 数据准备
  cell_metadata = seurat_obj@meta.data
  cell_metadata[,"cell_type"] = cell_metadata[,args$groupby]

  # 模块活动计算
  activity_matrix <- calc_csi_module_activity(
    clusters_df,
    regulonAUC,
    cell_metadata
  )
  
  # 热图可视化
  activity_heatmap <- pheatmap(
    activity_matrix,
    show_colnames = TRUE,
    color = viridis::viridis(n = 10),
    cellwidth = 24,
    cellheight = 24,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean"
  )
  
  save_plot(
    activity_heatmap,
    file.path(args$outdir, "3.3.csi_module_activity_heatmap"),
    width = 8,
    height = 8
  )
}

#
main <- function() {
  args <- handle_parameters()
  source(args$utils_path)
  
  dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
  
  data_objects <- load_data(args)
  
  print('regulon_activity 分析')
  analyze_regulon_activity(data_objects$seurat, data_objects$regulonAUC, args)

  print('rss 分析')
  perform_rss_analysis(data_objects$seurat, data_objects$regulonAUC, args)

  print('CSI分析')
  clusters_df <- perform_csi_analysis(data_objects$regulonAUC, args)
  
  print('CSI activity 分析')
  analyze_csi_activity(
    clusters_df,
    data_objects$seurat, 
    data_objects$regulonAUC, 
    args
  )
  
  if (!file.exists(file.path(args$outdir, "Regulon调控子分析说明.docx"))) {
    file.copy(
      "/public/scRNA_works/Documents/Regulon调控子分析说明.docx",
      args$outdir
    )
  }
}

# 执行主函数
if (!interactive()) {
  main()
}
