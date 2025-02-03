load_data <- function(args) {
  seurat_obj <- readSeurat(args$rds_filepath) %>% 
    subsetRDS(
      new_celltype = args$subnew_celltype,
      sampleid = args$subsampleid,
      group = args$subgroup,
      cluster = args$subcluster,
      predicate = args$predicate
    )
  if (is.null(levels(seurat_obj@meta.data[,args$groupby]))){
    seurat_obj@meta.data[,args$groupby] = reorderLevelsByFrequency(seurat_obj@meta.data[,args$groupby])
  }

  scenic_data <- open_loom(args$input)
  # 数据提取
  regulons_incidMat <- get_regulons(scenic_data, column.attr.name = 'Regulons')
  #   regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(scenic_data, column.attr.name = 'RegulonsAUC')
  # 数据对齐
  seurat_obj <- seurat_obj[, colnames(regulonAUC)]
  regulonAUC <- regulonAUC[, rownames(seurat_obj@meta.data)]
  regulonAUC <- regulonAUC[rowSums(regulonAUC) != 0, ]

  list(
    seurat = seurat_obj,
    regulonAUC = regulonAUC
  )
}

#RemoveOutlier 函数的作用是基于中位数和MAD（Median Absolute Deviation，中位数绝对偏差）来检测并移除异常值，支持 按批次处理数据，并提供不同的异常值筛选模式（上限、下限、或两者）
RemoveOutlier <- function(
  metric,
  nmads = 5,
  type = c("both", "lower", "higher"),
  log = FALSE,
  subset = NULL,
  batch = NULL,
  min_diff = NA
) {
    if (log) {
        metric <- log10(metric)
    }
    if (any(is.na(metric))) {
        warning("missing values ignored during outlier detection")
    }

    if (!is.null(batch)) {
        N <- length(metric)
        if (length(batch) != N) {
            stop("length of 'batch' must equal length of 'metric'")
        }

        # Coercing non-NULL subset into a logical vector.
        if (!is.null(subset)) {
            new.subset <- logical(N)
            names(new.subset) <- names(metric)
            new.subset[subset] <- TRUE
            subset <- new.subset
        }

        # Computing QC metrics for each batch.
        by.batch <- split(seq_len(N), batch)
        collected <- logical(N)
        all.threshold <- vector("list", length(by.batch))
        for (b in seq_along(by.batch)) {
            bdx <- by.batch[[b]]
            current <- Recall(metric[bdx], nmads = nmads, type = type, log = FALSE, subset = subset[bdx], batch = NULL, min_diff = min_diff)
            all.threshold[[b]] <- attr(x, "thresholds")
            collected[bdx] <- current
        }

        all.threshold <- do.call(cbind, all.threshold)
        colnames(all.threshold) <- names(by.batch)
        # return(.store_thresholds(collected, all.threshold, logged=log))
        if ( log ){ val <- 10^all.threshold }
        attr(collected, "thresholds") <- val
        return( collected )
    }
    # Computing median/MAD (possibly based on subset of the data).
    if (!is.null(subset)) {
        submetric <- metric[subset]
        if (length(submetric) == 0L) {
            warning("no observations remaining after subsetting")
        }
    } else {
        submetric <- metric
    }
    cur.med <- median(submetric, na.rm = TRUE)
    cur.mad <- mad(submetric, center = cur.med, na.rm = TRUE)

    diff.val <- max(min_diff, nmads * cur.mad, na.rm = TRUE)
    upper.limit <- cur.med + diff.val
    lower.limit <- cur.med - diff.val

    type <- match.arg(type)
    if (type == "lower") {
        upper.limit <- Inf
    } else if (type == "higher") {
        lower.limit <- -Inf
    }

    kx = metric < lower.limit | upper.limit < metric
    val = c(lower=lower.limit, higher=upper.limit)
    if ( log ){
      val <- 10^val
    }
    attr(kx, "thresholds") <- val
    return( kx )
}


process_regulon_activity <- function(cellInfo, regulonAUC_mat, groupby, nmads = 3, scale_by_group_count = TRUE) {
  group_sizes <- table(cellInfo[[groupby]])
  single_cell_groups <- names(group_sizes[group_sizes == 1])

  # Function to calculate mean after outlier removal
  calculate_mean_no_outliers <- function(values) {
    mean(values[!RemoveOutlier(values, nmads = nmads, type = "higher")])
  }

  if (length(single_cell_groups) == 0) {
    # No single-cell groups
    regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo[[groupby]]), function(cells) {
      apply(regulonAUC_mat[, cells, drop = FALSE], 1, calculate_mean_no_outliers)
    })
  } else {
    # Exclude single-cell groups temporarily
    cellInfo_filtered <- cellInfo[!cellInfo[[groupby]] %in% single_cell_groups, ]
    cellInfo_filtered[[groupby]] <- droplevels(cellInfo_filtered[[groupby]])

    regulonActivity_byCellType0 <- sapply(split(rownames(cellInfo_filtered), cellInfo_filtered[[groupby]]), function(cells) {
      apply(regulonAUC_mat[, cells, drop = FALSE], 1, calculate_mean_no_outliers)
    })

    # Add back single-cell groups
    single_cell_matrix <- regulonAUC_mat[, rownames(cellInfo[cellInfo[[groupby]] %in% single_cell_groups, ]), drop = FALSE]
    regulonActivity_byCellType <- cbind(regulonActivity_byCellType0, single_cell_matrix)
    colnames(regulonActivity_byCellType) <- c(colnames(regulonActivity_byCellType0), single_cell_groups)
  }

  # Scale and process the activity matrix
  scaled_data <- t(scale(t(regulonActivity_byCellType), center = TRUE, scale = ifelse(scale_by_group_count && length(unique(cellInfo[[groupby]])) == 2, FALSE, TRUE)))
  regulonActivity_byCellType_processed <- na.omit(scaled_data)

  result <- regulonActivity_byCellType_processed[which(rowSums(abs(regulonActivity_byCellType_processed)) > max(abs(regulonActivity_byCellType_processed))/4),]
  return(result)
}



RSSRanking <- function(
  rss_df,
   group.by,
   ggrepel_force = 1,
   ggrepel_point_padding = 0.2,
   top_genes = 4,
   plot_extended = FALSE
){
  require(ggrepel)
  require(cowplot)

  if(plot_extended == TRUE){
    rss_df <- rss_df %>%
      subset(grepl("extended",regulon))
  }else if(plot_extended == FALSE){
    rss_df <- rss_df %>%
      subset(!grepl("extended",regulon))
  }

  rss_df_sub <- rss_df %>% dplyr::group_by(.dots = group.by) %>%
    mutate("rank" = order(order(RSS, decreasing = TRUE)))

  rrs_ranking_plot <- ggplot(rss_df_sub,aes(rank,RSS,label = regulon)) +
    geom_point(color = "grey20",size = 2) +
    geom_point(data = subset(rss_df_sub,rank < top_genes),
               color = "red",size = 2) +
    geom_text_repel(data = subset(rss_df_sub,rank < top_genes),
                    force = ggrepel_force,point.padding = ggrepel_point_padding) +
    theme_bw() + theme(panel.grid =element_blank()) +
    labs(x = "Rank", y = "RSS", title = group.by) +
    facet_wrap(eval(expr(~!!ensym(group.by))), ncol = 2, scales = "free_y" )
  return(rrs_ranking_plot)
}


calculate_csi <- function(
  regulonAUC = regulonAUC,
  calc_extended = FALSE,
  verbose = FALSE
){
  if (calc_extended == FALSE){
    regulonAUC <- subset(regulonAUC,!grepl("extended",rownames(regulonAUC)))
  }

  pearson_cor <- cor(t(regulonAUC))
  pearson_cor_df <- as.data.frame(pearson_cor)
  pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
  pearson_cor_long <- pearson_cor_df %>%
    tidyr::gather(regulon_2,pcc,-regulon_1) %>%
    dplyr::mutate("regulon_pair" = paste(regulon_1,regulon_2,sep="_"))

  regulon_names <- unique(colnames(pearson_cor))
  num_of_calculations <- length(regulon_names)*length(regulon_names)
  csi_regulons <- data.frame(matrix(nrow=num_of_calculations,ncol = 3))
  colnames(csi_regulons) <- c("regulon_1", "regulon_2", "CSI")
  num_regulons <- length(regulon_names)

  f <- 0
  for(reg in regulon_names){
    for(reg2 in regulon_names){
      f <- f + 1
      # fraction_lower <- calc_csi(reg,reg2,pearson_cor)
      test_cor <- pearson_cor[reg,reg2]
      total_n <- ncol(pearson_cor)
      pearson_cor_sub <- subset(pearson_cor,rownames(pearson_cor) == reg | rownames(pearson_cor) == reg2)

      # sums <- apply(pearson_cor_sub,MARGIN = 2, FUN = compare_pcc, pcc = test_cor)
      sums <- apply(pearson_cor_sub, 2,function(m) ifelse( length(m[m>test_cor]) == length(m), 0, length(m)) )
      fraction_lower <- length(sums[sums == nrow(pearson_cor_sub)]) / total_n
      csi_regulons[f,] <- c(reg,reg2,fraction_lower)
    }
  }
  csi_regulons$CSI <- as.numeric(csi_regulons$CSI)
  return(csi_regulons)
}

plot_csi_modules <- function(
  csi_df,
  groupby,
  nclust = 4,
  row_anno = NULL,
  font_size_regulons = 6
){
  ## subset csi data frame based on threshold
  csi_test_mat <- csi_df %>% tidyr::spread(regulon_2,CSI)

  future_rownames <- csi_test_mat$regulon_1
  csi_test_mat <- as.matrix(csi_test_mat[,2:ncol(csi_test_mat)])
  rownames(csi_test_mat) <- future_rownames

  color_map = c("#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b",
    "#666666","#1b9e77","#d95f02","#7570b3","#d01b2a","#43acde",
    "#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff",
    "#ecb9e5","#813139","#743fd2","#434b7e","#e6908e","#214a00",
    "#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7",
    "#006874","#d2ad2c","#b5d7a5","#9e8442","#4e1737","#e482a7",
    "#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99",
    "#3fdacb","#bf5b17")[1:nclust]
    
  names(color_map) <- 1:nclust
  color_use = list()

  color_use[[groupby]] <- color_map  ##different with local script 
  
  pheatmap::pheatmap(csi_test_mat,
           show_colnames = TRUE,
           #border_color = NA,
           color = viridis::viridis(n = 10),
           fontsize_row = font_size_regulons,
           fontsize_col = font_size_regulons,
           angle_col = 90,
           cutree_cols = nclust,
           cutree_rows = nclust,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           annotation_row = row_anno,
           annotation_colors = color_use,
           treeheight_row = 20,
           treeheight_col = 20,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           width = 2000,
           height = 3200)
}

calc_csi_module_activity <- function(
  clusters_df,
  regulonAUC_mat,
  metadata
){
  cell_types<- unique(metadata$cell_type)
  regulons <- unique(clusters_df$regulon)

  regulonAUC_mat <- regulonAUC_mat[as.character(regulons),]

  csi_activity_matrix_list <- list()
  csi_cluster_activity <- data.frame("csi_module" = c(),
                                    "mean_activity" = c(),
                                    "cell_type" = c())

  cell_type_counter <- 0
  groupby_1 =c()   # select group with only one cell
  for(i in names(table(metadata$cell_type)) ){
	if( table(metadata$cell_type)[i]==1 ){
		groupby_1= i
	}
  }
  cell_types <- setdiff(cell_types, groupby_1)
  
  regulon_counter <-
    for(ct in cell_types) {
      cell_type_counter <- cell_type_counter + 1
      cell_type_aucs <- rowMeans(regulonAUC_mat[,rownames(subset(metadata,cell_type == ct))])
	  
      cell_type_aucs_df <- data.frame("regulon" = names(cell_type_aucs),
                                      "activtiy"= cell_type_aucs,
                                      "cell_type" = ct)
      csi_activity_matrix_list[[ct]] <- cell_type_aucs_df
    }
  if(!is.null(groupby_1)){
    regulonAUC_mat_1 = regulonAUC_mat[,rownames(subset(metadata,cell_type == groupby_1))]
    csi_activity_matrix_list_1=data.frame("regulon" = names(regulonAUC_mat_1),
                                      "activtiy"= as.vector(regulonAUC_mat_1),
                                      "cell_type" = groupby_1)
    csi_activity_matrix_list[[groupby_1]] <- csi_activity_matrix_list_1
  }
  for(ct in names(csi_activity_matrix_list)){
    for(cluster in unique(clusters_df$csi_module)){
      csi_regulon <- subset(clusters_df,csi_module == cluster)
      csi_regulon_activtiy <- subset(csi_activity_matrix_list[[ct]],regulon %in% csi_regulon$regulon)
      csi_activtiy_mean <- mean(csi_regulon_activtiy$activtiy)
      this_cluster_ct_activity <- data.frame("csi_module" = cluster,
                                             "mean_activity" = csi_activtiy_mean,
                                             "cell_type" = ct)
      csi_cluster_activity <- rbind(csi_cluster_activity,this_cluster_ct_activity)
    }
  }

  csi_cluster_activity[is.na(csi_cluster_activity)] <- 0

  csi_cluster_activity_wide <- csi_cluster_activity %>%
    spread(cell_type,mean_activity)

  rownames(csi_cluster_activity_wide) <- csi_cluster_activity_wide$csi_module
  csi_cluster_activity_wide <- as.matrix(csi_cluster_activity_wide[2:ncol(csi_cluster_activity_wide)])

  return(csi_cluster_activity_wide)
}