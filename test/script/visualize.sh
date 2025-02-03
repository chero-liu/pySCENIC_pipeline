
#!/bin/bash
module load OESingleCell/3.0.d

Rscript /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/scripts/script/visualize.r \
  --input /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/test1/pyscenic_test_all_modules_True_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather/data/sce_SCENIC.loom \
  --rds_filepath /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/test1/CD4T_20240805/Manualanno/data_ob_v3.rds \
  --outdir /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/test1/pyscenic_test_all_modules_True_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather/result \
  --subnew_celltype all \
  --subsampleid all  \
  --subgroup all  \
  --subcluster all \
  --predicate  all \
  --binmethod aucell \
  --groupby group_new_celltype \
  --threshold 0.0 \
  --topGenes 4 \
  --extended FALSE \
  --nclust  4 \
  --utils_path /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/scripts/tools/utils.r

