
#!/bin/bash
module load OESingleCell/3.0.d
Rscript /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/scripts/script/save_seurat_counts.r \
    --input /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/test1/CD4T_20240805/Manualanno/data_ob_v3.rds \
    --outdir /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/test1/pyscenic_test_all_modules_True_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather/data \
    --subnew_celltype all \
    --subsampleid all  \
    --subgroup all  \
    --subcluster all \
    --predicate  all
